

.onLoad <- function(libname, pkgname) {
  required_packages <- c("moments", "gramEvol", "dplyr", "odbc", "DBI")
  options(repos = c(CRAN = "https://cran.rstudio.com/"))
  # Function to check and install missing packages
  check_install_packages <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }

  # Check and load each package
  invisible(lapply(required_packages, check_install_packages))
}


#' Generate Distribution Metrics
#'
#' This function generates random distributions of a specified type and computes
#' key statistical metrics such as mean, standard deviation, skewness, and entropy.
#'
#' @param num_distributions An integer specifying the number of distributions to generate (default: 1000).
#' @param samples_per_dist An integer specifying the number of samples per distribution (default: 500).
#' @param dist_type A string specifying the type of distribution to generate.
#'        Options include `"log-normal"`, `"normal"`, `"uniform"`, and `"exponential"`.
#'
#' @return A data frame containing computed metrics for each generated distribution.
#'         Columns include:
#'         - `distNbr`: Distribution number (1 to `num_distributions`).
#'         - `distType`: Type of distribution generated.
#'         - `mu`: Mean parameter used for generating the distribution.
#'         - `sigma`: Standard deviation parameter used for generating the distribution.
#'         - `n`: Number of samples in each distribution.
#'         - `xbar`: Sample mean (log-transformed where applicable).
#'         - `std`: Sample standard deviation.
#'         - `skew`: Sample skewness.
#'         - `entropy`: Entropy (for log-normal and normal distributions).
#'
#' @examples
#' # Generate 100 distributions with 500 samples each (default settings)
#' metrics_df <- generate_distribution_metrics()
#'
#' # Generate 500 normal distributions with 1000 samples per distribution
#' normal_metrics <- generate_distribution_metrics(num_distributions = 500,
#'                                                 samples_per_dist = 1000,
#'                                                 dist_type = "normal")
#'
#' # Generate 200 exponential distributions with 250 samples per distribution
#' exp_metrics <- generate_distribution_metrics(num_distributions = 200,
#'                                              samples_per_dist = 250,
#'                                              dist_type = "exponential")
#'
#' @export
generate_distribution_metrics <- function(num_distributions = 1000,
                                          samples_per_dist = 500,
                                          dist_type = c("log-normal", "normal", "uniform", "exponential")) {


  # Ensure valid distribution type
  dist_type <- match.arg(dist_type)

  # Function to calculate entropy (for log-normal and normal distributions)
  calculate_entropy <- function(std_dev, mean_val) {
    log(std_dev * exp(mean_val + 1 / 2) * sqrt(2 * pi))
  }

  # Preallocate vectors for performance
  distNbr <- seq_len(num_distributions)
  mu <- log(runif(num_distributions, min = 1, max = 100))
  sigma <- log(runif(num_distributions, min = 1, max = 10))

  # Predefine storage for results
  results_df <- data.frame(
    distNbr = integer(num_distributions),
    distType = character(num_distributions),
    mu = numeric(num_distributions),
    sigma = numeric(num_distributions),
    n = integer(num_distributions),
    xbar = numeric(num_distributions),
    std = numeric(num_distributions),
    skew = numeric(num_distributions),
    entropy = numeric(num_distributions),
    stringsAsFactors = FALSE
  )

  # Vectorized generation of distributions
  set.seed(123456789)

  for (i in seq_len(num_distributions)) {
    # Generate the sample based on selected distribution type
    sample_data <- switch(
      dist_type,
      "log-normal" = rlnorm(samples_per_dist, meanlog = mu[i], sdlog = sigma[i]),
      "normal" = rnorm(samples_per_dist, mean = exp(mu[i]), sd = exp(sigma[i])),
      "uniform" = runif(samples_per_dist, min = exp(mu[i] - sigma[i]), max = exp(mu[i] + sigma[i])),
      "exponential" = rexp(samples_per_dist, rate = 1 / exp(mu[i]))
    )

    # Compute log-transformed statistics where applicable
    log_sample_data <- if (dist_type == "log-normal") log(sample_data) else sample_data
    sample_mean <- mean(log_sample_data)
    sample_std <- sd(log_sample_data)
    sample_skew <- skewness(log_sample_data)
    sample_entropy <- if (dist_type %in% c("log-normal", "normal")) calculate_entropy(sample_std, sample_mean) else NA

    # Assign values to the preallocated dataframe
    results_df[i, ] <- list(
      distNbr = i,
      distType = dist_type,
      mu = mu[i],
      sigma = sigma[i],
      n = samples_per_dist,
      xbar = sample_mean,
      std = sample_std,
      skew = sample_skew,
      entropy = sample_entropy
    )
  }

  return(results_df)
}


#' Perform Symbolic Regression using Grammatical Evolution
#'
#' This function applies symbolic regression to find an optimal mathematical expression
#' that approximates entropy based on input variables using grammatical evolution.
#'
#' @param data A data frame containing at least the following columns:
#'        `"mu"`, `"sigma"`, `"xbar"`, and `"entropy"`.
#' @param termination_cost A numeric value specifying the termination cost for the evolution process (default: 0.05).
#' @param optimizer A string specifying the optimization algorithm to use (default: `"es"`).
#' @param iterations An integer specifying the number of iterations for the evolution process (default: 3e6).
#' @param seed An integer used for random number generation to ensure reproducibility (default: 2).
#'
#' @return A list containing:
#'         - `best_expression`: The best mathematical expression found.
#'         - `runtime`: Execution time of the symbolic regression process.
#'         - `iterations`: The number of iterations performed.
#'         - `optimizer`: The optimization method used.
#'
#' @examples
#' # Generate sample data
#' sample_data <- data.frame(
#'   mu = runif(100, 1, 10),
#'   sigma = runif(100, 1, 5),
#'   xbar = rnorm(100, 5, 2),
#'   entropy = rexp(100, rate = 0.5)
#' )
#'
#' # Perform symbolic regression
#' result <- symbolic_regression(sample_data)
#'
#' # View best expression found
#' print(result$best_expression)
#'
#' @export
symbolic_regression <- function(data,
                                termination_cost = 0.05,
                                optimizer = "es",
                                iterations = 3e6,
                                seed = 2) {
  # Validate input data
  if (!all(c("mu", "sigma", "xbar", "entropy") %in% colnames(data))) {
    stop("Input data must contain columns: 'mu', 'sigma', 'xbar', 'entropy'")
  }

  # Extract relevant columns
  o_mu <- data$mu
  o_sigma <- data$sigma
  o_xbar <- data$xbar
  known_entropy <- data$entropy

  # Define the fitness function
  SymRegFitFunc <- function(expr) {
    result <- eval(expr)
    if (any(is.nan(result))) return(Inf)
    return(mean(log(1 + abs(known_entropy - result))))
  }

  # Define grammar rules
  ruleDef <- list(
    expr = grule(op(expr, expr), func(expr), var),
    func = grule(sin, cos, log, sqrt),
    op = grule("+", "-", "*", "/", "^"),
    var = grule(o_mu, o_sigma, o_xbar, n),
    n = grule(1, 2, 3)
  )

  # Create the grammar definition
  grammarDef <- CreateGrammar(ruleDef)

  # Set seed for reproducibility
  set.seed(seed)

  # Start runtime tracking
  start_time <- Sys.time()

  # Run Grammatical Evolution
  gramEvolution <- GrammaticalEvolution(
    grammarDef, SymRegFitFunc,
    terminationCost = termination_cost,
    optimizer = optimizer,
    iterations = iterations,
    monitorFunc = print
  )

  # End runtime tracking
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "secs")

  # Extract best expression
  best_expression <- gramEvolution$best$expression

  # Return results
  return(list(
    best_expression = best_expression,
    runtime = runtime,
    iterations = iterations,
    optimizer = optimizer
  ))
}


# # GIT CLONE HELPER:
# cd path\to\folder  # Navigate to the target directory
# git clone https://github.com/AaronRoseDA/srlib.git  # Clone the repository
# cd srlib  # Move into the cloned repo directory
#
#
# # GIT PUSH HELPER:
# cd srlib
# git remote set-url origin https://github.com/AaronRoseDA/srlib.git  # Set remote origin (only needed if not set)
# git pull origin main  # Ensure your local branch is up-to-date
# git add .  # Stage all changes
# git commit -m "Added symbolic_regression function"  # Commit changes
# git push origin main  # Push to the main branch
