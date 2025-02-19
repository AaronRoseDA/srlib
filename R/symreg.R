

.onLoad <- function(libname, pkgname) {
  required_packages <- c("moments", "gramEvol", "dplyr", "odbc", "DBI")

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

hello <- function(digit) {
  # Use R's built-in pi constant
  pi_str <- as.character(pi)
  pi_digits <- gsub("\\.", "", pi_str)  # Remove the decimal point

  # Ensure digit is within valid range
  if (digit < 1 || digit > nchar(pi_digits)) {
    return("Digit out of range. Try a smaller number.")
  }

  # Extract the requested digit
  selected_digit <- substr(pi_digits, digit, digit)

  # Return the result
  paste("The", digit, "digit of pi is:", selected_digit)
}


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



another_test_func <- function(n = 20){
  return(rnorm(n = n,mean = 0,sd = n/exp(1)))
}


# # GIT CLONE HELPER:
# cd path\to\folder  # Navigate to the target directory
# git clone https://github.com/AaronRoseDA/srlib.git  # Clone the repository
# cd srlib  # Move into the cloned repo directory
#
#
# # GIT PUSH HELPER:
# git remote set-url origin https://github.com/AaronRoseDA/srlib.git  # Set remote origin (only needed if not set)
# git pull origin main  # Ensure your local branch is up-to-date
# git add .  # Stage all changes
# git commit -m "Added symbolic_regression function"  # Commit changes
# git push origin main  # Push to the main branch
