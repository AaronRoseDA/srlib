

.onLoad <- function(libname, pkgname) {
  required_packages <- c("moments", "gramEvol", "dplyr", "odbc", "DBI", "e1071")
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

#####
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
#####
generate_distribution_metrics <- function(num_distributions = 1000,
                                          samples_per_dist = 500,
                                          dist_type, a = F, b = F) {

  generateNormalDistributions <- function(num_distributions = 1000,
                                          samples_per_dist = 500,
                                          mu_sd = 1, sigma_max = 10) {
    # Create an empty results data frame
    results_df <- data.frame(
      distNbr  = integer(num_distributions),
      distType = character(num_distributions),
      n        = integer(num_distributions),
      mu       = numeric(num_distributions),
      sigma    = numeric(num_distributions),
      skew     = numeric(num_distributions),
      entropy  = numeric(num_distributions),
      stringsAsFactors = FALSE
    )


    for (i in seq_len(num_distributions)) {

      # -- CHOOSE TRUE PARAMETERS (mu_i, sigma_i) --
      mu_i    <- rnorm(1, mean = 0, sd = mu_sd)
      sigma_i <- runif(1, min = 0.5, max = sigma_max)

      # -- GENERATE SAMPLES --
      x <- rnorm(n = samples_per_dist, mean = mu_i, sd = sigma_i)

      # -- CALCULATE SAMPLE STATISTICS --
      xbar_i <- mean(x)
      std_i  <- sd(x)
      skew_i <- skewness(x, type = 2)

      # -- ENTROPY FOR NORMAL DISTRIBUTION --
      entropy_i <- 0.5 * log(2 * pi * exp(1) * std_i^2)

      # -- FILL THE DATA FRAME --
      results_df$distNbr[i]  <- i
      results_df$distType[i] <- "normal"
      results_df$mu[i]       <- xbar_i
      results_df$sigma[i]    <- std_i
      results_df$n[i]        <- samples_per_dist
      results_df$skew[i]     <- skew_i
      results_df$entropy[i]  <- entropy_i
    }

    # Return the completed data frame
    return(results_df)
  }

  generateLognormalDistributions <- function(num_distributions = 1000,
                                             samples_per_dist = 500) {
    # Create an empty results data frame
    results_df <- data.frame(
      distNbr  = integer(num_distributions),
      distType = character(num_distributions),
      n        = integer(num_distributions),
      mu       = numeric(num_distributions),
      sigma    = numeric(num_distributions),
      skew     = numeric(num_distributions),
      entropy  = numeric(num_distributions),
      stringsAsFactors = FALSE
    )

    for (i in seq_len(num_distributions)) {
      # -- Choose mu and sigma --
      muLog_i    <- rnorm(1, mean = 0, sd = 1)
      sigmaLog_i <- runif(1, min = 0.5, max = 1)

      # samples_per_dist <- 10000000
      # -- Generate samples from the log-normal distribution --
      x <- rlnorm(n = samples_per_dist, meanlog = muLog_i, sdlog = sigmaLog_i)
      # hist(x,breaks = 1000)
      # hist(log(x),breaks = 1000)

      # -- find true mu and sigma per https://medium.com/towards-data-science/log-normal-distribution-a-simple-explanation-7605864fb67c#:~:text=The%20probability%20density%20function%20for,deviation%20from%20a%20normal%20distribution%20.--
      muLog_i    <- mean(log(x))
      sigmaLog_i <- sd(log(x))

      # -- Compute sample statistics --
      skew_i <- (exp(sigmaLog_i^2) + 2) * sqrt(exp(sigmaLog_i^2) - 1)

      # -- Calculate entropy  per https://en.wikipedia.org/wiki/Log-normal_distribution --
      entropy_i <- log2(sqrt(2 * pi) * sigmaLog_i * exp(muLog_i + 0.5))

      # -- Fill the data frame --
      results_df$distNbr[i]  <- i
      results_df$distType[i] <- "log-normal"
      results_df$mu[i]       <- muLog_i    # storing the "log" mean
      results_df$sigma[i]    <- sigmaLog_i # storing the "log" SD
      results_df$n[i]        <- samples_per_dist
      results_df$skew[i]     <- skew_i
      results_df$entropy[i]  <- entropy_i
    }

    return(results_df)
  }

  generateExponentialDistributions <- function(num_distributions = 1000,
                                               samples_per_dist = 500) {
    # Create an empty results data frame
    results_df <- data.frame(
      distNbr  = integer(num_distributions),
      distType = character(num_distributions),
      n        = integer(num_distributions),
      mu       = numeric(num_distributions),  # will store 1/lambda (true mean)
      sigma    = numeric(num_distributions),  # will store 1/lambda (true SD)
      skew     = numeric(num_distributions),
      lambda   = numeric(num_distributions),
      entropy  = numeric(num_distributions),
      stringsAsFactors = FALSE
    )

    for (i in seq_len(num_distributions)) {
      # -- Choose a rate parameter lambda randomly (customize as needed) --

      # samples_per_dist <- 1000000

      lambda_i <- runif(1, min = 0, max = 1)

      # -- Generate samples from the exponential distribution --
      x <- rexp(n = samples_per_dist, rate = lambda_i)

      # hist(x, breaks = 1000)
      # -- Compute sample statistics --
      mu <- 1/lambda_i
      skew_i <- 2

      # -- Calculate the differential entropy -- https://en.wikipedia.org/wiki/Exponential_distribution
      entropy_i <- 1 - log(lambda_i)

      # -- Fill the data frame --
      results_df$distNbr[i]  <- i
      results_df$distType[i] <- "exponential"
      results_df$mu[i]       <- mu
      results_df$sigma[i]    <- mu
      results_df$n[i]        <- samples_per_dist
      results_df$skew[i]     <- skew_i
      results_df$lambda[i]   <- lambda_i
      results_df$entropy[i]  <- entropy_i
    }

    return(results_df)
  }

  generateUniformDistributions <- function(num_distributions = 1000,
                                           samples_per_dist = 500, a = F, b = F) {
    # Create an empty results data frame
    results_df <- data.frame(
      distNbr  = integer(num_distributions),
      distType = character(num_distributions),
      n        = integer(num_distributions),
      mu       = numeric(num_distributions),  # (a + b)/2
      sigma    = numeric(num_distributions),  # sqrt((b - a)^2 / 12)
      skew     = numeric(num_distributions),
      a        = numeric(num_distributions),
      b        = numeric(num_distributions),
      entropy  = numeric(num_distributions),
      stringsAsFactors = FALSE
    )

    for (i in seq_len(num_distributions)) {
      # samples_per_dist <- 100000
      # -- Choose parameters a and b where a < b --
      if (is.numeric(a) & is.numeric(b) & a < b){
        a_i <- a
        b_i <- b
      } else{
        a_i <- runif(1, min = 0, max = 5)
        b_i <- runif(1, min = a_i, max = 10)  # ensures b_i > a_i
      }

      # -- Generate samples from the uniform distribution --
      x <- runif(n = samples_per_dist, min = a_i, max = b_i)

      # -- Compute sample statistics --
      xbar_i <- (a_i + b_i) / 2
      std_i  <- (b_i - a_i) / sqrt(12)
      skew_i <- 0

      # -- Compute the theoretical mean and std dev of the uniform distribution --
      mu_i    <- (a_i + b_i) / 2
      sigma_i <- (b_i - a_i) / sqrt(12)

      # -- Calculate the entropy https://en.wikipedia.org/wiki/Continuous_uniform_distribution
      entropy_i <- log(b_i - a_i)

      # -- Fill the data frame --
      results_df$distNbr[i]  <- i
      results_df$distType[i] <- "uniform"
      results_df$mu[i]       <- mu_i
      results_df$sigma[i]    <- sigma_i
      results_df$n[i]        <- samples_per_dist
      results_df$skew[i]     <- skew_i
      results_df$entropy[i]  <- entropy_i
      results_df$a[i]        <- a_i
      results_df$b[i]        <- b_i
    }

    return(results_df)
  }

  if (dist_type == "log-normal") {
    results_df <- generateLognormalDistributions(num_distributions,samples_per_dist)
  }
  if (dist_type == "normal") {
    results_df <- generateNormalDistributions(num_distributions,samples_per_dist)
  }
  if (dist_type == "uniform") {
    results_df <- generateUniformDistributions(num_distributions,samples_per_dist)
  }
  if (dist_type == "exponential") {
    results_df <- generateExponentialDistributions(num_distributions,samples_per_dist)
  }

  return(results_df)
}





{
#' Perform Symbolic Regression using Grammatical Evolution TEST TEST TEST
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
}
symbolic_regression <- function(data,
                                termination_cost = 0.05,
                                optimizer = "es",
                                iterations = 3e6,
                                seed = 2) {
  # dataN = generate_distribution_metrics(dist_type = "normal",num_distributions = 10,samples_per_dist = 10000000)
  # dataE = generate_distribution_metrics(dist_type = "exponential",num_distributions = 10,samples_per_dist = 10000000)
  # dataU = generate_distribution_metrics(dist_type = "uniform",num_distributions = 10,samples_per_dist = 10000000)

  # termination_cost = 0.05
  # optimizer = "ga"
  # iterations = 3e6
  # seed = 2

  # Set seed for reproducibility
  set.seed(seed)

  SymRegFitFunc <- function(expr) {
    result <- eval(expr)
    if (any(is.nan(result))) return(Inf)
    return(mean(log(1 + abs(known_entropy - result))))
  }

  # Validate input data
  if (!all(c("mu", "sigma", "entropy") %in% colnames(data))) {
    stop("Input data must contain columns: 'mu', 'sigma', 'entropy'")
  }

  ### REMOVE
  data <- dataU
  # Extract relevant columns
  if(data$distType[1] == "uniform"){
    a <- data$a
    b <- data$b
    known_entropy <- data$entropy

    # Define grammar rules
    ruleDef <- list(
      expr = grule(op(expr, expr), func(expr), var),
      func = grule(log, sqrt),
      op = grule("+", "-", "*", "/", "^"),
      var = grule(a, b),
      n = grule(1, 2, 3)
    )
    # Create the grammar definition
    grammarDef <- CreateGrammar(ruleDef)

  }


  ### REMOVE
  data <- dataE
  if(data$distType[1] == "exponential"){
    lambda <- data$lambda
    known_entropy <- data$entropy

    # Define grammar rules
    ruleDef <- list(
      expr = grule(op(expr, expr), func(expr), var),
      func = grule(log, sqrt),
      op = grule("+", "-", "*", "/", "^"),
      var = grule(lambda, 1,2),
      n = grule(1, 2, 3)
    )
    # Create the grammar definition
    grammarDef <- CreateGrammar(ruleDef)

  }


  ### REMOVE
  data <- dataN
  if (data$distType[1] == "normal") {
    sigma <- data$sigma
    known_entropy <- data$entropy

    # Define grammar rules
    ruleDef <- list(
      expr = grule(op(expr, expr), func(expr), var),
      func = grule(log, sqrt, exp),
      op   = grule("+", "-", "*", "/", "^"),
      var  = grule(sigma, 1, 2, pi)
    )

    # Create the grammar definition
    grammarDef <- CreateGrammar(ruleDef)
  }

  # termination_cost = 0.01
  # optimizer = "ga"
  # iterations = 3e6
  # seed = 2
  #
  # # Set seed for reproducibility
  # set.seed(seed)

  # ------------------------------------------------------------------
  # Cost (fitness) function
  SymRegFitFunc <- function(expr) {
    result <- eval(expr)
    if (any(is.nan(result))) {
      return(Inf)
    }
    mean(log(1 + abs(known_entropy - result)))
  }

  # ------------------------------------------------------------------
  # Parallel lapply function
  # myplapply <- function(X, FUN,...) {
  #   cl <- parallel::makeCluster(parallel::detectCores() - 1)
  #   on.exit(parallel::stopCluster(cl))
  #
  #   # Export only what the cost function needs on each worker
  #   parallel::clusterExport(
  #     cl,
  #     c("SymRegFitFunc", "sigma", "known_entropy"),
  #     envir = environment()
  #   )
  #
  #   # Load libraries on each worker
  #   parallel::clusterEvalQ(cl, {
  #     library(gramEvol)
  #     library(stats)
  #   })
  #
  #   # Apply in parallel
  #   parallel::parLapply(cl, X, FUN, ...)
  # }

  # ------------------------------------------------------------------

  # results <- list()
  # currentBest <- 1
  # gramEvolution <- GrammaticalEvolution(
  #   grammarDef      = grammarDef,
  #   evalFunc        = SymRegFitFunc,
  #   terminationCost = termination_cost,
  #   optimizer       = optimizer,
  #   iterations      = iterations,
  #   mutationChance  = 0.1,
  #
  #   monitorFunc = function(x) {
  #     if (x$best$cost < currentBest) {
  #       print(x)
  #       currentBest <<- x$best$cost
  #       results[[length(results) + 1]] <<- list(x)
  #       print(currentBest)
  #     }
  #   },
  #   plapply = myplapply
  # )




  results <- list()

  # Start runtime tracking
  start_time <- Sys.time()
  currentBest <- 1
  # Run Grammatical Evolution
  gramEvolution <- GrammaticalEvolution(
    grammarDef,
    SymRegFitFunc,
    terminationCost = termination_cost,
    optimizer = optimizer,
    iterations = iterations,mutationChance = .1,
    monitorFunc = function(x) {
      if (x$best$cost < currentBest) {
        print(x)
        currentBest <<- x$best$cost
        results[[length(results)+1]] <<- list(x)
        print(currentBest)
      }
    }
  )

  # resultsU <- results
   # resultsE <- results

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
