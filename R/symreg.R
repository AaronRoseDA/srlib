

.onLoad <- function(libname, pkgname) {
  required_packages <- c("moments", "gramEvol", "dplyr", "odbc", "DBI", "e1071", "parallel")
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

##### ---- ROXYGEN ---- #####
#' Generate Distribution Metrics
#'
#' This function generates synthetic data for different probability distributions and calculates
#' key metrics such as mean (`mu`), standard deviation (`sigma`), skewness, and entropy.
#'
#' @param num_distributions Integer. The number of distributions to generate (default: `1000`).
#' @param samples_per_dist Integer. The number of samples per generated distribution (default: `500`).
#' @param dist_type Character. Specifies the type of distribution to generate. Options include:
#'   - `"normal"`: Generates normal distributions.
#'   - `"log-normal"`: Generates log-normal distributions.
#'   - `"uniform"`: Generates uniform distributions.
#'   - `"exponential"`: Generates exponential distributions.
#' @param a Numeric or `FALSE`. The lower bound for uniform distributions. If `FALSE`, a random value is used.
#' @param b Numeric or `FALSE`. The upper bound for uniform distributions. If `FALSE`, a random value is used.
#'
#' @return A `data.frame` containing the following columns:
#'   \describe{
#'     \item{`distNbr`}{The index number of the generated distribution.}
#'     \item{`distType`}{The type of distribution generated (`"normal"`, `"log-normal"`, `"uniform"`, or `"exponential"`).}
#'     \item{`n`}{The number of samples in each distribution.}
#'     \item{`mu`}{The calculated mean of the generated distribution.}
#'     \item{`sigma`}{The standard deviation of the generated distribution.}
#'     \item{`skew`}{The skewness of the generated distribution.}
#'     \item{`entropy`}{The entropy of the distribution, calculated using theoretical formulas.}
#'     \item{`lambda`}{(For exponential distributions) The rate parameter.}
#'     \item{`a`}{(For uniform distributions) The lower bound.}
#'     \item{`b`}{(For uniform distributions) The upper bound.}
#'   }
#'
#' @details
#' This function internally calls specialized sub-functions for each distribution type:
#' \itemize{
#'   \item `generateNormalDistributions()`
#'   \item `generateLognormalDistributions()`
#'   \item `generateUniformDistributions()`
#'   \item `generateExponentialDistributions()`
#' }
#' These sub-functions generate random samples and compute distribution-specific statistics.
#'
#' @examples
#' \dontrun{
#'   # Generate metrics for 500 normal distributions with 1000 samples each
#'   normal_data <- generate_distribution_metrics(500, 1000, "normal")
#'
#'   # Generate 200 uniform distributions with predefined bounds
#'   uniform_data <- generate_distribution_metrics(200, 500, "uniform", a = 0, b = 10)
#' }
#'
#' @export
#####
generate_distribution_metrics <- function(num_distributions = 1000,
                                          samples_per_dist = 500,
                                          dist_type, a = F, b = F) {

  generateNormalDistributions <- function(num_distributions = 1000,
                                          samples_per_dist = 500,
                                          mu_max = 100000, sigma_max = 100000) {
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
      mu_i    <- runif(1, min = 0, max = mu_max)
      sigma_i <- runif(1, min = 0.5, max = sigma_max)

      # -- GENERATE SAMPLES --
      x <- rnorm(n = samples_per_dist, mean = mu_i, sd = sigma_i)

      ## KP CONT
      ########################
      # Perform KDE with default bandwidth
      # uses a gaussian kernel; no option to change
      # kde_result <- kde(x)

      # Extract bandwidth (assuming Gaussian kernel)
      # h <- kde_result$h ### SAVE AS kde_band
      #
      # # KDE-estimated mean (unbiased, no correction needed)
      # kde_mean <- sum(kde_result$eval.points * kde_result$estimate) * diff(kde_result$eval.points[1:2])
      #
      # # KDE-estimated variance (biased upward due to kernel smoothing)
      # kde_variance_raw <- sum((kde_result$eval.points - kde_mean)^2 * kde_result$estimate) * diff(kde_result$eval.points[1:2])
      #
      # # Apply theoretical bias-correction for Gaussian kernel
      # kde_variance_corrected <- kde_variance_raw - h^2
      #
      # # Calculate KDE-based entropy estimate
      # density_vals <- predict(kde_result, x = data)
      # entropy_kde <- -mean(log(density_vals))

      # Display all results clearly
      # entropy_theoretical <- uni_normal_entropy(theoretical_variance)
      # data_entropy_plugin <- uni_normal_entropy(data_variance)
      # kde_entropy_plugin <- uni_normal_entropy(kde_variance_corrected)


      ##########################

      # -- CALCULATE SAMPLE STATISTICS --
      xbar_i <- mean(x)
      std_i  <- sd(x)
      skew_i <- e1071::skewness(x, type = 2)

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


norm_df <-
  generate_distribution_metrics(
    num_distributions = 10,
    samples_per_dist = 100000,
    dist_type = "normal"
  )


##### ---- ROXYGEN ---- #####
#' Evaluate the Fitness of a Symbolic Regression Expression
#'
#' This function evaluates the fitness of a given expression in the context of symbolic regression.
#' The function calculates the mean log error between the predicted values (evaluated from the expression)
#' and the known entropy values in the dataset.
#'
#' @param expr A character string representing a mathematical expression to be evaluated.
#'             The expression should reference variables present in the dataset (`data`).
#'
#' @return A numeric value representing the fitness (cost) of the expression.
#'         Lower values indicate better fitness, while `Inf` is returned if the result contains `NaN` values.
#'
#' @details
#' The function first evaluates the given expression using `eval(parse(text = expr))`.
#' If the evaluation results in any `NaN` values, the function returns `Inf`, indicating an invalid solution.
#' Otherwise, the function computes the mean log error using:
#' \deqn{\text{cost} = \text{mean}(\log(1 + | \text{data\$entropy} - \text{result} |))}
#'
#' @examples
#' \dontrun{
#'   data <- generate_distribution_metrics(10, 1000000, "normal")
#'
#'   expr <- "log(data$b - log(exp(data$a)))"
#'   cost <- SymRegFitFunc(expr)
#'   print(cost)
#' }
#'
#' @export
#####
SymRegFitFunc <- function(expr) {
  # print(expr)
  result <-unlist(eval(parse(text = expr)))
  if (any(is.nan(result))) {
    return(Inf)
  }
  mean(log(1 + abs(data$entropy - result)))
}

SSE_eval_func <- function(expr, data) {
  # print(expr)
  result <-unlist(
    eval(parse(text = expr))
  )
  if (any(is.nan(result))) {
    return(Inf)
  }
  sum(abs(data$entropy - result)^2)
}



##### --- ROXYGEN ---- #####
#' Perform Symbolic Regression using Grammatical Evolution
#'
#' This function performs symbolic regression using Grammatical Evolution (GE).
#' It allows the user to define a grammar, an evaluation function, and various
#' optimization settings to evolve mathematical expressions that fit a given dataset.
#'
#' @param grammarDef A list containing the grammar definition, variable names, and additional metadata.
#'                   This should include:
#'                   - `name`: A character vector of variable names.
#'                   - `variables`: A vector of variables included in the regression.
#'                   - `grammarDef`: The formal grammar object used for GE.
#' @param evalFunc A function that evaluates the cost of an expression.
#'                 It should accept an expression as input and return a numeric cost.
#' @param termination_cost Numeric. The stopping criteria for evolution.
#'                         If the best cost reaches this value, evolution stops. Default is `NA` (no early stopping).
#' @param optimizer Character. The optimization method to use. Options are:
#'                  - `"es"`: Evolution Strategy (default).
#'                  - `"ga"`: Genetic Algorithm.
#'                  - `"random"`: Randomly selects between `"es"` and `"ga"`, with special handling for iterations.
#' @param iterations Integer. Number of iterations to run the evolutionary process (default: `200`).
#'                   If `optimizer = "random"` and `"es"` is chosen, iterations increase by 10x.
#' @param suggestions A matrix of suggested genomes to initialize the evolution process.
#'                    If `optimizer = "es"`, only the first row is used.
#' @param mutationChance Numeric. Mutation probability in the evolutionary algorithm.
#'                       If `"random"`, a random value in `[0,1]` is chosen.
#'                       Default is `-0.8`, which might need clarification or adjustment.
#' @param verbose Logical. If `TRUE`, prints real-time updates of the evolutionary process (default: `FALSE`).
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{`output_table`}{A `data.frame` containing the best discovered expression and associated metadata:
#'       \itemize{
#'         \item `distribution`: The name of the dataset distribution.
#'         \item `optimizer`: Optimization method used.
#'         \item `best_expression`: The best expression found.
#'         \item `best_cost`: The best cost achieved.
#'         \item `iterations`: The number of iterations performed.
#'         \item `seed`: The random seed used.
#'         \item `mutationChance`: The mutation probability used.
#'         \item `runtime`: The total execution time.
#'       }}
#'     \item{`results_table`}{A `data.frame` containing detailed iteration-wise results:
#'       \itemize{
#'         \item `distribution`: The dataset distribution.
#'         \item `optimizer`: The optimization method used.
#'         \item `expression`: The best expression at a given iteration.
#'         \item `cost`: The corresponding cost.
#'         \item `iterations`: The current iteration number.
#'         \item `mutationChance`: The mutation probability used.
#'         \item `seed`: The random seed used.
#'       }}
#'     \item{`results`}{A list containing detailed evolution process snapshots at each step.}
#'     \item{`best_genome`}{A numeric matrix containing the best genome representation.}
#'   }
#'
#' @examples
#' \dontrun{
#'
#'  data <- generate_distribution_metrics(100, 1000000, "normal")
#'  grammarDef <- create_grammar_wrapper(data, variables = c('sigma'))
#'
#'  SymRegFitFunc <- function(expr) {
#'    result <-unlist(eval(parse(text = expr)))
#'    if (any(is.nan(result))) {
#'      return(Inf)
#'    }
#'    mean(log(1 + abs(data$entropy - result)))
#'  }
#'
#' result <-
#'   symbolic_regression(
#'     grammarDef = grammarDef,
#'     evalFunc = SymRegFitFunc,
#'     iterations = 50,
#'     optimizer = 'ga',
#'     verbose = TRUE,
#'     termination_cost = .1
#'   )
#' }
#' @export
#####
symbolic_regression <- function(grammarDef,
                                evalFunc,
                                termination_cost = NA,
                                optimizer = "es",
                                iterations = 200,
                                suggestions = NULL,
                                mutationChance = NA,
                                verbose = FALSE) {

  if (optimizer == "random") {
    optimizer <- sample(c("es","ga"),size = 1)
    if (optimizer == "es") {
      iterations <- iterations * 10
      if (!is.null(suggestions)) {
        suggestions <- suggestions[1, ]

      }
    }
  }
  if (!is.na(mutationChance) & mutationChance == "random") {
    mutationChance <- runif(1,0,1)
  }


  variablesName <- grammarDef$name
  variables <- grammarDef$variables
  grammarDef <- grammarDef$grammarDef
  known_entropy <- data$entropy


  variables <- c(variables)

  seed <- round(runif(1,1,2147483647))
  set.seed(seed)

  results <- list()

  results_table <- data.frame(
    distribution = character(),
    optimizer =  character(),
    expression = character(),
    cost = numeric(),
    iterations = numeric(),
    mutationChance = numeric(),
    seed = numeric()
  )

  # Start runtime tracking
  start_time <- Sys.time()
  currentBest <- Inf
  # Run Grammatical Evolution
  gramEvolution <- GrammaticalEvolution(
    grammarDef,
    suggestions = suggestions,
    evalFunc = evalFunc,
    terminationCost = termination_cost,
    optimizer = optimizer,
    iterations = iterations,
    mutationChance = mutationChance,
    monitorFunc = function(x) {
      if (x$best$cost < currentBest) {
        x$best$expressions <- gsub("data\\$", "", x$best$expressions)
        if (verbose){
          print(x, show.genome = TRUE)
        }
        best_output <<- x
        currentBest <<- x$best$cost
        currentBestExpressions <<- x$best$expressions
        results[[length(results)+1]] <<- list(x)

        results_table <<-
          rbind(
            results_table,
            data.frame(
              distribution = data$distType[1],
              optimizer = optimizer,
              expression = x$best$expressions,
              cost = x$best$cost,
              iterations = x$population$currentIteration,
              mutationChance = mutationChance,
              seed = seed
            )
          )
        print(results_table)

      }
    }
  )

  # End runtime tracking
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "secs")

  # Extract best expression
  best_expression <- gramEvolution$best$expression

  output_table <- data.frame(
    distribution = data$distType[1],
    optimizer = optimizer,
    best_expression = currentBestExpressions,
    best_cost = currentBest,
    iterations = best_output$population$currentIteration,
    seed = seed,
    mutationChance = mutationChance,
    runtime = runtime
  )

  # Return results
  return(list(
    output_table = output_table,
    results_table = results_table,
    results = results,
    best_genome = matrix(c(currentBest, best_output[["best"]][["genome"]]),nrow = 1)
  ))

}


##### --- ROXYGEN ---- #####
#' Create a Grammar Definition for Symbolic Regression
#'
#' This function constructs a grammar definition for symbolic regression using Grammatical Evolution.
#' It dynamically includes variables from the provided dataset and allows customization of operators and functions.
#'
#' @param data A `data.frame` containing the dataset used for symbolic regression.
#' @param operators Character vector. Mathematical operators to include in the grammar (default: `c("+", "-", "*", "/", "^")`).
#' @param functions List of functions to include in the grammar (default: `c(log, sqrt, exp)`).
#' @param variables Character vector. Names of the columns in `data` to be included as variables in the grammar.
#'
#' @return A list containing:
#'   \describe{
#'     \item{`name`}{The symbolic representation of the selected variables.}
#'     \item{`known_entropy`}{A vector containing the entropy values from `data$entropy`.}
#'     \item{`grammarDef`}{The generated grammar definition for symbolic regression.}
#'     \item{Dynamic Variables}{Additional named elements in the list, where each variable in `variables` is stored as a separate list element with corresponding values from `data`.}
#'   }
#'
#' @details
#' The function constructs a grammar by defining symbolic rules for mathematical expressions.
#' It dynamically integrates variables from `data` while ensuring that constants (1, 2, π) are included in the grammar.
#'
#' @examples
#' \dontrun{
#'   grammar <- create_grammar_wrapper(
#'     data = my_data,
#'     variables = c("x", "y")
#'   )
#' }
#'
#' @export
#####
create_grammar_wrapper <- function(data,
                                   operators = c("+", "-", "*", "/", "^"),
                                   functions = c('log', 'sqrt', 'exp'),
                                   variables,
                                   constants = c(1, 2, pi)) {

  ivCol <- variables  # Preserve variable names
  variables <- syms(paste0("data$",ivCol))

  functions <- syms(functions)
  operators <- syms(operators)

  # Define grammar rules
  ruleDef <- list(
    expr = grule(op(expr, expr), func(expr), var),
    func = do.call(grule, functions),  # Ensure function set includes `exp`
    op = do.call(grule, operators),
    var = do.call(grule, c(variables, constants))
  )

  # Create the grammar definition
  grammarDef <- CreateGrammar(ruleDef)

  return(c(
    list(
      name = variables,
      known_entropy = data$entropy,
      grammarDef = grammarDef
    ),
    setNames(lapply(ivCol, function(col) data[[col]]), ivCol)  # Dynamically add multiple elements
  ))
}

# data <- generate_distribution_metrics(100, 100000,"normal")


# grammarDef <- create_grammar_wrapper(data, variables = "sigma",constants = c(.5,17.07947,2),functions = c('log','exp'),operators = c("^","*"))
#
#
# result <- par_sym_reg(
#   symbolic_regression = symbolic_regression,
#   evalFunc = SymRegFitFunc,
#   grammarDef = grammarDef,
#   data = data,
#   n_core = 14,iterations = 100,optimizer = 'ga'
# )
#
# result <- par_sym_reg(
#   symbolic_regression = symbolic_regression,
#   evalFunc = SymRegFitFunc,
#   grammarDef = grammarDef,
#   data = data,
#   n_core = 14,
#   iterations = 100,
#   optimizer = 'ga',
#   suggestions = result$genome_matrix,
#   mutationChance = .9
# )
#







##### --- ROXYGEN ---- #####
#' Perform Parallel Symbolic Regression
#'
#' This function runs symbolic regression in parallel across multiple cores, leveraging
#' `symbolic_regression()` for the evolutionary process on each core.
#'
#' @inheritParams symbolic_regression
#' @param n_core Integer. The number of CPU cores to use for parallel execution (default: `4`).
#'
#' @return A list containing:
#'   \describe{
#'     \item{`output_table`}{A `data.frame` summarizing the best symbolic expressions found across parallel runs.}
#'     \item{`results_table`}{A `data.frame` containing detailed iteration-wise results from all parallel runs.}
#'     \item{`genome_matrix`}{A matrix containing the top 3 best genome representations sorted by fitness.}
#'     \item{`output_raw`}{A list of raw outputs from each parallel execution of `symbolic_regression()`.}
#'   }
#'
#' @details
#' This function creates a parallel cluster and distributes the symbolic regression tasks
#' across multiple cores. Each core runs `symbolic_regression()` independently on a copy
#' of `grammarDef`. The results are then aggregated into structured output tables.
#'
#' @examples
#' \dontrun{
#'   result <- par_sym_reg(
#'     symbolic_regression = symbolic_regression,
#'     evalFunc = SymRegFitFunc,
#'     grammarDef = grammarDef,
#'     data = data,
#'     n_core = 4
#'   )
#' }
#'
#' @export
#####
par_sym_reg <- function(symbolic_regression = symbolic_regression,
                        evalFunc,
                        grammarDef,
                        data,  # Corrected argument order
                        termination_cost = NA,
                        optimizer = "es",
                        iterations = 50,
                        fit_func,
                        verbose = FALSE,
                        suggestions = NULL,
                        mutationChance = NA,
                        n_core = 4) {  # Moved n_core to the correct place
  # Detect the number of available cores and create cluster
  cl <- parallel::makeCluster(min(parallel::detectCores(), n_core))

  grammarDefLong <- replicate(n_core, grammarDef, simplify = FALSE)


  # Export necessary variables to the cluster
  parallel::clusterExport(
    cl,
    c("symbolic_regression",
      "evalFunc",
      "grammarDef",
      "termination_cost",
      "optimizer",
      "iterations",
      "verbose",
      "suggestions",
      "mutationChance",
      "data"),
    envir = environment()
  )

  # Load required packages on worker nodes
  parallel::clusterEvalQ(cl, {
    library(gramEvol)
  })

  # Run parallel symbolic regression
  opFromPar <- parallel::parLapply(
    cl = cl,
    fun = symbolic_regression,
    X = grammarDefLong,
    evalFunc = evalFunc,
    termination_cost = termination_cost,
    optimizer = optimizer,
    iterations = iterations,
    verbose = verbose,
    suggestions = suggestions,
    mutationChance = mutationChance
  )

  # Stop the cluster after execution
  parallel::stopCluster(cl)

  genome_list <- lapply(opFromPar, function(i) i$best_genome)

  # Convert list of vectors into a matrix with 3 rows
  genome_matrix <- as.matrix(
    as.data.frame(do.call(rbind, genome_list))    %>%
      arrange(V1)    %>%   # Sort by the first column
      distinct(.keep_all = TRUE)    %>%  # Remove duplicates based on V1
      slice_head(n = 3)    %>%
      select(-V1)       # Remove the first column if needed
  )

  output_table <- lapply(opFromPar, function(i)
    i$output_table)

  output_table <-
    as.data.frame(do.call(rbind, output_table) %>%
                    arrange(best_cost))

  results_table <- lapply(opFromPar, function(i)
    i$results_table)

  results_table <-
    as.data.frame(do.call(rbind, results_table) %>%
                    arrange(cost))


  return(list(
    output_table = output_table,
    results_table = results_table,
    genome_matrix = genome_matrix,
    output_raw = opFromPar
  ))
}

# data <- generate_distribution_metrics(100, 1000000, "normal")
#
# grammarDef <- create_grammar_wrapper(data, variables = c('sigma'))
#
# output <-
#   par_sym_reg(
#     symbolic_regression = symbolic_regression,
#     evalFunc = SymRegFitFunc,
#     optimizer = 'random',
#     grammarDef = grammarDef,
#     data = data,
#     iterations = 50,
#     termination_cost = .1,
#     mutationChance = "random",
#     n_core = 15
#   )
#
# output2 <-
#   par_sym_reg(
#     symbolic_regression = symbolic_regression,
#     evalFunc = SymRegFitFunc,
#     optimizer = 'random',
#     grammarDef = grammarDef,
#     data = data,
#     iterations = 50,
#     # termination_cost = .01,
#     suggestions = output$genome_matrix,
#     mutationChance = "random",
#     n_core = 15
#   )
#
#
#
#
#
# test <-
#   symbolic_regression(
#     grammarDef = grammarDef,
#     evalFunc = SymRegFitFunc,
#     iterations = 50,
#     optimizer = 'ga',
#     verbose = TRUE,
#     termination_cost = .2
#   )
#
# result <-
#   symbolic_regression(
#     grammarDef = grammarDef,
#     evalFunc = SymRegFitFunc,
#     iterations = 50,
#     optimizer = 'ga',
#     verbose = TRUE,
#     termination_cost = .1,
#     suggestions = matrix(test[["best_genome"]],nrow = 1)
#   )


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
