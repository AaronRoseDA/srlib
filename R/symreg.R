#these are super changes
rm(list = ls())
library(e1071)

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


generate_distribution_metrics <- function(num_distributions = 1000,
                                          samples_per_dist = 500,
                                          dist_type, a = F, b = F) {

  generateNormalDistributions <- function(num_distributions = 1000,
                                          samples_per_dist = 500,
                                          mu_max = 0, sigma_max = 10000) {
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


SymRegFitFunc <- function(expr) {
  result <-unlist(eval(parse(text = expr)))

  if (any(is.nan(result))) {
    return(Inf)
  }
  mean(log(1 + abs(grammar$target - result + runif(length(grammar$target),min = -target_bias,max = target_bias))))

}

SymRegFitFunc4th <- function(expr) {
  result <- unlist(eval(parse(text = expr)))

  if (any(is.nan(result))) {
    return(Inf)
  }
  mean((grammar$target - result + runif(length(grammar$target), min = -target_bias, max = target_bias))^4)
}

SymRegFitFuncTitr <- function(expr) {
  result <- unlist(eval(parse(text = expr)))

  if (any(is.nan(result))) {
    return(Inf)
  }
  mean(exp(abs(grammar$target - result + runif(length(grammar$target), min = -target_bias, max = target_bias))))
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


library(rlang)
library(rlang)
library(gramEvol)

create_grammar_wrapper <- function(data,
                                   target,
                                   target_bias,
                                   operators = c("+", "-", "*", "/", "^"),
                                   functions = c('log', 'sqrt', 'exp'),
                                   variables,
                                   constants = c(1, 2, pi),
                                   evalFunc,
                                   name = "distribution",
                                   termination_cost = NA,
                                   optimizer = "ga",
                                   iterations = 200,
                                   trials = 1,
                                   trial = 1,
                                   suggestions = NULL,
                                   maxDepth = NA,
                                   popSize = "auto" ,
                                   mutationChance = NA,
                                   verbose = FALSE) {

  # Capture the input data name
  data_name <- as_string(substitute(data))

  # Create expressions like data$var1 dynamically
  variables_expr <- syms(paste0("data$", variables))

  # Prepare functions, operators for grammar
  functions <- syms(functions)
  operators <- syms(operators)

  # Define grammar rules
  ruleDef <- list(
    expr = grule(op(expr, expr), func(expr), var),
    func = do.call(grule, functions),
    op = do.call(grule, operators),
    var = do.call(grule, c(variables_expr, constants))
  )

  # Create grammar
  grammarDef <- CreateGrammar(ruleDef)

  # Package everything into a list
  return(list(
    grammarDef = grammarDef,
    data = data,
    target_name = target,
    variables = variables_expr,
    operators = paste(operators,collapse=" "),
    functions = paste(functions,collapse=" "),
    constants = paste(round(constants,3),collapse=", "),
    target = data[[target]],
    target_bias = target_bias,
    evalFunc = evalFunc,
    name = toupper(name),
    termination_cost = termination_cost,
    optimizer = optimizer,
    iterations = iterations,
    trials = trials,
    trial = 1,
    suggestions = suggestions,
    maxDepth = ifelse(is.na(maxDepth),GrammarGetDepth(grammarDef),maxDepth),
    mutationChance = mutationChance,
    popSize = popSize,
    verbose = verbose,
    output = data.frame(),
    results = data.frame()
  ))
}




symbolic_regression <- function(grammar) {

  grammarDef <- grammar$grammarDef
  data <- grammar$data
  target_name <- grammar$target_name
  target <- grammar$target
  target_bias <- grammar$target_bias
  variables <- grammar$variables
  evalFunc <- grammar$evalFunc
  environment(evalFunc) <- environment()
  name <- grammar$name
  termination_cost <- grammar$termination_cost
  optimizer <- grammar$optimizer
  iterations <- grammar$iterations
  suggestions <- grammar$suggestions
  mutationChance <- grammar$mutationChance
  verbose <- grammar$verbose
  target <- grammar$target
  maxDepth <- grammar$maxDepth


  # Handle random optimizer
  if (optimizer == "random") {
    optimizer <- sample(c("es", "ga"), size = 1)
    if (optimizer == "es") {
      iterations <- iterations * 10
      if (!is.null(suggestions)) {
        suggestions <- suggestions[1, ]
      }
    }
  }

  # Handle random mutation chance
  if (!is.na(mutationChance) && mutationChance == "random") {
    mutationChance <- runif(1, (
      1 /
        GrammarMaxSequenceLen(
          grammar$grammarDef,
          GrammarGetDepth(grammar$grammarDef),
          GrammarStartSymbol(grammar$grammarDef)
        )
    ), 1)
  }

  variables <- c(variables)

  seed <- round(runif(1, 1, 2147483647))
  set.seed(seed)

  results <- list()
  results_table <- data.frame(
    distribution = character(),
    optimizer = character(),
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
      max.depth = maxDepth,
      popSize = grammar$popSize,
      monitorFunc = function(x) {
        if (x$best$cost < currentBest) {
          x$best$expressions <- gsub("data\\$", "", x$best$expressions)
          if (verbose) {
            print(x, show.genome = TRUE)
          }
          best_output <<- x
          currentBest <<- x$best$cost
          currentBestExpressions <<- x$best$expressions
          results[[length(results) + 1]] <<- list(x)

          results_table <<- rbind(
            results_table,
            data.frame(
              distribution = name,
              optimizer = optimizer,
              target = target_name,
              expression = x$best$expressions,
              cost = x$best$cost,
              iterations = x$population$currentIteration,
              trial = grammar$trial,
              mutationChance = mutationChance,
              `target_bias` = grammar$target_bias,
              possible_funcs = grammar$functions,
              possible_opers = grammar$operators,
              possible_vars = paste(gsub("data\\$", "", grammar$variables),collapse=" "),
              possible_cons = paste(grammar$constants),
              maxDepth = maxDepth,
              seed = seed
            )
          )
          print(tail(results_table[2:6],1))
        }
      }
    )

    # End runtime tracking
    end_time <- Sys.time()
    runtime <- difftime(end_time, start_time, units = "secs")

    # Extract best expression
    best_expression <- gramEvolution$best$expression

    output_table <- data.frame(
      distribution = name,
      optimizer = optimizer,
      target = target_name,
      best_expression = best_output$best$expressions,
      best_cost = currentBest,
      iterations = best_output$population$currentIteration,
      trial = grammar$trial,
      mutationChance = mutationChance,
      `target_bias` = grammar$target_bias,
      possible_funcs = grammar$functions,
      possible_opers = grammar$operators,
      possible_vars = paste(gsub("data\\$", "", grammar$variables),collapse=" "),
      possible_cons = paste(grammar$constants),
      maxDepth = maxDepth,
      seed = seed,
      runtime = runtime
    )

    suggestions = matrix(c(currentBest, best_output[["best"]][["genome"]]), nrow = 1)

  # Return results
  return(list(
    output_table = output_table,
    results_table = results_table,
    results = results,
    best_genome = suggestions
  ))
}






par_sym_reg <- function(symbolic_regression = symbolic_regression,
                        grammar,
                        n_core = 4) {

  start_trial <- grammar$trial

  # Loop over trials
  for (trial in start_trial:(start_trial + grammar$trials - 1)) {

    print(trial)

    grammar$trial <- trial  # Update grammar$trial to current trial number

    # Detect number of available cores and create cluster
    cl <- parallel::makeCluster(min(parallel::detectCores(), n_core))

    # Replicate the grammar object across cores
    grammar_list <- replicate(n_core, grammar, simplify = FALSE)

    # Export necessary functions to the cluster
    parallel::clusterExport(
      cl,
      varlist = c("symbolic_regression"),
      envir = environment()
    )

    # Load required packages on worker nodes
    parallel::clusterEvalQ(cl, {
      library(gramEvol)
      library(rlang)
    })

    # Run symbolic_regression in parallel
    opFromPar <- parallel::parLapply(
      cl = cl,
      X = grammar_list,
      fun = function(g) symbolic_regression(g)
    )

    # Stop the cluster after execution
    parallel::stopCluster(cl)

    # Extract best genomes
    genome_list <- lapply(opFromPar, function(i) i$best_genome)

    # Convert list of vectors into a matrix with 3 rows
    genome_matrix <- as.matrix(
      as.data.frame(do.call(rbind, genome_list)) %>%
        arrange(V1) %>%
        distinct(.keep_all = TRUE) %>%
        slice_head(n = 5) %>%
        select(-V1)
    )

    # Combine output tables
    output_table <- do.call(rbind, lapply(opFromPar, function(i) i$output_table)) %>%
      arrange(best_cost)

    grammar$output <<- rbind(grammar$output, output_table)

    # Combine results tables
    results_table <- do.call(rbind, lapply(opFromPar, function(i) i$results_table)) %>%
      arrange(cost)

    grammar$results <<- rbind(grammar$results, results_table)

    grammar$suggestions <<- genome_matrix
  }

  # After the loop, increment grammar$trial so next call knows it ran 4 trials
  grammar$trial <<- start_trial + grammar$trials


  # Return all outputs
  return(list(
    output_table = grammar$output,
    results_table = grammar$result,
    genome_matrix = genome_matrix,
    output_raw = opFromPar
  ))
}



normal <- generate_distribution_metrics(1000,10000,dist_type = "normal")
normal$var <- normal$sigma^2
grammar <- create_grammar_wrapper(
  data = normal,
  target = "entropy",
  target_bias = 0,
  variables = c('var'),
  evalFunc = SymRegFitFunc4th,
  name = 'normal',
  iterations = 100,
  trials = 1,
  operators = c("-","+","/","*","^"),
  functions = c("log","exp","sqrt"),
  constants = c(1,2,pi),
  optimizer = "ga",
  mutationChance = "random",
  maxDepth = 5,popSize = 500
)

# output <- symbolic_regression(grammar = grammar)
# View(output$output_table)

grammar$popSize <- 500
grammar$target_bias <- 0
parallel_results <- par_sym_reg(
  symbolic_regression = symbolic_regression,
  grammar = grammar,
  n_core = 12
)
View(grammar[["output"]])

SymRegFitFuncTitr
# grammar <- create_grammar_wrapper(
#   data = temp,
#   target = "entropy_kde",
#   target_bias = .01,
#   variables = c('data_mean', 'data_var', 'data_min', 'data_max'),
#   evalFunc = SymRegFitFunc,
#   name = 'uniform',
#   iterations = 100,
#   trials = 4,
#   operators = c("-","+","/","*","^"),
#   functions = c("log","exp","sqrt"),
#   constants = c(1,2,pi),
#   optimizer = "ga",
#   mutationChance = "random"
# )
#
# output <- symbolic_regression(grammar = grammar)
# View(output$output_table)
#
#
#
# parallel_results <- par_sym_reg(
#   symbolic_regression = symbolic_regression,
#   grammar = grammar,
#   n_core = 12
# )


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
