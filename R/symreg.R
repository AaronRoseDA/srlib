#these are super changes
rm(list = ls())
library(e1071)
library(rlang)
library(gramEvol)
library(dplyr)


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


SymRegFitFuncExp <- function(expr) {
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

  with(grammar, {
    evalFunc_name <- sample(names(evalFunc), 1)
    evalFunc <- evalFunc[[evalFunc_name]]
    environment(evalFunc) <- environment()

    if (optimizer == "random") {
      optimizer <- sample(c("es", "ga"), 1)
      iterations <- if (optimizer == "es") iterations * 10 else iterations
      if (!is.null(suggestions) && optimizer == "es") suggestions <- suggestions[1, ]
    }

    if (!is.na(mutationChance) && mutationChance == "random") {
      mutationChance <- runif(1, 1 / GrammarMaxSequenceLen(grammarDef, GrammarGetDepth(grammarDef), GrammarStartSymbol(grammarDef)), 1)
    }

    seed <- sample(1:2147483647, 1)
    set.seed(seed)

    results <- list()
    results_table <- data.frame(
      distribution = character(), optimizer = character(), expression = character(),
      cost = numeric(), iterations = numeric(), mutationChance = numeric(), seed = numeric()
    )

    start_time <- Sys.time()
    currentBest <- Inf

    gramEvolution <- GrammaticalEvolution(
      grammarDef,
      suggestions = suggestions,
      evalFunc = evalFunc,
      terminationCost = termination_cost,
      optimizer = optimizer,
      iterations = iterations,
      mutationChance = mutationChance,
      max.depth = maxDepth,
      popSize = popSize,
      monitorFunc = function(x) {
        if (x$best$cost < currentBest) {
          x$best$expressions <- gsub("data\\$", "", x$best$expressions)
          if (verbose) print(x, show.genome = TRUE)
          best_output <<- x
          currentBest <<- x$best$cost
          results[[length(results) + 1]] <<- list(x)

          results_table <<- rbind(
            results_table,
            data.frame(
              distribution = name,
              optimizer = optimizer,
              target = target_name,
              expression = x$best$expressions,
              cost = x$best$cost,
              cost_funciton = evalFunc_name,
              iterations = x$population$currentIteration,
              trial = trial,
              mutationChance = mutationChance,
              target_bias = target_bias,
              possible_funcs = functions,
              possible_opers = operators,
              possible_vars = paste(gsub("data\\$", "", variables), collapse = " "),
              possible_cons = paste(constants),
              maxDepth = maxDepth,
              genome = paste(x$best$genome, collapse = ", "),
              seed = seed,
              date = Sys.time()
            )
          )
          print(tail(results_table[2:6], 1))
        }
      }
    )

    runtime <- difftime(Sys.time(), start_time, units = "secs")

    output_table <- data.frame(
      distribution = name,
      optimizer = optimizer,
      target = target_name,
      best_expression = best_output$best$expressions,
      best_cost = currentBest,
      cost_funciton = evalFunc_name,
      iterations = best_output$population$currentIteration,
      trial = trial,
      mutationChance = mutationChance,
      target_bias = target_bias,
      possible_funcs = functions,
      possible_opers = operators,
      possible_vars = paste(gsub("data\\$", "", variables), collapse = " "),
      possible_cons = paste(constants),
      maxDepth = maxDepth,
      genome = paste(best_output$best$genome, collapse = ", "),
      seed = seed,
      date = Sys.time(),
      runtime = runtime
    )

    suggestions <- matrix(c(currentBest, best_output$best$genome), nrow = 1)

    list(output_table = output_table, results_table = results_table, results = results, best_genome = suggestions)
  })
}


par_sym_reg <- function(symbolic_regression, grammar, n_core = 4) {

  pad_suggestions <- function(grammar) {
    reqSeqLen <- GrammarMaxSequenceLen(grammar$grammarDef, grammar$maxDepth, GrammarStartSymbol(grammar$grammarDef))
    if (!is.null(grammar$suggestions) && ncol(grammar$suggestions) < reqSeqLen) {
      padding <- matrix(0, nrow = nrow(grammar$suggestions), ncol = reqSeqLen - ncol(grammar$suggestions))
      grammar$suggestions <- cbind(grammar$suggestions, padding)
    }
    grammar
  }

  start_trial <- grammar$trial

  for (trial in start_trial:(start_trial + grammar$trials - 1)) {

    grammar <- pad_suggestions(grammar)
    grammar$trial <- trial

    cl <- parallel::makeCluster(min(parallel::detectCores(), n_core))
    grammar_list <- replicate(n_core, grammar, simplify = FALSE)

    parallel::clusterExport(cl, varlist = "symbolic_regression", envir = environment())
    parallel::clusterEvalQ(cl, { library(gramEvol); library(rlang) })

    opFromPar <- parallel::parLapply(cl, grammar_list, symbolic_regression)
    parallel::stopCluster(cl)

    genome_matrix <- do.call(rbind, lapply(opFromPar, function(i) i$best_genome)) %>%
      as.data.frame() %>%
      arrange(V1) %>%
      distinct(.keep_all = TRUE) %>%
      slice_head(n = 5) %>%
      select(-V1) %>%
      as.matrix()

    grammar$output <- rbind(grammar$output, do.call(rbind, lapply(opFromPar, `[[`, "output_table")) %>% arrange(best_cost))
    grammar$results <- rbind(grammar$results, do.call(rbind, lapply(opFromPar, `[[`, "results_table")) %>% arrange(cost))
    grammar$suggestions <- genome_matrix
  }

  grammar$trial <- start_trial + grammar$trials
  grammar$best_output <- grammar$output %>%
    group_by(trial) %>%
    slice_min(best_cost, with_ties = FALSE) %>%
    ungroup()

  last_10 <- grammar$best_output %>%
    arrange(desc(trial)) %>%
    slice_head(n = 10)

  if (nrow(last_10) == 10 && length(unique(last_10$best_expression)) == 1) {
    if (grammar$maxDepth < 7) {
      grammar$maxDepth <- grammar$maxDepth + 1
      message("maxDepth increased to ", grammar$maxDepth)
    } else {
      message("maxDepth already at maximum (6). No change.")
    }
  }


  return(grammar)
}


grammar <- create_grammar_wrapper(
  data = tri,
  target = "entropy_kde",
  target_bias = 0,
  variables = c('kde_var','data_var', 'data_skew', 'data_min', 'data_max'),
  evalFunc <- list(
    SymRegFitFunc = SymRegFitFunc
  ),
  name = 'triangular',
  iterations = 250,
  trials = 20,
  operators = c("-", "+", "/", "*", "^"),
  functions = c("log", "exp", "sqrt"),
  constants = c(.5, 1, 2, pi),
  optimizer = "ga",
  mutationChance = "random",
  popSize = 250
)

# output <- symbolic_regression(grammar = grammar)
# View(output$output_table)
grammar <- par_sym_reg(
  symbolic_regression = symbolic_regression,
  grammar = grammar,
  n_core = 12
)

test <- grammar$output %>% arrange(best_cost, )

test <- grammar$output %>%
  group_by(across(-c(date, runtime,iterations,seed))) %>%
  summarize(
    date = min(date),
    runtime = min(runtime),
    iterations = min(iterations),
    .groups = "drop"
  )
library(dplyr)

test <- grammar$output %>%
  group_by(trial) %>%
  slice_min(best_cost, with_ties = FALSE) %>%
  ungroup()



with(grammar, GrammarGetDepth(grammarDef))

grammar$maxDepth <- 5
parallel_results <- par_sym_reg(
  symbolic_regression = symbolic_regression,
  grammar = grammar,
  n_core = 12
)

export <- grammar$output %>% distinct(best_expression, .keep_all = TRUE)


View(grammar[["output"]])

