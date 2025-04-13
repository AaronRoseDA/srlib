rm(list=ls())
# Install and load the necessary package
library(ks)
library(R.utils)
library(MASS)
library(sn)
library(dplyr)
library(e1071)

myFolder<-"C:\\Users\\kpflugho\\Downloads"
setwd(myFolder)

#dist<-"UNIFORM"
#dist<-"TRIANGULAR"
#dist<-"EXPONENTIAL" #log transformation and reflection
#dist<-"LAPLACE"
#dist<-"LOGISTIC" #seems to be some problems
#dist<-"HALFNORMAL" #seems to be some problems
#dist<-"SKEWNORMAL"
#dist<-"NORMAL"
#dist<-"NEWNORMAL"
#dist<-"NEWUNIFORM"
#dist<-"NEWTRIANGULAR"
#dist<-"NEWSKEWNORMAL"
dist<-"NEWNORMAL2D"
log_transform<-F
reflection<-F

safe_row <- function(template_names, actual_values, mismatch_flag = NA, context = NULL) {
  # Coerce list to named vector if needed
  if (is.list(actual_values) && !is.null(names(actual_values))) {
    actual_values <- unlist(actual_values)
  }

  # Initialize output vector
  out <- setNames(vector("list", length(template_names)), template_names)

  # Track missing or extra fields
  missing_fields <- setdiff(template_names, names(actual_values))
  extra_fields <- setdiff(names(actual_values), template_names)

  # Fill in values where names match
  for (name in template_names) {
    if (name %in% names(actual_values)) {
      out[[name]] <- actual_values[[name]]
    } else {
      out[[name]] <- mismatch_flag
    }
  }

  if (length(missing_fields) > 0 || length(extra_fields) > 0) {
    message("[safe_row] Warning:")
    if (!is.null(context)) {
      message("  Context: ", context)
    }
    if (length(missing_fields) > 0) {
      message("  Missing fields: ", paste(missing_fields, collapse = ", "))
    }
    if (length(extra_fields) > 0) {
      message("  Extra fields: ", paste(extra_fields, collapse = ", "))
    }
  }

  df_out <- as.data.frame(out, stringsAsFactors = FALSE)

  # Attempt robust numeric conversion
  for (col in names(df_out)) {
    if (col %in% c("axis", "dist")) next  # Keep known character columns untouched

    # Convert explicit "NaN", "Inf", and "NA" strings properly
    df_out[[col]][df_out[[col]] %in% c("NaN", "nan")] <- NaN
    df_out[[col]][df_out[[col]] %in% c("Inf", "inf")] <- Inf
    df_out[[col]][df_out[[col]] %in% c("-Inf", "-inf")] <- -Inf
    df_out[[col]][df_out[[col]] %in% c("NA", "na")] <- NA

    suppressWarnings({
      maybe_numeric <- as.numeric(df_out[[col]])
      # Convert if any non-NA value is numeric
      if (any(!is.na(maybe_numeric))) {
        df_out[[col]] <- maybe_numeric
      } else {
        # Ensure fully NA numeric columns are numeric type
        if (all(is.na(maybe_numeric))) {
          df_out[[col]] <- as.numeric(df_out[[col]])
        }
      }
    })
  }

  return(df_out)
}


# Template of all possible columns expected for a full row
template_row_fields <- c(
  "seed", "dim", "count", "start_time", "end_time", "duration_sec",
  "axis", "mu", "sigma", "gamma", "min", "max", "loc", "scale", "shape",
  "data_mean", "data_var", "data_std", "data_skew", "data_min", "data_max",
  "kde_mean", "kde_var", "kde_var_correct", "kde_skew", "kde_skew_correct",
  "entropy_kde", "kde_bandW", "entropy_plugin_data", "entropy_plugin_kde",
  "entropy_theoretical", "corr_theor", "cov_theor", "corr_sample", "cov_sample",
  "corr_kde", "cov_kde", "corr_kde_correct", "cov_kde_correct"
)

# Row construction using safe_row
make_row <- function(row_values, context = NULL) {
  safe_row(template_row_fields, row_values, mismatch_flag = NaN, context = context)
}

# Example usage for univariate case:
make_univariate_row <- function(seed_used, aN, time_start, time_end, duration_sec,
                                axis, theor_moments, data_moments, kde_moments,
                                entropy_plugin_data, entropy_plugin_kde, entropy_theor) {
  make_row(list(
    seed = seed_used,
    dim = 1,
    count = aN,
    start_time = time_start,
    end_time = time_end,
    duration_sec = duration_sec,
    axis = axis,
    mu = theor_moments[["mean"]],
    sigma = theor_moments[["stdev"]],
    gamma = theor_moments[["skew"]],
    min = theor_moments[["min"]],
    max = theor_moments[["max"]],
    loc = theor_moments[["loc"]],
    scale = theor_moments[["scale"]],
    shape = theor_moments[["shape"]],
    data_mean = data_moments[["mean"]],
    data_var = data_moments[["var"]],
    data_std = data_moments[["stdev"]],
    data_skew = data_moments[["skew"]],
    data_min = data_moments[["min"]],
    data_max = data_moments[["max"]],
    kde_mean = kde_moments[["mean"]],
    kde_var = kde_moments[["var"]],
    kde_var_correct = kde_moments[["var_correct"]],
    kde_skew = kde_moments[["skew"]],
    kde_skew_correct = kde_moments[["skew_correct"]],
    entropy_kde = kde_moments[["entropy"]],
    kde_bandW = kde_moments[["bandW"]],
    entropy_plugin_data = entropy_plugin_data,
    entropy_plugin_kde = entropy_plugin_kde,
    entropy_theoretical = entropy_theor,
    corr_theor = NA,
    cov_theor = NA,
    corr_sample = NA,
    cov_sample = NA,
    corr_kde = NA,
    cov_kde = NA,
    corr_kde_correct = NA,
    cov_kde_correct = NA
  ), context = "make_univariate_row")
}

# Generalized wrapper for all univariate row makers
make_skewnormal_row <- make_univariate_row
make_normal_row <- make_univariate_row
make_uniform_row <- make_univariate_row
make_triangular_row <- make_univariate_row

# Row construction for bivariate case
make_bivariate_row <- function(seed_used, aN, time_start, time_end, duration_sec, axis,
                               theor_moments, data_moments, kde_moments,
                               entropy_plugin_data, entropy_theor,
                               corr_theor = NA, cov_theor = NA,
                               corr_sample = NA, cov_sample = NA,
                               corr_kde = NA, cov_kde = NA,
                               corr_kde_correct = NA, cov_kde_correct = NA,
                               context = NULL) {
  print("hello")
  str(kde_moments)
  make_row(list(
    seed = seed_used,
    dim = 2,
    count = aN,
    start_time = time_start,
    end_time = time_end,
    duration_sec = duration_sec,
    axis = axis,
    mu = theor_moments[["mean"]],
    sigma = theor_moments[["stdev"]],
    gamma = theor_moments[["skew"]],
    min = theor_moments[["min"]],
    max = theor_moments[["max"]],
    loc = theor_moments[["loc"]],
    scale = theor_moments[["scale"]],
    shape = theor_moments[["shape"]],
    data_mean = data_moments[["mean"]],
    data_var = data_moments[["var"]],
    data_std = data_moments[["stdev"]],
    data_skew = data_moments[["skew"]],
    data_min = data_moments[["min"]],
    data_max = data_moments[["max"]],
    kde_mean = kde_moments[[paste0("mean_", tolower(axis))]],
    kde_var = kde_moments[[paste0("var_", tolower(axis))]],
    kde_var_correct = kde_moments[[paste0("var_", tolower(axis), "_correct")]],
    kde_skew = kde_moments[[paste0("skew_correct_", tolower(axis))]],
    kde_skew_correct = kde_moments[[paste0("skew_correct_", tolower(axis))]],
    entropy_kde = kde_moments[[paste0("entropy_", tolower(axis))]],
    kde_bandW = kde_moments[[paste0("bandwidth_", tolower(axis))]],
    entropy_plugin_data = entropy_plugin_data,
    entropy_theoretical = entropy_theor,
    corr_theor = corr_theor,
    cov_theor = cov_theor,
    corr_sample = corr_sample,
    cov_sample = cov_sample,
    corr_kde = corr_kde,
    cov_kde = cov_kde,
    corr_kde_correct = corr_kde_correct,
    cov_kde_correct = cov_kde_correct
  ), context = paste0("make_bivariate_", axis, "_row"))
}

make_bivariate_x_row <- function(...) make_bivariate_row(..., axis = "X")
make_bivariate_y_row <- function(...) make_bivariate_row(..., axis = "Y")
make_bivariate_xy_row <- function(seed_used, aN, time_start, time_end, duration_sec, kde_moments,
                                  entropy_plugin_joint, entropy_theor_joint,
                                  rho, cov_xy, data_corr, data_cov) {
  make_row(list(
    seed = seed_used,
    dim = 2,
    count = aN,
    start_time = time_start,
    end_time = time_end,
    duration_sec = duration_sec,
    axis = "XY",
    # explicitly set non-applicable fields to NA
    mu = NA, sigma = NA, gamma = NA, min = NA, max = NA, loc = NA, scale = NA, shape = NA,
    data_mean = NA, data_var = NA, data_std = NA, data_skew = NA, data_min = NA, data_max = NA,
    kde_mean = NA, kde_var = NA, kde_var_correct = NA, kde_skew = NA, kde_skew_correct = NA,
    kde_bandW = NA, entropy_plugin_kde = NA,
    entropy_kde = kde_moments[["entropy_joint"]],
    entropy_plugin_data = entropy_plugin_joint,
    entropy_theoretical = entropy_theor_joint,
    corr_theor = rho,
    cov_theor = cov_xy,
    corr_sample = data_corr,
    cov_sample = data_cov,
    corr_kde = kde_moments[["corr_xy"]],
    cov_kde = kde_moments[["cov_xy"]],
    corr_kde_correct = kde_moments[["corr_xy_corrected"]],
    cov_kde_correct = kde_moments[["cov_xy_corrected"]]
  ), context = "make_bivariate_xy_row")
}

# Objective function
# par = (xi, omega, alpha)
obj_fn_skew_normal <- function(par) {
  xi    <- par[1]
  omega <- par[2]
  alpha <- par[3]
  d_mean <- par[4]
  d_variance <- par[5]
  d_skew <- par[6]
  delta <- alpha / sqrt(1 + alpha^2)

  mu     <- xi + omega * delta * sqrt(2 / pi)
  sigma2 <- omega^2 * (1 - (2 * delta^2 / pi))
  skew   <- ((4 - pi) / 2) * (delta * sqrt(2 / pi))^3 / ((1 - 2 * delta^2 / pi)^(3/2))

  # Squared error loss
  #loss <- (mu - d_mean)^2 + (sigma2 - d_variance)^2 + (skew - d_skew)^2
  loss <- (mu - d_mean)^2 + (sigma2 - d_variance)^2 + 0.1 * (skew - d_skew)^2
  return(loss)
}

# Function to generate random numbers from a Triangular Distribution
rtriangle <- function(n, a, b, c) {
  u <- runif(n)  # Generate n uniform random numbers

  # Apply inverse transform sampling
  x <- ifelse(u < (c - a) / (b - a),
              a + sqrt(u * (b - a) * (c - a)),  # Left of the peak
              b - sqrt((1 - u) * (b - a) * (b - c)))  # Right of the peak
  return(x)
}

#generate random univariate normal data

gen_uni_uniform<-function(aMinimum,aMaximum,aN,aSeed) {
  set.seed(aSeed)
  return(runif(aN, min=aMinimum, max=aMaximum))
}

gen_uni_triangle<-function(aMinimum,aMaximum,aPeak,aN,aSeed) {
  set.seed(aSeed)
  return(rtriangle(aN, a = aMinimum, b = aMaximum, c = aPeak))
}

gen_uni_normal<-function(aMean,aSD,aN,aSeed) {
  set.seed(aSeed)
  return(rnorm(aN, mean=aMean, sd=aSD))
}

gen_uni_half_normal<-function(aMean,aSD,aN,aSeed) {
  set.seed(aSeed)
  return(abs(rnorm(aN, mean=aMean, sd=aSD)))
}

gen_uni_exponential<-function(aLambda,aN,aSeed) {
  set.seed(aSeed)
  return(rexp(aN, rate=aLambda))
}

gen_uni_laplace<-function(aLocation,aScale,aN,aSeed) {
  set.seed(aSeed)
  # u <- runif(aN, -0.5, 0.5)  # Generate n uniform random numbers in (-0.5, 0.5)
  # x <- aLocation - aScale * sign(u) * log(1 - 2 * abs(u))  # Apply inverse transform
  x<-rlaplace(aN, location = aLocation, scale = aScale)
  return(x)
}

gen_uni_logistic<-function(aLocation,aScale,aN,aSeed) {
  set.seed(aSeed)
  return(aN, location = aLocation, scale = aScale)
}

gen_uni_skew_normal<-function(aLocation,aScale,aShape,aN,aSeed) {
  set.seed(aSeed)
  #xi Location parameter
  #omega Scale parameter
  #alpha Shape (skewness)
  return(rsn(aN, xi = aLocation, omega = aScale, alpha = aShape))
}

uni_triangle_skew<-function(aMinimum,aMaximum,aPeak) {
  t_skew<-(sqrt(2)*(aMinimum+aMaximum-2*aPeak) *
             (2*aMinimum-aMaximum-aPeak) * (aMinimum-2*aMaximum+aPeak)) /
    (5*(aMinimum^2+aMaximum^2+aPeak^2-aMinimum*aMaximum-aMinimum*aPeak -
          aMaximum*aPeak)^1.5)
  return(t_skew)
}

uni_skew_normal_skew<-function(aShape) {
  delta<-aShape/sqrt(1+aShape^2)
  return ((4-pi)/2) * (delta*sqrt(2/pi))^3/((1-2*delta^2/pi)^1.5)
}

uni_uniform_entropy<-function(aMinimum,aMaximum) {
  return(log(aMaximum - aMinimum))
}

uni_triangle_entropy<-function(aMinimum,aMaximum) {
  return(1/2+log((aMaximum - aMinimum)/2))
}

uni_normal_entropy<-function(aVariance) {
  return(0.5 * log(2 * pi * exp(1) * aVariance))
}

uni_half_normal_entropy<-function(aVariance) {
  return(log(2 * pi * exp(1) * aVariance)/(2*log(2))-1)
}

uni_exponential_entropy<-function(aLambda) {
  return(1 - log(aLambda))
}

uni_laplace_entropy<-function(aScale) {
  return(log(2*aScale*exp(1)))
}

uni_logistic_entropy<-function(aScale) {
  return(log(aScale)+2)
}


uni_skew_normal_entropy<-function(aLocation,aScale,aShape) {
  #this formula is an approximation Azzalini (1985)
  #xi Location parameter epsilon
  #omega Scale parameter
  #alpha Shape (skewness)
  delta<-aShape/sqrt(1+aShape^2)
  #entropy<-log(aScale*sqrt(2*pi*exp(1))) - delta^2
  entropy <- log(aScale * sqrt(2 * pi * exp(1))) -
    0.5 * log(1 - (2 * delta^2) / pi)
  return(entropy)
}

bi_normal_entropy<-function(aVarX,aVarY,aCorr) {

  return(1+log(2*pi)+1/2*log(aVarX*aVarY*(1-aCorr^2)))
}

set_theor_moments<-function(aMean,aVar,aSkew,aMin,aMax,aLoc=NA,aScale=NA,aShape=NA) {
  pMom<-c(mean=aMean,stdev=sqrt(aVar),var=aVar,skew=aSkew,
          min=aMin,max=aMax,loc=aLoc,scale=aScale,shape=aShape)
  return(pMom)
}

calc_sample_moments<-function(aData) {
  sMom<-c(mean=0,stdev=0,var=0,skew=0,min=0,max=0)
  sMom["mean"]<-mean(aData)
  sMom["var"]<-var(aData)
  sMom["stdev"]<-sd(aData)
  sMom["skew"] <- e1071::skewness(aData, type = 2)
  sMom["min"]<-min(aData)
  sMom["max"]<-max(aData)
  return(sMom)
}

calc_kde_moments <- function(aInput, margin = NULL, aGrid = 1024) {
  if (is.vector(aInput)) {
    # Univariate raw data input
    kde_result <- kde(aInput, gridsize = aGrid)
    grid <- kde_result$eval.points
    f_k <- kde_result$estimate
    dx <- diff(grid[1:2])
  } else if (inherits(aInput, "kde")) {
    # Project marginal from multivariate KDE
    if (is.null(margin)) stop("Specify margin index when using a multivariate KDE object.")
    kde_result <- aInput
    grid <- kde_result$eval.points[[margin]]
    dx <- diff(grid[1:2])

    # Integrate out other dimensions to get marginal
    f_k <- apply(kde_result$estimate, margin, sum)
    for (d in seq_along(kde_result$eval.points)) {
      if (d != margin) {
        dx_other <- diff(kde_result$eval.points[[d]][1:2])
        f_k <- f_k * dx_other
      }
    }
  } else {
    stop("Input must be either a univariate vector or a 'kde' object with a specified margin.")
  }

  # Mean
  kde_mean <- sum(grid * f_k) * dx

  # Variance and skewness
  var_raw <- sum((grid - kde_mean)^2 * f_k) * dx
  skew_raw <- sum((grid - kde_mean)^3 * f_k) * dx

  # Bias correction
  #h <- if (is.null(margin)) kde_result$h else kde_result$H[margin, margin]
  if (is.null(margin)) {
    # Univariate: scalar bandwidth
    h <- kde_result$h
  } else {
    # Multivariate: bandwidth matrix
    if (is.matrix(kde_result$H)) {
      h <- kde_result$H[margin, margin]
    } else if (length(kde_result$H) == 1) {
      h <- kde_result$H
    } else {
      stop("Bandwidth H has unexpected structure.")
    }
  }

  var_correct <- var_raw - h^2
  skew_correct <- if (var_correct > 0) skew_raw / (sqrt(var_correct)^3) else NA

  # Entropy
  entropy <- tryCatch({
    if (is.vector(aInput)) {
      density_vals <- predict(kde_result, x = aInput)
      -mean(log(density_vals + 1e-12))
    } else {
      -sum(f_k * log(f_k + 1e-12)) * dx
    }
  }, error = function(e) NA)

  return(c(
    mean = kde_mean,
    var = var_raw,
    var_correct = var_correct,
    skew = skew_raw,
    skew_correct = skew_correct,
    entropy = entropy,
    bandW = h
  ))
}


calc_kde_moments_2d <- function(aData, aGrid = 1024) {
  if (!is.matrix(aData) || ncol(aData) != 2) {
    stop("Input data must be a matrix with two columns (bivariate).")
  }

  # Run 2D KDE
  kde_result <- kde(x = aData, gridsize = rep(aGrid, 2))

  # Marginals via projection
  kde_x <- calc_kde_moments(kde_result, margin = 1)
  kde_y <- calc_kde_moments(kde_result, margin = 2)

  # Extract grid and density
  x_grid <- kde_result$eval.points[[1]]
  y_grid <- kde_result$eval.points[[2]]
  f_xy <- kde_result$estimate
  dx <- diff(x_grid[1:2])
  dy <- diff(y_grid[1:2])
  dA <- dx * dy

  x_mat <- matrix(rep(x_grid, each = aGrid), nrow = aGrid)
  y_mat <- matrix(rep(y_grid, times = aGrid), nrow = aGrid)

  mean_x <- sum(x_mat * f_xy) * dA
  mean_y <- sum(y_mat * f_xy) * dA
  x_centered <- x_mat - mean_x
  y_centered <- y_mat - mean_y

  # Raw covariances
  var_xx <- sum(x_centered^2 * f_xy) * dA
  var_yy <- sum(y_centered^2 * f_xy) * dA
  cov_xy <- sum(x_centered * y_centered * f_xy) * dA
  cov_matrix_raw <- matrix(c(var_xx, cov_xy, cov_xy, var_yy), nrow = 2)
  corr_xy <- cov_xy / sqrt(var_xx * var_yy)

  # Bias correction
  if (!is.null(kde_result$H)) {
    bias_matrix <- kde_result$H %*% t(kde_result$H)
    cov_matrix_corrected <- cov_matrix_raw - bias_matrix
    var_x_corrected <- cov_matrix_corrected[1, 1]
    var_y_corrected <- cov_matrix_corrected[2, 2]
    cov_xy_corrected <- cov_matrix_corrected[1, 2]
    corr_xy_corrected <- if (var_x_corrected > 0 && var_y_corrected > 0) {
      cov_xy_corrected / sqrt(var_x_corrected * var_y_corrected)
    } else { NA }
  } else {
    cov_matrix_corrected <- cov_matrix_raw
    var_x_corrected <- var_y_corrected <- cov_xy_corrected <- corr_xy_corrected <- NA
  }

  # Entropy
  density_vals_joint <- predict(kde_result, x = aData)
  entropy_joint <- -mean(log(density_vals_joint + 1e-12))
  entropy_x <- kde_x[["entropy"]]
  entropy_y <- kde_y[["entropy"]]
  mutual_info <- entropy_x + entropy_y - entropy_joint

  return(list(
    mean_x = kde_x[["mean"]],
    mean_y = kde_y[["mean"]],
    var_x = var_xx,
    var_y = var_yy,
    var_x_correct = var_x_corrected,
    var_y_correct = var_y_corrected,
    skew_correct_x = kde_x[["skew_correct"]],
    skew_correct_y = kde_y[["skew_correct"]],
    entropy_x = entropy_x,
    entropy_y = entropy_y,
    bandwidth_x = kde_x[["bandW"]],
    bandwidth_y = kde_y[["bandW"]],
    cov_xy = cov_xy,
    cov_xy_corrected = cov_xy_corrected,
    corr_xy = corr_xy,
    corr_xy_corrected = corr_xy_corrected,
    mean_vec = c(mean_x, mean_y),
    cov_matrix_raw = cov_matrix_raw,
    cov_matrix_corrected = cov_matrix_corrected,
    entropy_joint = entropy_joint,
    mutual_information = mutual_info
  ))
}


showStats<-function(aN,aTheor_moments,aData_moments,aKde_moments,
                    aData_entropy_plugin,aKde_entropy_plugin,
                    aEntropy_theor) {
  cat("\nNormal Distribution\n")
  cat(sprintf("\nData N: %7d\n",aN))
  cat("Theoretical Moments\n")
  print(stack(theor_moments))
  cat("Data Moments\n")
  print(stack(data_moments))
  cat("KDE Moments\n")
  print(stack(kde_moments))

  cat("Data plug-in Entropy:", data_entropy_plugin, "\n")
  cat("KDE plug-in Entropy:", kde_entropy_plugin, "\n")
  cat("Theoretical Entropy:", entropy_theor, "\n")
}


# Refactored exper_* and run_simulations* functions using make_*_row

exper_normal <- function(aN, aGrid, aM1, aM2, showAll) {
  seed_used <- sample.int(.Machine$integer.max, 1)
  set.seed(seed_used)
  time_start <- Sys.time()

  theor_moments <- set_theor_moments(aM1, aM2, 0, -Inf, Inf)
  data <- gen_uni_normal(theor_moments["mean"], theor_moments["stdev"], aN, seed_used)
  data_moments <- calc_sample_moments(data)
  kde_moments <- calc_kde_moments(data, aGrid)

  entropy_theor <- uni_normal_entropy(theor_moments[["var"]])
  entropy_plugin_data <- uni_normal_entropy(data_moments[["var"]])
  entropy_plugin_kde <- uni_normal_entropy(kde_moments[["var_correct"]])

  if (showAll) showStats(aN, theor_moments, data_moments, kde_moments,
                         entropy_plugin_data, entropy_plugin_kde, entropy_theor)

  time_end <- Sys.time()
  duration_sec <- as.numeric(difftime(time_end, time_start, units = "secs"))
  make_normal_row(seed_used, aN, time_start, time_end, duration_sec, "X",
                  theor_moments, data_moments, kde_moments,
                  entropy_plugin_data, entropy_plugin_kde, entropy_theor)
}

exper_uniform <- function(aN, aGrid, aM4, aM5, showAll) {
  seed_used <- sample.int(.Machine$integer.max, 1)
  set.seed(seed_used)
  time_start <- Sys.time()

  theor_mean <- (aM4 + aM5) / 2
  theor_variance <- (aM5 - aM4)^2 / 12
  theor_moments <- set_theor_moments(theor_mean, theor_variance, 0, aM4, aM5)
  data <- gen_uni_uniform(aM4, aM5, aN, seed_used)
  data_moments <- calc_sample_moments(data)
  kde_moments <- calc_kde_moments(data, aGrid)

  entropy_theor <- uni_uniform_entropy(aM4, aM5)
  entropy_plugin_data <- uni_uniform_entropy(data_moments[["min"]], data_moments[["max"]])
  entropy_plugin_kde <- entropy_plugin_data

  if (showAll) showStats(aN, theor_moments, data_moments, kde_moments,
                         entropy_plugin_data, entropy_plugin_kde, entropy_theor)

  time_end <- Sys.time()
  duration_sec <- as.numeric(difftime(time_end, time_start, units = "secs"))
  make_uniform_row(seed_used, aN, time_start, time_end, duration_sec, "X",
                   theor_moments, data_moments, kde_moments,
                   entropy_plugin_data, entropy_plugin_kde, entropy_theor)
}

exper_triangular <- function(aN, aGrid, aM1, aM4, aM5, showAll) {
  seed_used <- sample.int(.Machine$integer.max, 1)
  set.seed(seed_used)
  time_start <- Sys.time()

  theor_mean <- (aM4 + aM5 + aM1) / 3
  theor_variance <- (aM5^2 + aM4^2 + aM1^2 - aM4 * aM5 - aM4 * aM1 - aM5 * aM1) / 18
  theor_skew <- uni_triangle_skew(aM4, aM5, aM1)
  theor_moments <- set_theor_moments(theor_mean, theor_variance, theor_skew, aM4, aM5, aM1)
  data <- gen_uni_triangle(aM4, aM5, aM1, aN, seed_used)
  data_moments <- calc_sample_moments(data)
  kde_moments <- calc_kde_moments(data, aGrid)

  entropy_theor <- uni_triangle_entropy(aM4, aM5)
  entropy_plugin_data <- uni_triangle_entropy(data_moments[["min"]], data_moments[["max"]])
  entropy_plugin_kde <- entropy_plugin_data

  if (showAll) showStats(aN, theor_moments, data_moments, kde_moments,
                         entropy_plugin_data, entropy_plugin_kde, entropy_theor)

  time_end <- Sys.time()
  duration_sec <- as.numeric(difftime(time_end, time_start, units = "secs"))
  make_triangular_row(seed_used, aN, time_start, time_end, duration_sec, "X",
                      theor_moments, data_moments, kde_moments,
                      entropy_plugin_data, entropy_plugin_kde, entropy_theor)
}

exper_skew_normal <- function(aN, aGrid, aM1, aM2, aM3, showAll) {
  seed_used <- sample.int(.Machine$integer.max, 1)
  set.seed(seed_used)
  time_start <- Sys.time()

  theor_delta <- aM3 / sqrt(1 + aM2^2)
  theor_mean <- aM1 + aM2 * theor_delta * sqrt(2 / pi)
  theor_variance <- aM2^2 * (1 - 2 * theor_delta^2 / pi)
  theor_skew <- uni_skew_normal_skew(aM3)
  theor_moments <- set_theor_moments(theor_mean, theor_variance, theor_skew, -Inf, Inf, aM1, aM2, aM3)
  data <- gen_uni_skew_normal(aM1, aM2, aM3, aN, seed_used)
  data_moments <- calc_sample_moments(data)
  kde_moments <- calc_kde_moments(data, aGrid)

  entropy_theor <- uni_skew_normal_entropy(aM1, aM2, aM3)
  #uni_normal_entropy called as an upper bound
  entropy_plugin_data <- uni_normal_entropy(data_moments[["var"]])
  entropy_plugin_kde <- NA

  if (showAll) showStats(aN, theor_moments, data_moments, kde_moments,
                         entropy_plugin_data, entropy_plugin_kde, entropy_theor)

  time_end <- Sys.time()
  duration_sec <- as.numeric(difftime(time_end, time_start, units = "secs"))
  make_skewnormal_row(seed_used, aN, time_start, time_end, duration_sec, "X",
                      theor_moments, data_moments, kde_moments,
                      entropy_plugin_data, entropy_plugin_kde, entropy_theor)
}


exper_bivariate_normal <- function(aN, aGrid, mu_vec, sigma_vec, rho = 0, showAll = TRUE) {
  seed_used <- sample.int(.Machine$integer.max, 1)
  set.seed(seed_used)
  time_start <- Sys.time()

  # Theoretical values
  theor_moments_x <- set_theor_moments(mu_vec[1], sigma_vec[1]^2, 0, -Inf, Inf)
  theor_moments_y <- set_theor_moments(mu_vec[2], sigma_vec[2]^2, 0, -Inf, Inf)
  cov_xy <- rho * sigma_vec[1] * sigma_vec[2]
  cov_matrix <- matrix(c(sigma_vec[1]^2, cov_xy, cov_xy, sigma_vec[2]^2), nrow = 2)

  # Sample data and empirical stats
  data <- MASS::mvrnorm(n = aN, mu = mu_vec, Sigma = cov_matrix)
  data_corr <- cor(data[, 1], data[, 2])
  data_cov <- cov(data)[1, 2]
  data_moments_x <- calc_sample_moments(data[, 1])
  data_moments_y <- calc_sample_moments(data[, 2])
  kde_moments <- calc_kde_moments_2d(data, aGrid)

  # Entropies
  entropy_theor_x <- uni_normal_entropy(theor_moments_x[["var"]])
  entropy_theor_y <- uni_normal_entropy(theor_moments_y[["var"]])
  entropy_theor_joint <- log(2 * pi * exp(1)) + 0.5 * log(det(cov_matrix))

  entropy_plugin_x <- uni_normal_entropy(data_moments_x[["var"]])
  entropy_plugin_y <- uni_normal_entropy(data_moments_y[["var"]])
  entropy_plugin_joint <- log(2 * pi * exp(1)) + 0.5 * log(det(cov(data)))

  time_end <- Sys.time()
  duration_sec <- as.numeric(difftime(time_end, time_start, units = "secs"))

  # Assemble rows
  #print("x_row")
  #str(kde_moments)
  x_row <- make_bivariate_x_row(
    seed_used, aN, time_start, time_end, duration_sec,
    theor_moments_x, data_moments_x, kde_moments,
    entropy_plugin_x, entropy_theor_x
  )
  #print("y_row")
  #str(kde_moments)
  y_row <- make_bivariate_y_row(
    seed_used, aN, time_start, time_end, duration_sec,
    theor_moments_y, data_moments_y, kde_moments,
    entropy_plugin_y, entropy_theor_y
  )
  #print("xy_row")
  #str(kde_moments)
  xy_row <- make_bivariate_xy_row(
    seed_used, aN, time_start, time_end, duration_sec, kde_moments,
    entropy_plugin_joint, entropy_theor_joint,
    rho, cov_xy, data_corr, data_cov
  )

  return(rbind(x_row, y_row, xy_row))
}


run_simulations <- function(param_grid, n_sim = 10, aDist) {
  do.call(rbind, lapply(1:nrow(param_grid), function(i) {
    m1 <- param_grid$aM1[i]
    m2 <- param_grid$aM2[i]
    m3 <- param_grid$aM3[i]
    m4 <- param_grid$aM4[i]
    m5 <- param_grid$aM5[i]

    sim_df <- do.call(rbind, lapply(1:n_sim, function(j) {
      if (aDist == "NEWNORMAL") {
        result <- exper_normal(100000, 1024, m1, m2, FALSE)
      } else if (aDist == "NEWUNIFORM") {
        result <- exper_uniform(100000, 1024, m4, m5, FALSE)
      } else if (aDist == "NEWTRIANGULAR") {
        result <- exper_triangular(100000, 1024, m1, m4, m5, FALSE)
      } else if (aDist == "NEWSKEWNORMAL") {
        result <- exper_skew_normal(100000, 1024, m1, m2, m3, FALSE)
      }

      result$group_id <- i
      result$sim_id <- j
      result$dist <- aDist
      result
    }))

    sim_df
  }))
}


run_simulations_bivariate <- function(param_grid, n_sim = 10, aDist = "NEWNORMAL") {
  do.call(rbind, lapply(1:nrow(param_grid), function(i) {
    m1 <- param_grid$aM1[i]
    m2 <- param_grid$aM2[i]
    m3 <- param_grid$aM3[i]
    m4 <- param_grid$aM4[i]
    m5 <- param_grid$aM5[i]
    m6 <- param_grid$aM6[i]
    m7 <- param_grid$aM7[i]
    m8 <- param_grid$aM8[i]
    m9 <- param_grid$aM9[i]
    m10 <- param_grid$aM10[i]
    m11 <- param_grid$aM11[i]

    do.call(rbind, lapply(1:n_sim, function(j) {
      sim_result <- exper_bivariate_normal(
        aN = 1e5, aGrid = 1024,
        mu_vec = c(m1, m6),
        sigma_vec = c(sqrt(m2), sqrt(m7)),  # If m2/m7 are variances
        rho = m11, showAll = FALSE
      )

      sim_result$group_id <- i
      sim_result$sim_id <- j
      sim_result$dist <- aDist

      return(sim_result)
    }))
  }))
}


n_sim <- 30

set_normal_param<-function() {
  m1_vals <- c(0, 1, 2, 3, 4, 5) # mean
  m2_vals <- c(1^2, 2^2, 3^2, 4^2, 5^2, 6^2) # standard deviation
  m3_vals <- c(0) # skewness
  m4_vals <- c(-Inf) # minimum
  m5_vals <- c(Inf) # maximum
  return(list(m1_vals,m2_vals,m3_vals,m4_vals,m5_vals))
}

set_skew_normal_param<-function() {
  #solving the mean, variance and skewness from
  #xi, omega, and alpha is straightforward and
  #known as the method of moments; the inverse is not straightforward
  #it requires non-linear optimization and recovery can be iffy
  m1_vals <- c(0, 2, 4, 5, 6, 8) # loc (xi)
  m2_vals <- c(2, 3, 4, 6, 8, 10) # scale (omega)
  m3_vals <- c(-2,-1,0,0.5,1.5,2) # shape (alpha)
  m4_vals <- c(-Inf) # minimum
  m5_vals <- c(Inf) # maximum
  return(list(m1_vals,m2_vals,m3_vals,m4_vals,m5_vals))
}

set_uniform_param<-function() {
  m1_vals <- c(NA) # mean
  m2_vals <- c(NA) # standard deviation
  m3_vals <- c(0) # skewness
  m4_vals <- c(-1,-2,-3,-4,-5,-6) # minimum
  m5_vals <- c(1,2,4,6,8,10) # maximum
  return(list(m1_vals,m2_vals,m3_vals,m4_vals,m5_vals))
}

set_triangular_param<-function() {
  # the peak must be within the min/max values
  m1_vals <- c(1,2,3) # location/peak
  m2_vals <- c(NA) # standard deviation
  m3_vals <- c(NA) # skewness
  m4_vals <- c(-1,-2,-4) # minimum
  m5_vals <- c(5,7,9) # maximum
  return(list(m1_vals,m2_vals,m3_vals,m4_vals,m5_vals))
}

set_normal_param_bivariate <- function() {
  m1_vals <- c(0, 1, 2) # x means
  m2_vals <- c(1, 2) # x standard deviations
  m3_vals <- c(0) # x skewness
  m4_vals <- c(-Inf) # x minimum
  m5_vals <- c(Inf) # x maximum
  m6_vals <- c(0, 1, 2) # y means
  m7_vals <- c(1, 2) # y standard deviations
  m8_vals <- c(0) # y skewness
  m9_vals <- c(-Inf) # y minimum
  m10_vals <- c(Inf) # y maximum
  m11_vals <- c(0, 0.5, 0.9)  # Correlation values

  # m1_vals <- c(0) # x means
  # m2_vals <- c(1) # x standard deviations
  # m3_vals <- c(0) # x skewness
  # m4_vals <- c(-Inf) # x minimum
  # m5_vals <- c(Inf) # x maximum
  # m6_vals <- c(3) # y means
  # m7_vals <- c(4) # y standard deviations
  # m8_vals <- c(0) # y skewness
  # m9_vals <- c(-Inf) # y minimum
  # m10_vals <- c(Inf) # y maximum
  # m11_vals <- c(0.5)  # Correlation values

  return(list(m1_vals,m2_vals,m3_vals,m4_vals,m5_vals,
              m6_vals,m7_vals,m8_vals,m9_vals,m10_vals,m11_vals))

}

#####################################################################
### UNIVARIATE UNIVARIATE UNIVARIATE UNIVARIATE UNIVARIATE UNIVARIATE
#####################################################################

if(dist=="NEWNORMAL") {
  m<-set_normal_param()
  myFileName<-"normal.csv"
} else if(dist=="NEWUNIFORM") {
  m<-set_uniform_param()
  myFileName<-"uniform.csv"
} else if(dist=="NEWTRIANGULAR") {
  m<-set_triangular_param()
  myFileName<-"triangular.csv"
} else if(dist=="NEWSKEWNORMAL") {
  m<-set_skew_normal_param()
  myFileName<-"skewnormal.csv"
}

param_grid <- expand.grid(aM1 = m[[1]], aM2 = m[[2]],
                          aM3 = m[[3]], aM4 = m[[4]],
                          aM5 = m[[5]])

all_results <- run_simulations(param_grid, n_sim, dist)
results_df <- as.data.frame(all_results)
results_df <- results_df %>%
  select(group_id,sim_id, everything())
results_df <- results_df %>%
  mutate(
    dist = dist,
    dim = 1,
    axis = "X"
  ) %>%
  relocate(dist, dim, axis,.after = sim_id)
results_df$duration_sec <- as.numeric(results_df$duration_sec)



#  View(results_df)


cols_to_exclude <- c("group_id", "sim_id", "dist", "dim", "axis", "seed", "start_time", "end_time", "duration_sec")

summary_df <- results_df %>%
  group_by(group_id) %>%
  summarise(
    dist = first(dist),
    dim = first(dim),
    axis = first(axis),
    count = n(),
    start_time = first(start_time),
    end_time = last(end_time),
    duration_sec = sum(duration_sec),
    across(
      .cols = setdiff(names(results_df), cols_to_exclude),
      .fns = ~ mean(., na.rm = TRUE),
      .names = "{.col}"
    ),
    .groups = "drop"
  )
#View(summary_df)
write.csv(summary_df,file=myFileName,row.names=F)

#####################################################################
### BIVARIATE BIVARIATE BIVARIATE BIVARIATE BIVARIATE BIVARIATE
#####################################################################

if (dist == "NEWNORMAL2D") {
  m <- set_normal_param_bivariate()
  myFileName <- "bivariate_normal.csv"
  param_grid <- expand.grid(
    aM1 = m[[1]],  # mu_x
    aM2 = m[[2]],  # sigma_x
    aM3 = m[[3]],  # skew_y
    aM4 = m[[4]],  # min_y
    aM5 = m[[5]],   # max_y
    aM6 = m[[6]],  # mu_y
    aM7 = m[[7]],  # sigma_y
    aM8 = m[[8]],  # skew_y
    aM9 = m[[9]],  # min_y
    aM10 = m[[10]],   # max_y
    aM11 = m[[11]]   # rho
  )

  all_results <- run_simulations_bivariate(param_grid, n_sim = 10, dist)
  results_df <- as.data.frame(all_results)

  results_df <- results_df %>%
    select(group_id, sim_id, everything()) %>%
    mutate(
      dist = dist,
      dim = 2
      #axis = "XY"  # distinguishes joint result from marginal ones
    ) %>%
    relocate(dist, dim, axis, .after = sim_id)
  rownames(results_df) <- NULL
  # View(results_df)
  # Columns to exclude from numeric averaging
  cols_to_exclude <- c(
    "group_id", "sim_id", "dist", "dim", "axis",
    "seed", "start_time", "end_time", "duration_sec"
  )
  #View(results_df)
  summary_df <- results_df %>%
    group_by(group_id, axis) %>%  # now groups by both group_id and axis
    summarise(
      dist = first(dist),
      dim = first(dim),
      count = n(),
      start_time = first(start_time),
      end_time = last(end_time),
      duration_sec = sum(duration_sec),
      across(
        .cols = setdiff(names(results_df), cols_to_exclude),
        .fns =  ~mean(.x, na.rm = TRUE),
        .names = "{.col}"
      ),
      .groups = "drop"
    )
  #View(summary_df)
  write.csv(summary_df, file = myFileName, row.names = FALSE)
}




sample_size<-100000

if(dist == "NORMAL") {
  log_transform<-FALSE
  reflection<-FALSE
  theor_mean<-0
  theor_sd<-1
  theor_variance<-theor_sd^2
  theor_skew<-0
  data<-gen_uni_normal(theor_mean,theor_sd,
                       sample_size,8675309)
} else if(dist == "HALFNORMAL") {
  log_transform<-FALSE
  reflection<-FALSE
  theor_center<-0
  theor_center_sd<-1
  theor_mean<-theor_center_sd * sqrt(2) / sqrt(pi)
  theor_variance<-theor_center_sd^2*(1-2/pi)
  theor_sd <-theor_variance^0.5

  data<-gen_uni_half_normal(theor_center,theor_center_sd,
                            sample_size,8675309)
} else if(dist == "UNIFORM") {
  log_transform<-FALSE
  reflection<-FALSE
  theor_minimum<-0
  theor_maximum<-10
  theor_mean<-(theor_minimum+theor_maximum)/2
  theor_variance<-(theor_maximum - theor_minimum)^2/12
  theor_sd<-sqrt(theor_variance)
  theor_skew<-0
  data<-gen_uni_uniform(theor_minimum,theor_maximum,
                        sample_size,8675309)
} else if(dist == "TRIANGULAR") {
  log_transform<-FALSE
  reflection<-FALSE
  theor_minimum<-0
  theor_maximum<-10
  theor_peak<-3
  theor_mean<-(theor_minimum+theor_maximum+theor_peak)/3
  theor_variance<-(theor_maximum^2 + theor_minimum^2 + theor_peak^2 -
                     theor_minimum*theor_maximum -
                     theor_minimum*theor_peak -
                     theor_maximum*theor_peak)/18
  theor_sd<-sqrt(theor_variance)
  theor_skew<-uni_triangle_skew(theor_minimum,theor_maximum,theor_peak)
  data<-gen_uni_triangle(theor_minimum,theor_maximum,theor_peak,
                         sample_size,8675309)

} else if(dist == "EXPONENTIAL") {
  log_transform<-TRUE
  reflection<-TRUE
  theor_lambda<-2.0
  theor_mean<-1/theor_lambda
  theor_variance<-1/theor_lambda^2
  theor_sd<-sqrt(theor_variance)
  data<-gen_uni_exponential(theor_lambda,
                            sample_size,42)
} else if(dist == "LAPLACE") {
  log_transform<-FALSE
  reflection<-FALSE
  theor_location<-2
  theor_scale<-4
  theor_mean<-theor_location
  theor_variance<-2*theor_scale^2
  theor_sd<-sqrt(theor_variance)
  # data<-gen_uni_laplace(theor_location,theor_scale,
  #                           sample_size,8675309)
  data<-gen_uni_laplace(theor_location,theor_scale,
                        sample_size,8675309)
} else if(dist == "LOGISTIC") {
  log_transform<-FALSE
  reflection<-FALSE
  theor_location<-2
  theor_scale<-4
  theor_mean<-theor_location
  theor_variance<-theor_scale^2*pi^2/3
  theor_sd<-sqrt(theor_variance)
  data<-gen_uni_laplace(theor_location,theor_scale,
                        sample_size,8675309)
} else if(dist == "SKEWNORMAL") {
  log_transform<-FALSE
  reflection<-FALSE
  theor_location<-0
  theor_scale<-1
  theor_shape<-1
  #the theoretical location, scale and shape are not
  #recoverable from the estimated moments.
  #However, the solved theoretical moments should equal data moments.
  theor_delta<-theor_shape/sqrt(1+theor_shape^2)
  theor_mean<-theor_location + theor_scale*theor_delta*sqrt(2/pi)
  theor_variance<-theor_scale^2*(1-2*theor_delta^2/pi)
  theor_sd<-sqrt(theor_variance)
  theor_skew<-uni_skew_normal_skew(theor_shape)
  data<-gen_uni_skew_normal(theor_location,theor_scale,theor_shape,
                            sample_size,8675309)
  #data <- 'attributes<-'(as.vector(data), NULL)
  data <- as.vector(unclass(data))
} else if(dist == "BINORMAL") {
  log_transform<-FALSE
  reflection<-FALSE
  theor_means<- c(50, 80)
  theor_sigmas<-c(3,10)
  theor_correlation<-0.6
  theor_variances<-theor_sigmas^2
  data<-gen_uni_laplace(theor_location,theor_scale,
                        sample_size,8675309)
}

if(dist != "NEWNORMAL") {
  if(log_transform==FALSE) {
    data_mean<-mean(data)
    data_variance<-var(data)
    data_skew <- skewness(data, type = 2)
    data_minimum<-min(data)
    data_maximum<-max(data)

    # Perform KDE with default bandwidth
    # uses a gaussian kernel; no option to change
    kde_result <- kde(data, gridsize = 1024)
    x_vals <- kde_result$eval.points
    density_vals <- kde_result$estimate
  } else {
    log_data<-log(data)
    data_mean<-mean(data)
    data_variance<-var(data)
    data_skew <- skewness(data, type = 2)
    data_minimum<-min(data)
    data_maximum<-max(data)
    # hscv is not always stable for large samples
    h_pilot <- hscv(log_data)
    h_pilot<-0.1013111
    h_opt <- h_pilot * 1.05
    #KDE on Log-Data
    kde_result <- kde(x = log_data, h = h_opt, adj = 1.05, gridsize = 2048)
    # Extract and back-transform KDE results
    x_vals <- kde_result$eval.points
    density_vals<-kde_result$estimate
    exp_x_vals <- exp(x_vals)
    exp_density_vals <- density_vals / exp_x_vals
  }
}
if(dist != "NEWNORMAL") {
  if(reflection==TRUE) {
    # Step 4b: Reflection around 0 to handle boundary for exponential distribution
    reflected_x_vals <- -exp_x_vals[exp_x_vals < 1]
    reflected_density_vals <- exp_density_vals[exp_x_vals < 1]

    # Combine original and reflected values
    combined_x_vals <- c(reflected_x_vals, exp_x_vals)
    combined_density_vals <- c(reflected_density_vals, exp_density_vals)

    # Sort combined values
    sorted_indices <- order(combined_x_vals)
    combined_x_vals <- combined_x_vals[sorted_indices]
    combined_density_vals <- combined_density_vals[sorted_indices]
    #x_vals<-combined_density_vals
    #density_vals<-combined_density_vals
  }
  # Extract bandwidth (assuming Gaussian kernel)

  h <- kde_result$h
}
# KDE-estimated mean (unbiased, no correction needed)
#kde_mean <- sum(x_vals * density_vals) * diff(kde_result$eval.points[1:2])
# Compute midpoints and bin widths for numerical integration



#midpoints <- (x_vals[-1] + x_vals[-length(x_vals)]) / 2
#mid_density_vals <- (combined_density_vals[-1] + combined_density_vals[-length(combined_density_vals)]) / 2

# Calculate KDE mean robustly (works for constant or varying spacing)
if(dist != "NEWNORMAL") {
  if(reflection == FALSE) {
    #midpoints <- (x_vals[-1] + x_vals[-length(x_vals)]) / 2
    kde_mean <- sum(kde_result$eval.points * kde_result$estimate) * diff(kde_result$eval.points[1:2])
    #grid_x <- kde_result$eval.points[[1]]
    #kde_mode <- grid_x[which.max(kde_result$estimate)]
    # KDE-estimated variance (biased upward due to kernel smoothing)
    kde_variance_raw <- sum((kde_result$eval.points - kde_mean)^2 * kde_result$estimate) * diff(kde_result$eval.points[1:2])

    # Apply theoretical bias-correction for Gaussian kernel
    kde_variance_corrected <- kde_variance_raw - h^2
    # Calculate KDE-based entropy estimate
    density_vals <- predict(kde_result, x = data)
    entropy_kde <- -mean(log(density_vals))
  } else {
    # Use only positive x-values and densities for mean calculation
    positive_indices <- combined_x_vals >= 0
    positive_x_vals <- combined_x_vals[positive_indices]
    positive_density_vals <- combined_density_vals[positive_indices]

    # Compute midpoints for positive densities only
    positive_midpoints <- (positive_x_vals[-1] + positive_x_vals[-length(positive_x_vals)]) / 2
    positive_mid_density_vals <- (positive_density_vals[-1] + positive_density_vals[-length(positive_density_vals)]) / 2

    # Bin widths for positive values only
    positive_dx <- diff(positive_x_vals)

    # Correct KDE mean using only physically meaningful (positive) values
    kde_mean <- sum(positive_midpoints * positive_mid_density_vals * positive_dx)
    # KDE-estimated variance (biased upward due to kernel smoothing)
    kde_variance_raw <- sum((positive_midpoints - kde_mean)^2 * positive_mid_density_vals * positive_dx)

    # Apply theoretical bias-correction for Gaussian kernel
    kde_variance_corrected <- kde_variance_raw - h^2
    midpoints <- (combined_x_vals[-1] + combined_x_vals[-length(combined_x_vals)]) / 2
    mid_density_vals <- (combined_density_vals[-1] + combined_density_vals[-length(combined_density_vals)]) / 2
    combined_dx <- diff(combined_x_vals)
    epsilon <- .Machine$double.eps
    entropy_kde <- -sum(mid_density_vals * log(pmax(mid_density_vals,epsilon)) * combined_dx)
  }
}
#test<-function() {
if(dist == "NORMAL") {
  # Display all results clearly
  entropy_theor <- uni_normal_entropy(theor_variance)
  data_entropy_plugin <- uni_normal_entropy(data_variance)
  kde_entropy_plugin <- uni_normal_entropy(kde_variance_corrected)
  cat("Normal Distribution\n")
  cat("Theoretical Mean:", theor_mean, "\n")
  cat("Theoretical Variance:", theor_variance, "\n")
  cat("Theoretical Skewness:", theor_skew, "\n")
  print(sprintf("Data N: %7d",sample_size))
  cat("Data Mean:", data_mean, "\n")
  cat("Data Variance:", data_variance, "\n")
  cat("Data Minimum:", data_minimum, "\n")
  cat("Data Maximum:", data_maximum, "\n")
  cat("Mean (KDE):", kde_mean, "\n")
  cat("Variance (uncorrected KDE):", kde_variance_raw, "\n")
  cat("Variance (bias-corrected KDE):", kde_variance_corrected, "\n")
  cat("KDE-based Differential Entropy:", entropy_kde, "\n")
  cat("Data plug-in Entropy:", data_entropy_plugin, "\n")
  cat("KDE plug-in Entropy:", kde_entropy_plugin, "\n")
  cat("Theoretical Entropy:", entropy_theor, "\n")
} else if (dist == "HALFNORMAL") {
  # Display all results clearly
  entropy_theor <- uni_half_normal_entropy(theor_variance)
  data_entropy_plugin <- uni_half_normal_entropy(data_variance)
  kde_entropy_plugin <- uni_half_normal_entropy(kde_variance_corrected)
  cat("HALF Normal Distribution\n")
  cat("Theoretical Center",theor_center,"\n")
  cat("Theoretical Mean:", theor_mean, "\n")
  cat("Theoretical Variance:", theor_variance, "\n")
  #cat("Theoretical Skewness:", theor_skew, "\n")
  print(sprintf("Data N: %7d",sample_size))
  cat("Data Mean:", data_mean, "\n")
  cat("Data Variance:", data_variance, "\n")
  cat("Data Skewness:", data_skew, "\n")
  cat("Data Minimum:", data_minimum, "\n")
  cat("Data Maximum:", data_maximum, "\n")
  cat("Mean (KDE):", kde_mean, "\n")
  cat("Variance (uncorrected KDE):", kde_variance_raw, "\n")
  cat("Variance (bias-corrected KDE):", kde_variance_corrected, "\n")
  cat("KDE-based Differential Entropy:", entropy_kde, "\n")
  cat("Data plug-in Entropy:", data_entropy_plugin, "\n")
  cat("KDE plug-in Entropy:", kde_entropy_plugin, "\n")
  cat("Theoretical Entropy:", entropy_theor, "\n")
} else if(dist == "UNIFORM") {
  # a non Gaussian kernel like rectangular would be better
  entropy_theor <- uni_uniform_entropy(theor_minimum,theor_maximum)
  data_entropy_plugin <- uni_uniform_entropy(data_minimum,data_maximum)
  #kernel does not alter minimum and maximum values
  kde_entropy_plugin <- uni_uniform_entropy(data_minimum,data_maximum)
  cat("Uniform Distribution\n")
  cat("Theoretical Mean:", theor_mean, "\n")
  cat("Theoretical Variance:", theor_variance, "\n")
  cat("Theoretical Skewness:", theor_skew, "\n")
  cat("Theoretical Minimum:", theor_minimum, "\n")
  cat("Theoretical Maximum:", theor_maximum, "\n")
  print(sprintf("Data N: %7d",sample_size))
  cat("Data Mean:", data_mean, "\n")
  cat("Data Variance:", data_variance, "\n")
  cat("Data Skewness:", data_skew, "\n")
  cat("Data Minimum:", data_minimum, "\n")
  cat("Data Maximum:", data_maximum, "\n")

  cat("Mean (KDE):", kde_mean, "\n")
  cat("Variance (uncorrected KDE):", kde_variance_raw, "\n")
  cat("Variance (bias-corrected KDE):", kde_variance_corrected, "\n")
  cat("KDE-based Differential Entropy:", entropy_kde, "\n")
  cat("Data plug-in Entropy:", data_entropy_plugin, "\n")
  cat("KDE plug-in Entropy:", kde_entropy_plugin, "\n")
  cat("Theoretical Entropy:", entropy_theor, "\n")
} else if(dist == "TRIANGULAR") {
  # a non Gaussian kernel like rectangular would be better
  entropy_theoretical <- uni_triangle_entropy(theor_minimum,theor_maximum)
  data_entropy_plugin <- uni_triangle_entropy(data_minimum,data_maximum)
  #kernel does not alter minimum and maximum values
  kde_entropy_plugin <- uni_triangle_entropy(data_minimum,data_maximum)
  cat("Triangular Distribution\n")
  cat("Theoretical Mean:", theor_mean, "\n")
  cat("Theoretical Variance:", theor_variance, "\n")
  cat("Theoretical Skewness:", theor_skew, "\n")
  cat("Theoretical Minimum:", theor_minimum, "\n")
  cat("Theoretical Maximum:", theor_maximum, "\n")
  cat("Theoretical Peak:", theor_peak, "\n")
  print(sprintf("Data N: %7d",sample_size))
  cat("Data Mean:", data_mean, "\n")
  cat("Data Variance:", data_variance, "\n")
  cat("Data Skewness:", data_skew, "\n")
  cat("Data Minimum:", data_minimum, "\n")
  cat("Data Maximum:", data_maximum, "\n")

  cat("Mean (KDE):", kde_mean, "\n")
  cat("Variance (uncorrected KDE):", kde_variance_raw, "\n")
  cat("Variance (bias-corrected KDE):", kde_variance_corrected, "\n")
  cat("KDE-based Differential Entropy:", entropy_kde, "\n")
  cat("Data plug-in Entropy:", data_entropy_plugin, "\n")
  cat("KDE plug-in Entropy:", kde_entropy_plugin, "\n")
  cat("Theoretical Entropy:", entropy_theoretical, "\n")
} else if(dist == "EXPONENTIAL") {
  # Display all results clearly
  entropy_theor <- uni_exponential_entropy(theor_lambda)
  data_entropy_plugin <- uni_exponential_entropy(1/data_mean)
  kde_entropy_plugin <- uni_exponential_entropy(1/kde_mean)
  cat("Exponential Distribution\n")
  cat("Theoretical Lambda:", theor_lambda, "\n")
  cat("Theoretical Mean:", theor_mean, "\n")
  cat("Theoretical Variance:", theor_variance, "\n")
  #cat("Theoretical Skewness:", theor_skew, "\n")
  print(sprintf("Data N: %7d",sample_size))
  cat("Data Mean:", data_mean, "\n")
  cat("Data Variance:", data_variance, "\n")
  cat("Data Skewness:", data_skew, "\n")
  cat("Data Minimum:", data_minimum, "\n")
  cat("Data Maximum:", data_maximum, "\n")
  cat("Mean (KDE):", kde_mean, "\n")
  cat("Variance (uncorrected KDE):", kde_variance_raw, "\n")
  cat("Variance (bias-corrected KDE):", kde_variance_corrected, "\n")
  cat("KDE-based Differential Entropy:", entropy_kde, "\n")
  cat("Data plug-in Entropy:", data_entropy_plugin, "\n")
  cat("KDE plug-in Entropy:", kde_entropy_plugin, "\n")
  cat("Theoretical Entropy:", entropy_theor, "\n")
} else if(dist == "LAPLACE") {
  # Display all results clearly
  data_scale<-sqrt(data_variance/2)
  kde_scale_corrected<-sqrt(kde_variance_corrected/2)
  entropy_theor <- uni_laplace_entropy(theor_scale)
  data_entropy_plugin <- uni_laplace_entropy(data_scale)
  kde_entropy_plugin <- uni_laplace_entropy(kde_scale_corrected)
  cat("Laplace Distribution\n")
  cat("Theoretical Mean/Location:", theor_mean, "\n")
  cat("Theoretical Variance:", theor_variance, "\n")
  cat("Theoretical Scale:", theor_scale, "\n")
  #cat("Theoretical Skewness:", theor_skew, "\n")
  print(sprintf("Data N: %7d",sample_size))
  cat("Data Mean/Location:", data_mean, "\n")
  cat("Data Variance:", data_variance, "\n")
  cat("Data Skewness:", data_skew, "\n")
  cat("Data Scale:", data_scale, "\n")
  cat("Data Minimum:", data_minimum, "\n")
  cat("Data Maximum:", data_maximum, "\n")
  cat("Mean/location (KDE):", kde_mean, "\n")
  cat("Variance (uncorrected KDE):", kde_variance_raw, "\n")
  cat("Variance (bias-corrected KDE):", kde_variance_corrected, "\n")
  cat("Scale (bias-corrected KDE):", kde_scale_corrected, "\n")
  cat("KDE-based Differential Entropy:", entropy_kde, "\n")
  cat("Data plug-in Entropy:", data_entropy_plugin, "\n")
  cat("KDE plug-in Entropy:", kde_entropy_plugin, "\n")
  cat("Theoretical Entropy:", entropy_theor, "\n")
} else if(dist == "LOGISTIC") {
  # Display all results clearly
  data_scale<-sqrt(3*data_variance/pi^2)
  kde_scale_corrected<-sqrt(3*kde_variance_corrected/pi^2)
  entropy_theoretical <- uni_logistic_entropy(theor_scale)
  data_entropy_plugin <- uni_logistic_entropy(data_scale)
  kde_entropy_plugin <- uni_logistic_entropy(kde_scale_corrected)
  cat("Logistic Distribution\n")
  cat("Theoretical Mean/Location:", theor_mean, "\n")
  cat("Theoretical Variance:", theor_variance, "\n")
  cat("Theoretical Scale:", theor_scale, "\n")
  #cat("Theoretical Skewness:", theor_skew, "\n")
  print(sprintf("Data N: %7d",sample_size))
  cat("Data Mean/Location:", data_mean, "\n")
  cat("Data Variance:", data_variance, "\n")
  cat("Data Skewness:", data_skew, "\n")
  cat("Data Scale:", data_scale, "\n")
  cat("Data Minimum:", data_minimum, "\n")
  cat("Data Maximum:", data_maximum, "\n")
  cat("Mean/location (KDE):", kde_mean, "\n")
  cat("Variance (uncorrected KDE):", kde_variance_raw, "\n")
  cat("Variance (bias-corrected KDE):", kde_variance_corrected, "\n")
  cat("Scale (bias-corrected KDE):", kde_scale_corrected, "\n")
  cat("KDE-based Differential Entropy:", entropy_kde, "\n")
  cat("Data plug-in Entropy:", data_entropy_plugin, "\n")
  cat("KDE plug-in Entropy:", kde_entropy_plugin, "\n")
  cat("Theoretical Entropy:", entropy_theor, "\n")
} else if(dist == "SKEWNORMAL") {
  fit <- selm(data ~ 1, family = "SN")
  mle<-coef(fit, "DP")  # ξ, ω, α  assumes data is skew normal
  mle_location<-mle[1]
  mle_scale<-mle[2]
  mle_shape<-mle[3]
  mle_delta<-mle_shape/sqrt(1+mle_shape^2) #intermediate value
  mle_mean<-mle_location + mle_scale*mle_delta*sqrt(2/pi)
  mle_variance<-mle_scale^2*(1-2*mle_scale^2/pi)

  # Initial guess: (xi, omega, alpha)
  data_init <- c(data_mean, sqrt(data_variance), 1,data_mean,data_variance,data_skew)

  # Constrain omega > 0
  data_res <- optim(par = data_init, fn = obj_fn_skew_normal, method = "L-BFGS-B",
                    lower = c(-Inf, 1e-6, -10), upper = c(Inf, Inf, 10))

  # Extract estimated parameters
  data_location <- data_res$par[1]
  data_scale <- data_res$par[2]
  data_shape <- data_res$par[3]

  # Display all results clearly
  entropy_theor <- uni_skew_normal_entropy(theor_location,theor_scale,theor_shape)
  entropy_mle <- uni_skew_normal_entropy(mle_location,mle_scale,mle_shape)
  #data_entropy_plugin <- uni_skew_normal_entropy(data_location,data_scale,data_shape)
  #kde_entropy_plugin <- uni_skew_normal_entropy(kde_location,kde_scale,kde_shape)
  cat("Skew Normal Distribution\n")
  cat("Theoretical Location:", theor_location, "\n")
  cat("Theoretical Scale:", theor_scale, "\n")
  cat("Theoretical Shape:", theor_shape, "\n")
  cat("Theoretical Mean:", theor_mean, "\n")
  cat("Theoretical Variance:", theor_variance, "\n")
  cat("Theoretical Skewness:", theor_skew, "\n")
  print(sprintf("Data N: %7d",sample_size))
  cat("Data Location:", data_location, "\n")
  cat("Data Scale:", data_scale, "\n")
  cat("Data Shape:", data_shape, "\n")
  cat("Data Mean:", data_mean, "\n")
  cat("Data Variance:", data_variance, "\n")
  cat("Data Skewness:", data_skew, "\n")
  cat("Data Minimum:", data_minimum, "\n")
  cat("Data Maximum:", data_maximum, "\n")
  cat("MLE Location:", mle_location, "\n")
  cat("MLE Scale:", mle_scale, "\n")
  cat("MLE Shape:", mle_shape, "\n")
  cat("MLE Mean:", mle_mean, "\n")
  cat("MLE Variance:", mle_variance, "\n")
  cat("Mean (KDE):", kde_mean, "\n")
  cat("Variance (uncorrected KDE):", kde_variance_raw, "\n")
  cat("Variance (bias-corrected KDE):", kde_variance_corrected, "\n")
  # cat("KDE-based Differential Entropy:", entropy_kde, "\n")
  # cat("Data plug-in Entropy:", data_entropy_plugin, "\n")
  # cat("KDE plug-in Entropy:", kde_entropy_plugin, "\n")
  cat("Theoretical Entropy:", entropy_theor, "\n")
  cat("MLE Entropy:", entropy_mle, "\n")
}
#}

#test()
