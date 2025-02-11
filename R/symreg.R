# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function(digit) {
  pi <- 3.141529
  print(paste("The ", digit, " of pi is: ", pi, sep = ""))
}

generate_distribution_metrics <- function(num_distributions = 1000,
                                          samples_per_dist = 500,
                                          dist_type = c("log-normal", "normal", "uniform", "exponential")) {
  library(moments)  # For skewness calculation

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
