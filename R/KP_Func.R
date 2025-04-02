
rm(list=ls())
# Install and load the necessary package
library(ks)
library(R.utils)
library(MASS)
library(sn)

#dist<-"UNIFORM"
#dist<-"NORMAL"
#dist<-"TRIANGULAR"
#dist<-"EXPONENTIAL" #log transformation and reflection
#dist<-"LAPLACE"
#dist<-"LOGISTIC" #seems to be some problems
#dist<-"HALFNORMAL" #seems to be some problems
dist<-"SKEWNORMAL"

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
  entropy<-log(aScale*sqrt(2*pi*exp(1))) - delta^2
  return(entropy)
}


bi_normal_entropy<-function(aVarX,aVarY,aCorr) {

  return(1+log(2*pi)+1/2*log(aVarX*aVarY*(1-aCorr^2)))
}



sample_size<-5000000

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

# KDE-estimated mean (unbiased, no correction needed)
#kde_mean <- sum(x_vals * density_vals) * diff(kde_result$eval.points[1:2])
# Compute midpoints and bin widths for numerical integration



#midpoints <- (x_vals[-1] + x_vals[-length(x_vals)]) / 2
#mid_density_vals <- (combined_density_vals[-1] + combined_density_vals[-length(combined_density_vals)]) / 2

# Calculate KDE mean robustly (works for constant or varying spacing)
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
