
fx <- function(x, mu, priorPars, w, tau, n) {
  
  # Check that priorPars is a matrix and w is a vector
  if (!is.matrix(priorPars) || !is.vector(w)) {
    stop("priorPars must be a matrix and w must be a vector")
  }
  
  # Check if the number of columns in priorPars matches the length of w
  if (ncol(priorPars) != length(w)) {
    stop("Number of columns in priorPars must match the length of w")
  }
  
  # Calculate post.var
  post.var <- 1 / (1 / priorPars[2, ]^2 + n / tau^2)
  

  post.mean <- matrix(NA, nrow = length(mu), ncol = length(w))
  
  # Calculate gamma and post.mean
  gamma <- (tau^2 / n) / ((tau^2 / n) + priorPars[2, ]^2)
  for (i in seq_along(w)) {
    post.mean[, i] <- gamma[i] * priorPars[1, i] + (1 - gamma[i]) * mu
  }
  

  tauPred <- matrix(NA, nrow = length(mu), ncol = length(w))
  
  # Calculate tauPred
  for (i in seq_along(w)) {
    tauPred[, i] <- sqrt(priorPars[2, i]^2 + tau^2 / n)
  }
  
  # Calculate margT and individual marginals
  log_marginals <- matrix(NA, nrow = length(mu), ncol = length(w))
  for (i in seq_along(w)) {
    log_marginals[, i] <- log(w[i]) + dnorm(mu, priorPars[1, i], tauPred[, i], log = TRUE)
  }
  
  margT <- rowSums(exp(log_marginals))
  marginals <- exp(log_marginals)
  
  # Calculate post.weights
  post.w <- exp(log_marginals - log(margT))
  
  # Calculate postmean
  postmean <- rowSums(post.w * post.mean)
  
  # Calculate postvar
  postvar <- rowSums(post.w * (post.var + post.mean^2)) - postmean^2
  
  return(list(postmean = postmean, postvar = postvar))
# 
#   f2 <- dnorm(x, mean = postmean, sd = 1/sqrt(postvar))
# 
#   return(f2)
}


## Example usage
# Declare the parameters outside the function
x <- seq(-5, 20, by = 0.01)
priorPars <- matrix(c(7, 9, 11, 1/sqrt(1), 1/sqrt(0.25), 1/sqrt(1)), nrow = 2, byrow = TRUE)
w <- c(0.62, 0.03, 0.35)
mu <- 206.3/20
tau <- 1/sqrt(0.4) 
n <- 20

w1 <- c(0.25,0.5,0.25)

# Call the function with the declared parameters
fx(x,mu, priorPars, w1, tau, n)
fx(x,mu, priorPars, w, tau, n)



