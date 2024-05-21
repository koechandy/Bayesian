# updated function for generating posterior distributions
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
  
  # Calculate posterior weights
  post.w <- exp(log_marginals - log(margT))
  
  # Calculate posterior mean
  postmean <- rowSums(post.w * post.mean)
  
  # Calculate posterior variance
  postvar <- rowSums(post.w * (post.var + post.mean^2)) - postmean^2
  
  # return(list(postmean = postmean, postvar = postvar))
  
  f2 <- dnorm(x, mean = postmean, sd = 1/sqrt(postvar))
  
  return(f2)
}


## Example usage
# Declare the parameters outside the function
x <- seq(-5, 30, by = 0.01)
priorPars <- matrix(c(3,9,1/sqrt(.25),sqrt(1)),nrow =2,byrow = T)
w <- c(0.75,0.25)
mu <- 206.3/20
tau <- 1/sqrt(0.4) 
n <- 20
alpha_level <- 0.90

# # Call the function with the declared parameters
# fx(x,mu, priorPars, w, tau, n)


y <- seq(1, 0, -0.001)
rt_list <- list()
alpha <- 0

for (i in 1:length(y)) {
  roots <- rootSolve::uniroot.all(function(x) fx(x, mu, priorPars, w, tau, n) - y[i], lower = min(x), upper = max(x))
  
  
  rt_list[[i]] <- list(y.val = y[i], roots = roots)
  
  num_roots <- length(rt_list[[i]]$roots)
  alpha_sum <- 0
  
  if (num_roots == 2 || num_roots == 4) {
    for (j in seq(1, num_roots, by = 2)) {
      alpha_sum <- alpha_sum + sum(w * (pnorm(rt_list[[i]]$roots[j + 1], priorPars[1,], priorPars[2,]) - pnorm(rt_list[[i]]$roots[j], priorPars[1,], priorPars[2,])))
    }
    alpha[i] <- alpha_sum
  } else {
    alpha[i] <- NA  # when number of roots is neither 2 or 4
  }
  
  # Storing the corresponding alphas in a list
  rt_list[[i]] <- list(y.val = y[i], roots = roots, alpha = alpha[i])
}

alphal <- do.call(c, lapply(rt_list, function(y) y$alpha))

hpd <- rt_list[[which.min(abs(alphal - alpha_level))]]$roots

# return the roots
if (length(hpd) == 2) {
  cat(alpha_level*100,"% HPD is:[",hpd[1],",",hpd[2],"]\n")
} else if (length(hpd) == 4) {
  cat(alpha_level*100,"% HPD Interval is:[",hpd[1],",",hpd[2],"] and [",hpd[3],",",hpd[4], "]\n")
} else {
  cat("No HPD Interval for ",alpha_level*100,"%\n")
}  
