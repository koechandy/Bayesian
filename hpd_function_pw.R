if (!requireNamespace("rootSolve", quietly = TRUE)) {
  install.packages("rootSolve")
}
library(rootSolve)

# Updated function for generating posterior distributions
fx <- function(x, priorPars, w, mu, tau, n) {
  if (!is.matrix(priorPars) || !is.vector(w)) {
    stop("priorPars must be a matrix and w must be a vector")
  }
  
  if (ncol(priorPars) != length(w)) {
    stop("Number of columns in priorPars must match the length of w")
  }
  
  post.var <- 1 / (1 / priorPars[2, ]^2 + n / tau^2)
  post.mean <- matrix(NA, nrow = length(mu), ncol = length(w))
  
  gamma <- (tau^2 / n) / ((tau^2 / n) + priorPars[2, ]^2)
  for (i in 1:length(w)) {
    post.mean[, i] <- gamma[i] * priorPars[1, i] + (1 - gamma[i]) * mu
  }
  
  tauPred <- matrix(NA, length(mu), length(w))
  for (i in 1:length(w)) {
    tauPred[, i] <- sqrt(priorPars[2, i]^2 + tau^2 / n)
  }
  
  margT <- 0
  log_marginals <- matrix(NA, length(mu), length(w))
  for (i in 1:length(w)) {
    log_marginals[, i] <- log(w[i]) + dnorm(mu, priorPars[1, i], tauPred[, i], log = TRUE)
    margT <- margT + exp(log_marginals[, i])
  }
  
  marginals <- exp(log_marginals)
  
  post.w <- matrix(NA, length(mu), length(w))
  for (i in 1:length(w)) {
    post.w[, i] <- exp(log_marginals[, i] - log(margT))
  }
  
  postmean <- 0
  for (i in 1:length(w)) {
    postmean <- postmean + post.w[, i] * post.mean[, i]
  }
  
  postvar <- 0
  for (i in 1:length(w)) {
    postvar <- postvar + post.w[, i] * (post.var[i] + post.mean[, i]^2)
  }
  postvar <- postvar - postmean^2
  
  ##---posterior densities for each component i aggregated------------------
  posterior_vals <- 0
  for (i in 1:length(w)) {
    posterior_vals <- posterior_vals + post.w[, i]*dnorm(x, mean = post.mean[, i], sd = sqrt(post.var[i]))
  }
  return(list(f = posterior_vals, post.w=post.w))
}

##----------- Declaring the parameters ------------------------------#
x <- seq(-2, 2, by = 0.01)
priorPars <- matrix(c(0,0,(1/sqrt(50)),10),nrow =2,byrow = T)
w <- c(0.5,0.5)
mu <- 0.78
tau <- 1
n <- 20

alpha_level <- 0.60

# Call the function with the declared parameters
y <- seq(max(fx(x, priorPars, w, mu, tau, n)$f), 0, -0.001)

rt_list <- list()
for (i in 1:length(y)) {
  roots <- rootSolve::uniroot.all(function(x) fx(x, priorPars, w, mu, tau, n)$f - y[i], lower = min(x), upper = max(x))
  
  rt_list[[i]] <- list(y.val = y[i], roots = roots)
  
  num_roots <- length(rt_list[[i]]$roots)
  alpha <- 0
  alpha_sum <- 0
  
  if (num_roots == 2 || num_roots == 4) {
    for (j in seq(1, num_roots, by = 2)) {
      for (k in 1:length(w)) {

        alpha_sum <- alpha_sum + sum(fx(x, priorPars, w, mu, tau, n)$post.w *
                                       (pnorm(rt_list[[i]]$roots[j + 1], priorPars[1, k], priorPars[2, k]) - 
                                          pnorm(rt_list[[i]]$roots[j], priorPars[1, k], priorPars[2, k])))
      }
    }
    alpha<- alpha_sum
  } else {
    alpha<- NA  # when the number of roots is neither 2 nor 4
  }
  
  rt_list[[i]] <- list(y.val = y[i], roots = roots, alpha = alpha)
}

alphal <- do.call(c, lapply(rt_list, function(y) y$alpha))

hpd <- rt_list[[which.min(abs(alphal - alpha_level))]]$roots

# Combined mixture prior
cmp <- 0
for (i in 1:length(w)) {
  cmp <- cmp + w[i] * dnorm(x, priorPars[1, i], priorPars[2, i])
}

posterior_density <- fx(x, priorPars, w, mu, tau, n)$f

# Determine ylim upper limit
up.lim <- max(c(cmp, posterior_density))

# Plot the posterior distribution and the HPD interval
plot(x, cmp, type = "l", col = "darkorange", lty = 1, lwd = 2,
     main = "Mixture prior components and posterior density plots",
     xlab = expression(theta), ylab = "Density",ylim=c(0,up.lim))

colors <- c("blue", "purple", "darkgreen", "darkcyan", "magenta")
for (i in 1:length(w)) {
  lines(x, dnorm(x, priorPars[1, i], priorPars[2, i]), type = "l", 
        col = colors[i], lty = 2, lwd = 2)
}

# Posterior mixture density
lines(x, posterior_density, type = "l", lwd = 2, lty = 1, col = "black")

# Adding vertical and horizontal lines to the HPD interval
hpd_segments <- function(hpd) {
  for (i in seq(1, length(hpd), by = 2)) {
    segments(hpd[i], 0, hpd[i], fx(hpd[i], priorPars, w, mu, tau, n)$f, col = "red", lwd = 2, lty = 2)
    segments(hpd[i + 1], 0, hpd[i + 1], fx(hpd[i + 1], priorPars, w, mu, tau, n)$f, col = "red", lwd = 2, lty = 2)
    segments(hpd[i], fx(hpd[i], priorPars, w, mu, tau, n)$f, hpd[i + 1], fx(hpd[i], priorPars, w, mu, tau, n)$f, col = "red", lwd = 2, lty = 2)
  }
}

# HPD interval
if (length(hpd) %in% c(2, 4)) {
  hpd_segments(hpd)
}

legend_labels <- c(paste0("Component ", 1:length(w), ": N(", priorPars[1, ], ",", round(1/priorPars[2, ], 2), ")"),
                   "Mixture prior density", "Posterior mixture density", "HPD Interval")

legend("topleft", lty = c(rep(2, length(w)), 1, 1), 
       legend = legend_labels, lwd = 2,
       col = c(colors[1:length(w)], "darkorange", "black", "red"), bty = "n")

# Return the roots equivalent to alpha level
if (length(hpd) == 2) {
  cat(alpha_level * 100, "% HPDI is:[", hpd[1], ",", hpd[2], "]\n")
} else if (length(hpd) == 4) {
  cat(alpha_level * 100, "% HPDI is:[", hpd[1], ",", hpd[2], "] and [", hpd[3], ",", hpd[4], "]\n")
} else {
  cat("No HPDI for", alpha_level * 100, "%\n")
}


