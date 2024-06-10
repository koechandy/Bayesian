if (!requireNamespace("rootSolve", quietly = TRUE)) {
  install.packages("rootSolve")
}
library(rootSolve)


hpdi_f <- function(priorPars, w, mu, tau, n, alpha_level = 0.80) {
  # Function for generating posterior distributions
  fx <- function(x, priorPars, w, mu, tau, n) {
    if (!is.matrix(priorPars) || !is.vector(w)) {
      stop("priorPars must be a matrix and w must be a vector")
    }
    
    if (ncol(priorPars) != length(w)) {
      stop("Number of columns in priorPars must match the length of w")
    }
    
    post.var <- 1 / (1 / priorPars[2, ]^2 + n / tau^2)
    
    gamma <- (tau^2 / n) / ((tau^2 / n) + priorPars[2, ]^2)
    
    post.mean <- gamma * priorPars[1, ] + (1 - gamma) * mu
    
    tauPred <- sqrt(priorPars[2,]^2 + tau^2 / n)
    
    log_marginals <- log(w) + dnorm(mu, priorPars[1,], tauPred, log = TRUE)
    
    margT <- sum(exp(log_marginals))
    
    marginals <- exp(log_marginals)
    
    post.w <- exp(log_marginals - log(margT))
    
    # Posterior densities for each component i aggregated
    posterior_vals <- 0
    for (i in 1:length(post.w)) {
      posterior_vals <- posterior_vals + post.w[i] * dnorm(x, mean = post.mean[i], sd = sqrt(post.var[i]))
    }
    return(list(f = posterior_vals, post.w = post.w, post.m = post.mean, post.s = sqrt(post.var)))
  }
  
  x.lo<-floor(min(c(priorPars[1,],mu))) + ceiling(max(c(priorPars[2,],tau)))
  x.up<-ceiling(max(c(priorPars[1,],mu))) + ceiling(max(c(priorPars[2,],tau)))
  x <- seq(-x.lo, x.up, by = 0.001)
  y <- seq(max(fx(x, priorPars, w, mu, tau, n)$f), 0, -0.001)
  
  
  post.w <- fx(x, priorPars, w, mu, tau, n)$post.w
  post.m <- fx(x, priorPars, w, mu, tau, n)$post.m
  post.s <- fx(x, priorPars, w, mu, tau, n)$post.s
  postPars <- matrix(c(post.m, post.s), nrow = 2, byrow = TRUE)
  
  # Function to calculate alpha
  my_alpha <- function(y, priorPars, w, mu, tau, n, postPars, alpha_level) {
    roots <- rootSolve::uniroot.all(function(x) fx(x, priorPars, w, mu, tau, n)$f - y, lower = min(x), upper = max(x), tol = 1e-10)
    
    num_roots <- length(roots)
    alpha <- 0
    alpha_sum <- 0
    
    if (num_roots == 2) {
      alpha <- sum(sapply(1:length(w), function(j) {
        post.w <- fx(roots, priorPars, w, mu, tau, n)$post.w
        post.w[j] * pnorm(roots[2], postPars[1, j], postPars[2, j]) - post.w[j] * pnorm(roots[1], postPars[1, j], postPars[2, j])
      }))
    } else if (num_roots == 4) {
      for (k in seq(1, num_roots, by = 2)) {
        z <- (k + 1) / 2
        post.w <- fx(roots, priorPars, w, mu, tau, n)$post.w
        alpha_sum <- alpha_sum + (post.w[z] * pnorm(roots[k + 1], postPars[1, z], postPars[2, z]) -
                                    post.w[z] * pnorm(roots[k], postPars[1, z], postPars[2, z]))
      }
      alpha <- alpha_sum
    } else {
      alpha <- NA  # when the number of roots is neither 2 nor 4
    }
    
    return(alpha)
  }
  
  # Objective function for optimization
  objectivefx <- function(y, ...) { 
    alpha <- my_alpha(y, ...)
    abs(alpha - alpha_level)
  }
  
  opt.res <- optim(0, objectivefx, method = "Brent", lower = min(y), upper = max(y), 
                   priorPars = priorPars, w = w, mu = mu, tau = tau, n = n, 
                   postPars = postPars, alpha_level = alpha_level)
  
  #--------------- roots corresponding to the minimum absolute difference--------
  hpd <- rootSolve::uniroot.all(function(x) fx(x, priorPars, w, mu, tau, n)$f - opt.res$par, lower = min(x), upper = max(x),tol = 1e-10)
  
  ###---------combined mixture prior distribution-------------------------------
  cmp <- 0
  for (i in 1:length(w)) {
    cmp <- cmp + w[i] * dnorm(x, priorPars[1, i], priorPars[2, i])
  }
  
  posterior_density <- fx(x, priorPars, w, mu, tau, n)$f
  
  #--------------------Determine y-axis upper limit---------------------------
  up.lim <- max(c(cmp, posterior_density))
  
  x.ll<-floor(range(density(post.m)$x)[1])-2
  x.ul<-ceiling(range(density(post.m)$x)[2])+2
  
  #--------------Plot the posterior distribution and the HPD interval---------
  plot(x, cmp, type = "l", col = "darkorange", lty = 1, lwd = 2,
       xlab = expression(theta), ylab = "Density",xlim=c(x.ll,x.ul), ylim = c(0, up.lim))
  
  colors <- c("blue", "purple", "darkgreen", "darkcyan", "magenta")
  
  for (i in 1:length(post.w)) {
    lines(x, post.w[i]*dnorm(x, postPars[1, i], postPars[2, i]), type = "l",
          col = colors[i], lty = 2, lwd = 2)
  }
  
  #--------------------Posterior mixture density--------------------------------
  lines(x, posterior_density, type = "l", lwd = 2, lty = 1, col = "black")
  legend_labels <- c(paste0("Posterior Component ", 1:length(w), ": N(", round(postPars[1, ], 3), ",", round(postPars[2, ], 3), ")"),
                     "Mixture prior density", "Posterior mixture density", "HPD Interval")
  
  legend("topleft", lty = c(rep(2, length(w)), 1, 1),
         legend = legend_labels, lwd = 2,
         col = c(colors[1:length(w)], "darkorange", "black", "red"), bty = "n")
  
  #--------Adding vertical and horizontal lines to the HPD interval-------------
  hpd_segments <- function(hpd) {
    for (i in seq(1, length(hpd), by = 2)) {
      segments(hpd[i], 0, hpd[i], fx(hpd[i], priorPars, w, mu, tau, n)$f, col = "red", lwd = 2, lty = 2)
      segments(hpd[i + 1], 0, hpd[i + 1], fx(hpd[i + 1], priorPars, w, mu, tau, n)$f, col = "red", lwd = 2, lty = 2)
      segments(hpd[i], fx(hpd[i], priorPars, w, mu, tau, n)$f, hpd[i + 1], fx(hpd[i], priorPars, w, mu, tau, n)$f, col = "red", lwd = 2, lty = 2)
    }
  }
  
  #--------------HPD interval-----------------------------------------
  if (length(hpd) %in% c(2, 4)) {
    hpd_segments(hpd)
  }
  
  
  # Return the roots equivalent to alpha level
  hpd_string <-if (length(hpd) == 2) {
    cat(alpha_level * 100, "% HPDI is:[",round(hpd[1],4),",",round(hpd[2],4),"]\n")
  } else if (length(hpd) == 4) {
    cat(alpha_level * 100, "% HPDI is:[",round(hpd[1],4),",",round(hpd[2],4),"] and [",round(hpd[3],4),",",round(hpd[4],4),"]\n")
  } else {
    cat("No HPDI for", alpha_level * 100, "%\n")
  }
  
  cat(hpd_string)
  
  return(hpd_string = hpd_string)
}

mu.vals <- seq(-2, 2, by = 0.1)
priorPars_list <- lapply(mu.vals, function(prior) {
  var.vals <- c((1 / sqrt(50)), 10)
  matrix(c(prior, prior,var.vals[1],var.vals[2]), nrow = 2, byrow = TRUE)
})

seq_w <- seq(0, 1, by = 0.01)
w_m <- lapply(seq_w, function(val) c(val, 1 - val))

mu <- 0.78
tau <- 1
n <- 20
alpha_level_seq <- seq(0.80, 0.95,by=0.05)


param_grid <- expand.grid(
  priorPars = I(priorPars_list),
  w = I(w_m),
  mu = mu,
  tau = tau,
  n = n,
  alpha_level = alpha_level_seq
)

hpd <- apply(param_grid, 1, function(param) {
  hpdi_f(param$priorPars, param$w, param$mu, param$tau, param$n, param$alpha_level)
})

hpd
