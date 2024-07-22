if (!requireNamespace("shiny", quietly = TRUE)) {
  install.packages("shiny")
}
if (!requireNamespace("bslib", quietly = TRUE)) {
  install.packages("bslib")
}
if (!requireNamespace("rootSolve", quietly = TRUE)) {
  install.packages("rootSolve")
}

if (!requireNamespace("DT", quietly = TRUE)) {
  install.packages("DT")
}

library(shiny)
library(bslib)
library(rootSolve)
library(DT)

# Define the hpdi function
options(digits = 5)
rhpdi <- function(priorPars, w, mu, tau, n, alpha_level = 0.85) {
  options(scipen = 999)
  
  
  if (!is.matrix(priorPars) || !is.vector(w)) {
    stop("priorPars must be a matrix and w must be a vector")
  }
  if (ncol(priorPars) != length(w)) {
    stop("Number of columns in priorPars must match the length of w")
  }
  post.var <- 1 / (1 / priorPars[2, ]^2 + n / tau^2)
  
  gamma <- (tau^2 / n) / ((tau^2 / n) + priorPars[2, ]^2)
  
  post.m <- gamma * priorPars[1, ] + (1 - gamma) * mu
  
  tauPred <- sqrt(priorPars[2,]^2 + tau^2 / n)
  
  log_marginals <- log(w) + dnorm(mu, priorPars[1,], tauPred, log = TRUE)
  
  margT <- sum(exp(log_marginals))
  
  marginals <- exp(log_marginals)
  
  post.w <- exp(log_marginals - log(margT))
  
  fx <- function(x){
    
    posterior_vals <- 0
    for (i in 1:length(post.w)) {
      posterior_vals <- posterior_vals + 
        post.w[i] * dnorm(x, mean = post.m[i], sd = sqrt(post.var[i]))
    }
    
    return(f = posterior_vals)
  }
  
  # min_pr <- sum(w * (priorPars[1, ] - 3 * priorPars[2, ]))
  # max_pr <- sum(w * (priorPars[1, ] + 3 * priorPars[2, ]))
  
  min_post <- sum(post.w * (post.m - 3 * sqrt(post.var)))
  max_post <- sum(post.w * (post.m + 3 * sqrt(post.var)))
  
  min_lik <- mu - 3*tau
  max_lik <- mu + 3*tau
  
  x <- seq(min(c(min_post,min_lik)),max(c(max_post,max_lik)),by = 0.001)
  y <- seq(max(fx(x)), 0, by = -0.001)
  post.s <- sqrt(post.var)
  postPars <- matrix(c(post.m, post.s), nrow = 2, byrow = TRUE)
  
  my_alpha <- function(y) {
    
    roots <- rootSolve::uniroot.all(function(x) fx(x) - y, lower = min(x), 
                                    upper = max(x),tol = 1e-10,maxiter = 10, 
                                    trace = 0, n = 1000)
    
    num_roots <- length(roots)
    alpha <- 0
    alpha <- sum(sapply(1:length(w), function(j) {
      sum(diff(post.w[j] * pnorm(roots, postPars[1, j], 
                                 postPars[2, j]))[seq(1,length(roots),2)])
    }))
    return(c(alpha,roots))
  }
  
  objectivefx <- function(y,...) {
    alpha <- my_alpha(y)[1]
    abs(alpha - alpha_level)
  }
  
  opt.res <- optim(0, objectivefx, method = "Brent", lower = min(y), 
                   upper = max(y),priorPars = priorPars, w = w, mu = mu, 
                   tau = tau, n = n,postPars = postPars, 
                   alpha_level = alpha_level)
  
  hpd <- my_alpha(opt.res$par)[-1]
  
  
  cmp <- 0
  
  for (i in 1:length(w)) {
    cmp <- cmp + w[i] * dnorm(x, priorPars[1, i], priorPars[2, i])
  }
  
  posterior_density <- fx(x)
  lik<- dnorm(x, mu,tau)
  up.lim <- max(c(cmp, posterior_density,lik))
  
  
  x.min <- ifelse(min(x) > 0, floor(min(x)), ceiling(min(x)))
  
  x.max <- ifelse(max(x) > 0, ceiling(max(x)), floor(max(x)))
  
  #---plotting of the posterior and respective likelihood and mixture priors---
  plot(x, cmp, type = "l", col = "darkorange", lty = 2, lwd = 1,
       xlab = expression(theta), ylab = "Density", 
       xlim = c(x.min, x.max), ylim = c(0, up.lim),
       cex.lab = 1.5, font.lab = 2, col.lab = "black", cex.axis = 1.5, font.axis = 2, col.axis = "black")
  
  
  axis(1, col.axis = "black", cex.axis = 1.5, font.axis = 2)
  axis(2, col.axis = "black", cex.axis = 1.5, font.axis = 2)
  
  #add the likelihood
  lines(x, lik, type = "l", lwd = 1, lty = 2, col = "darkblue")
  
  lines(x, posterior_density, type = "l", lwd = 2, lty = 1, col = "black")
  
  legend_labels <- c("prior mixture density","Likelihood", 
                     "Posterior mixture density","HPD Interval")
  legend("topleft", lty = c(2,2,1,2),
         legend = legend_labels, lwd = c(1,1,2,1),
         col = c("darkorange","darkblue", "black","red"), bty = "n", cex = 1.5, text.font = 2)
  
  
  hpd_segments <- function(hpd) {
    
    for (i in seq(1, length(hpd), by = 2)) {
      
      segments(hpd[i], 0, hpd[i], fx(hpd[i]), col = "red", lwd = 1, lty = 2)
      segments(hpd[i + 1], 0, hpd[i + 1], fx(hpd[i + 1]), col = "red", 
               lwd = 1, lty = 2)
      segments(hpd[i], fx(hpd[i]), hpd[i + 1], fx(hpd[i + 1]), col = "red", 
               lwd = 1, lty = 2)
    }
  }
  
  if (length(hpd) %in% c(2, 4)) {
    hpd_segments(hpd)
  }
  
  hpd_string <- if (length(hpd) == 2) {
    paste0(alpha_level * 100, "% HPDI is:[", round(hpd[1], 4), ",",
           round(hpd[2], 4), "]")
  } else if (length(hpd) == 4) {
    paste0(alpha_level * 100, "% HPDI is:[", round(hpd[1], 4), ",", 
           round(hpd[2], 4), "] and [", round(hpd[3], 4), ",", 
           round(hpd[4], 4), "]")
  } else if (length(hpd) == 6) {
    paste0(alpha_level * 100, "% HPDI is:[", round(hpd[1], 4), ",", 
           round(hpd[2], 4), "] and [", round(hpd[3], 4), ",", 
           round(hpd[4], 4), "] and [", round(hpd[5], 4), ",", 
           round(hpd[6], 4), "]")
  }else {
    paste("No HPDI for", alpha_level * 100, "%")
  }
  return(list(hpd_string,post.w,post.m,post.s))
}

priorPars <- matrix(c(0, 0, (1 / sqrt(50)), 10), nrow = 2, byrow = TRUE)
w <- c(0.5, 0.5)
mu <- 3
tau <- 7
n <- 40


hpdi(priorPars, w, mu, tau, n, alpha_level = 0.85)





