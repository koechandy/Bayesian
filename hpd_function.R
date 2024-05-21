#Code snippet one trial
  f <- function(x) {
    priorPars <- matrix(c(7, 9, 11, 1/sqrt(1), 1/sqrt(0.25), sqrt(1)), nrow = 2, byrow = TRUE)
    w <- c(0.62, 0.03, 0.35)
    
    f2 <- w[1] * dnorm(x, priorPars[1, 1], priorPars[2, 1]) +
      w[2] * dnorm(x, priorPars[1, 2], priorPars[2, 2]) +
      w[3] * dnorm(x, priorPars[1, 3], priorPars[2, 3])
    return(f2)
  }
  

  x <- seq(-5, 20, by = 0.01)
  
  y <- seq(1, 0, -0.001)
  priorPars <- matrix(c(7, 9, 11, 1/sqrt(1), 1/sqrt(0.25), sqrt(1)), nrow = 2, byrow = TRUE)
  w <- c(0.62, 0.03, 0.35)

  rt_list <- list()
  alpha <- 0
  alpha_level <- 0.90
  
  for (i in 1:length(y)) {
    roots <- rootSolve::uniroot.all(function(x) f(x) - y[i], lower = min(x), upper = max(x))
        
    
    rt_list[[i]] <- list(y.val = y[i], roots = roots)
    
    #Calculate alpha based on the number of roots
    num_roots <- length(rt_list[[i]]$roots)
    
    
    if (num_roots == 2) {
      
      alpha[i] <- w[1]*pnorm(rt_list[[i]]$roots[2], priorPars[1,1], priorPars[2,1]) +  
        w[2]*pnorm(rt_list[[i]]$roots[2], priorPars[1,2], priorPars[2,2]) +  
        w[3]*pnorm(rt_list[[i]]$roots[2], priorPars[1,3], priorPars[2,3]) -  
        w[1]*pnorm(rt_list[[i]]$roots[1], priorPars[1,1], priorPars[2,1]) - 
        w[2]*pnorm(rt_list[[i]]$roots[1], priorPars[1,2], priorPars[2,2]) - 
        w[3]*pnorm(rt_list[[i]]$roots[1], priorPars[1,3], priorPars[2,3])
      
    } else if (num_roots == 4) {
      alpha[i] <- w[1]*pnorm(rt_list[[i]]$roots[2], priorPars[1,1], priorPars[2,1]) +  
        w[2]*pnorm(rt_list[[i]]$roots[2], priorPars[1,2], priorPars[2,2]) +  
        w[3]*pnorm(rt_list[[i]]$roots[2], priorPars[1,3], priorPars[2,3]) -  
        w[1]*pnorm(rt_list[[i]]$roots[1], priorPars[1,1], priorPars[2,1]) - 
        w[2]*pnorm(rt_list[[i]]$roots[1], priorPars[1,2], priorPars[2,2]) - 
        w[3]*pnorm(rt_list[[i]]$roots[1], priorPars[1,3], priorPars[2,3]) +
        w[1]*pnorm(rt_list[[i]]$roots[4], priorPars[1,1], priorPars[2,1]) +  
        w[2]*pnorm(rt_list[[i]]$roots[4], priorPars[1,2], priorPars[2,2]) +  
        w[3]*pnorm(rt_list[[i]]$roots[4], priorPars[1,3], priorPars[2,3]) -  
        w[1]*pnorm(rt_list[[i]]$roots[3], priorPars[1,1], priorPars[2,1]) - 
        w[2]*pnorm(rt_list[[i]]$roots[3], priorPars[1,2], priorPars[2,2]) - 
        w[3]*pnorm(rt_list[[i]]$roots[3], priorPars[1,3], priorPars[2,3])
      
    } else {
      alpha[i] <- NA  # when number of roots is neither 2 or 4
    }
    
    # storing the corresponding alphas in a list
    rt_list[[i]] <- list(y.val = y[i], roots = roots,alpha=alpha[i])
    
}

  alphal = do.call(c,lapply(rt_list, function(y) y$alpha))
  
  hpd <- rt_list[[which.min(abs(alphal - alpha_level))]]$roots
  
  # return the roots
  if (length(hpd) == 2) {
    cat(alpha_level*100,"% HPD is:[",hpd[1],",",hpd[2],"]\n")
  } else if (length(hpd) == 4) {
    cat(alpha_level*100,"% HPD Interval is:[",hpd[1],",",hpd[2],"] and [",hpd[3],",",hpd[4], "]\n")
  } else {
    cat("No HPD Interval for ",alpha_level*100,"%\n")
  }  
  