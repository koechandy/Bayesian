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


options(digits = 5)
#defining the hpdi function
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
  
  #plotting of the posterior and respective likelihood and mixture priors---
  plot(x, cmp, type = "l", col = "darkorange", lty = 2, lwd = 2,
       xlab = expression(theta), ylab = "Density", 
       xlim = c(x.min, x.max), ylim = c(0, up.lim),
       cex.lab = 1.5, font.lab = 1.8, col.lab = "black", cex.axis = 1.5, font.axis = 2, col.axis = "black")
  
  
  axis(1, col.axis = "black", cex.axis = 1.5, font.axis = 2)
  axis(2, col.axis = "black", cex.axis = 1.5, font.axis = 2)
  
  #add the likelihood
  lines(x, lik, type = "l", lwd = 2, lty = 2, col = "darkblue")
  
  lines(x, posterior_density, type = "l", lwd = 2, lty = 1, col = "black")
  
  legend_labels <- c("prior mixture density","Likelihood", 
                     "Posterior mixture density","HPD Interval")
  legend("topleft", lty = c(2,2,1,2),
         legend = legend_labels, lwd = c(2,2,2,2),
         col = c("darkorange","darkblue", "black","red"), bty = "n", cex = 1.5, text.font = 1.8)
  
  hpd_segments <- function(hpd) {
    
    for (i in seq(1, length(hpd), by = 2)) {
      
      segments(hpd[i], 0, hpd[i], fx(hpd[i]), col = "red", lwd = 2, lty = 2)
      segments(hpd[i + 1], 0, hpd[i + 1], fx(hpd[i + 1]), col = "red", 
               lwd = 2, lty = 2)
      segments(hpd[i], fx(hpd[i]), hpd[i + 1], fx(hpd[i + 1]), col = "red", 
               lwd = 2, lty = 2)
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



ui <- page_sidebar(
  title = "HPD Interval Computation",
  theme = bslib::bs_theme(bootswatch = "united"),
  
  sidebar = sidebar(
    width = 500,
    tags$head(
      tags$style(HTML("
        .sidebar {
          background-color: #f0f0f0;
          padding: 15px;
          border-radius: 5px;
        }
        .well {
          background-color: #f0f0f0;
          border: 1px solid #ddd;
        }
      "))
    ),
    tags$style(type = "text/css", "#hpd_text {
      color: blue;
      text-align: center;
      font-size: 22px;
      font-weight: bold;
    }"),
    h4("Likelihood parameters"),
    numericInput("mu", HTML("mean (&#956;)"), 3, min = 0, step = 1),
    numericInput("tau", HTML("standard deviation (&#964;)"), 7, min = 0, step = 1),
    numericInput("n", "sample size (n):", 40),
    h4("Prior parameters"),
    selectInput("num_prior_params", "Number of prior parameters:", choices = 1:10, selected = 2),
    uiOutput("dynamic_prior_inputs"),
    numericInput("alpha_level", "alpha level:", 0.85, min = 0, max = 1, step = 0.05),
    actionButton("compute", "Compute", class = "btn-lg btn-success"),
    class = "sidebar"
  ),
  
  navset_card_underline(
    # Panel with intro ----
    nav_panel("About",
              p(tags$span(style = "color: blue;","The WebApp rHPDI is a web application that computes the highest posterior density interval (HPDI) for a given set of informative and vague prior parameters, weights, likelihood parameters and desired alpha level.")),
              h4("Inputs:"),
              tags$ul(
                tags$li(HTML("prior means (&#956;)"), ": Informative and vague prior means."),
                tags$li(HTML("prior standard deviation (&#964;)"), ": respective informative and vague prior standard deviations"),
                tags$li("weights (w): a vector of weights corresponding to the prior components"),
                tags$li(HTML("mean (&#956;)"), ": The mean of the likelihood."),
                tags$li(HTML("tau (&#964;)"), ": The standard deviation of the likelihood."),
                tags$li("n: The sample size from the likelihood."),
                tags$li("alpha level: The desired level of the HPDI, typically set to 0.90.")
              ),
              h4("Outputs:"),
              tags$ul(
                tags$li("plot representing the HPD interval together with the likelihood, prior & posterior mixture densities."),
                tags$li("values representing the Highest Posterior Density interval (HPDI) at desired alpha level."),
                tags$li("posterior weights"),
                tags$li("posterior means."),
                tags$li("posterior variances")
              )
    ),
    # Panel with plot ----
    nav_panel("Visualization",
              plotOutput("plot_hpd", width = "95%", height = "900px"),
              verbatimTextOutput("hpd_text")
    ),
    
    # Panel with summary ----
    nav_panel("Summary",
              DT::dataTableOutput("stats")
    )
  )
)


server <- function(input, output, session) {
  
  output$dynamic_prior_inputs <- renderUI({
    num <- as.integer(input$num_prior_params)
    
    prior_inputs <- lapply(1:num, function(i) {
      fluidRow(
        column(4, numericInput(paste0("prior_mu_", i), paste("prior mean", i), 0, step = 1)),
        column(4, numericInput(paste0("prior_sd_", i), paste("prior sd", i), 0.14142,min = 0.001, step = 1)),
        column(4, numericInput(paste0("weight_", i), paste("weight", i), 1/i,min = 0, max = 1, step = 0.05))
      )
    })
    
    do.call(tagList, prior_inputs)
  })
  
  observeEvent(input$compute, {
    prior_mus <- sapply(1:input$num_prior_params, function(i) input[[paste0("prior_mu_", i)]])
    prior_sds <- sapply(1:input$num_prior_params, function(i) input[[paste0("prior_sd_", i)]])
    w <- sapply(1:input$num_prior_params, function(i) input[[paste0("weight_", i)]])
    
    if (length(prior_mus) != length(prior_sds) || length(prior_mus) != length(w) || length(prior_sds) != length(w)) {
      showNotification("Error: Check that the length of prior means and standard deviations parameters must be EQUAL",
                       type = "error", duration = 5)
      return(NULL)
    }
    
    
    if (sum(w) != 1) {
      showNotification("Error: Check your inputs, the weights MUST add up to 1.",
                       type = "error", duration = 5)
      return(NULL)
    }
    
    
    mu <- input$mu
    tau <- input$tau
    n <- input$n
    priorPars <- matrix(c(prior_mus, prior_sds), nrow = 2, byrow = TRUE)
    alpha_level <- input$alpha_level
    
    if (tau==0) {
      showNotification("Error: Check your inputs, the standard deviation cannot be 0.",
                       type = "error", duration = 5)
      return(NULL)
    }
    
    results <- rhpdi(priorPars, w, mu, tau, n, alpha_level)
    
    if (!is.null(results)) {
      output$plot_hpd <- renderPlot({
        rhpdi(priorPars, w, mu, tau, n, alpha_level)
      })
      
      output$hpd_text <- renderText({
        results[[1]]
      })
      
      output$stats <- DT::renderDataTable({
        tb <- data.frame(
          "Posterior Weights" = round(results[[2]], 4),
          "Posterior Means" = round(results[[3]], 4),
          "Posterior SDs" = round(results[[4]], 4)
        )
        rownames(tb) <- paste0("Posterior Component ", 1:nrow(tb))
        colnames(tb) <- c("Posterior Weights", "Posterior Means", "Posterior SDs")
        tb
      }, rownames = TRUE, options = list(pageLength = 3, scrollX = TRUE, autoWidth = TRUE))
    }
  })
}

shinyApp(ui = ui, server = server)