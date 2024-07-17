if (!requireNamespace("shiny", quietly = TRUE)) {
  install.packages("shiny")
}
if (!requireNamespace("bslib", quietly = TRUE)) {
  install.packages("bslib")
}
if (!requireNamespace("rootSolve", quietly = TRUE)) {
  install.packages("rootSolve")
}

library(shiny)
library(bslib)
library(rootSolve)

# Define the hpdi function (same as your existing function)
options(digits = 5)
hpdi <- function(priorPars, w, mu, tau, n, alpha_level = 0.85) {
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
       xlim = c(x.min, x.max), ylim = c(0, up.lim))
  
  #add the likelihood
  lines(x, lik, type = "l", lwd = 1, lty = 2, col = "darkblue")
  
  lines(x, posterior_density, type = "l", lwd = 2, lty = 1, col = "black")
  
  legend_labels <- c("Mixture prior density","Likelihood", 
                     "Posterior mixture density","HPD Interval")
  legend("topleft", lty = c(2,2,1,2),
         legend = legend_labels, lwd = c(1,1,2,1),
         col = c("darkorange","darkblue", "black","red"), bty = "n")
  
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


ui <- page_sidebar(
  title = "HPD Interval Computation",
  theme = bslib::bs_theme(bootswatch = "united"),
  
  sidebar = sidebar(
    width = 400,
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
    numericInput("mu", HTML("mean (&#956;)"), 0.78),
    numericInput("tau", HTML("standard deviation (&#964;)"), 1),
    numericInput("n", "sample size (n):", 20),
    h4("Prior parameters"),
    p("(comma-separated)"),
    textInput("prior_mus",  HTML("Prior means (&#956;)"), "0, 0, 0"),
    textInput("prior_sds",  HTML("Prior standard deviations (&#964;)"), "0.1414214, 10, 1"),
    textInput("weights", "weights:", "0.50, 0.46,0.04"),
    numericInput("alpha_level", "alpha level:", 0.90, min = 0, max = 1, step = 0.05),
    actionButton("compute", "Compute", class = "btn-lg btn-success"),
    class = "sidebar"
  ),
  
  navset_card_underline(
    # Panel with intro ----
    nav_panel("About",
              p("The WebApp rHPDI is a web application that computes the highest posterior density interval (HPDI) for a given set of informative and vague prior parameters, weights, likelihood parameters and desired alpha level."),
              h4("Parameters:"),
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
    nav_panel("Plot",         
              plotOutput("plot_hpd", width = "95%", height = "900px"),
              verbatimTextOutput("hpd_text")
    ),
    
    # Panel with summary ----
    nav_panel("Summary", 
              verbatimTextOutput("summary")
    )
  )
)


server <- function(input, output, session) {
  
  parse_input <- function(input) {
    as.numeric(unlist(strsplit(input, ",")))
  }
  
  mu <- reactive(input$mu)
  tau <- reactive(input$tau)
  n <- reactive(input$n)
  
  observeEvent(input$compute, {
    prior_mus <- parse_input(input$prior_mus)
    prior_sds <- parse_input(input$prior_sds)
    weights   <- parse_input(input$weights)
    weights <- weights / sum(weights)
    
    
    if (length(prior_mus) != length(prior_sds)) {
      showNotification("Error: The number of prior means and standard deviations must be the same.",
                       type = "error", col = "red")
      return()
    }
    
    priorPars <- matrix(c(prior_mus, prior_sds), nrow = 2, byrow = TRUE)
    alpha_level <- input$alpha_level
    
    output$plot_hpd <- renderPlot({
      tryCatch({
        hpdi(priorPars, weights, mu(), tau(), n(), alpha_level)
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, "Error in computation: Check your inputs", col = "red")
      })
    })
    
    output$hpd_text <- renderText({
      tryCatch({
        hpdi(priorPars, weights, mu(), tau(), n(), alpha_level)[[1]]
      }, error = function(e) {
        showNotification("Error in computation: Check your inputs", type = "error", col = "red")
        NULL
      })
    })
  })
}

shinyApp(ui = ui, server = server)