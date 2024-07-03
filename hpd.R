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

# Define the hpdi function
hpdi <- function(priorPars, w, mu, tau, n, alpha_level = 0.85, x_range = c(-1, 2)) {
  options(scipen = 999)
  
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
    
    posterior_vals <- 0
    for (i in 1:length(post.w)) {
      posterior_vals <- posterior_vals + post.w[i] * dnorm(x, mean = post.mean[i], sd = sqrt(post.var[i]))
    }
    
    return(list(f = posterior_vals, post.w = post.w, post.m = post.mean, post.s = sqrt(post.var)))
  }
  
  x <- seq(x_range[1], x_range[2], by = 0.0001)
  y <- seq(max(fx(x, priorPars, w, mu, tau, n)$f), 0, by = -0.0001)
  post.w <- fx(x, priorPars, w, mu, tau, n)$post.w
  post.m <- fx(x, priorPars, w, mu, tau, n)$post.m
  post.s <- fx(x, priorPars, w, mu, tau, n)$post.s
  postPars <- matrix(c(post.m, post.s), nrow = 2, byrow = TRUE)
  
  my_alpha <- function(y) {
    
    roots <- rootSolve::uniroot.all(function(x) fx(x, priorPars, w, mu, tau, n)$f - y, lower = min(x), upper = max(x), tol = 1e-10)
    
    num_roots <- length(roots)
    alpha <- 0
    if (num_roots == 2) {
      
      alpha <- sum(sapply(1:length(w), function(j) {
        
        post.w[j] * pnorm(roots[2], postPars[1, j], postPars[2, j]) - post.w[j] * pnorm(roots[1], postPars[1, j], postPars[2, j])
      }))
      
    } else if (num_roots == 4) {
      
      alpha <- sum(sapply(1:length(w), function(j) {
        
        (post.w[j] * pnorm(roots[2], postPars[1, j], postPars[2, j]) - post.w[j] * pnorm(roots[1], postPars[1, j], postPars[2, j])) +
          (post.w[j] * pnorm(roots[4], postPars[1, j], postPars[2, j]) - post.w[j] * pnorm(roots[3], postPars[1, j], postPars[2, j]))
      }))
    } else {
      alpha <- NA  # when the number of roots is neither 2 nor 4
    }
    return(alpha)
  }
  
  objectivefx <- function(y, ...) { 
    alpha <- my_alpha(y)
    abs(alpha - alpha_level)
  }
  
  opt.res <- optim(0, objectivefx, method = "Brent", lower = min(y), upper = max(y), 
                   priorPars = priorPars, w = w, mu = mu, tau = tau, n = n, 
                   postPars = postPars, alpha_level = alpha_level)
  
  hpd <- rootSolve::uniroot.all(function(x) fx(x, priorPars, w, mu, tau, n)$f - opt.res$par, lower = min(x), upper = max(x),tol = 1e-10)
  
  cmp <- 0
  
  for (i in 1:length(w)) {
    cmp <- cmp + w[i] * dnorm(x, priorPars[1, i], priorPars[2, i])
  }
  
  posterior_density <- fx(x, priorPars, w, mu, tau, n)$f
  
  up.lim <- max(c(cmp, posterior_density))
  
  x.ll <- floor(range(density(post.m)$x)[1]) - 2
  
  x.ul <- ceiling(range(density(post.m)$x)[2]) + 2
  
  plot(x, cmp, type = "l", col = "darkorange", lty = 1, lwd = 2,
       xlab = expression(theta), ylab = "Density", xlim = c(x.ll, x.ul), ylim = c(0, up.lim))
  
  colors <- c("blue", "purple", "darkgreen", "darkcyan", "magenta")
  
  for (i in 1:length(post.w)) {
    lines(x, post.w[i]*dnorm(x, postPars[1, i], postPars[2, i]), type = "l",
          col = colors[i], lty = 2, lwd = 2)
  }
  
  lines(x, posterior_density, type = "l", lwd = 2, lty = 1, col = "black")
  
  legend_labels <- c(paste0("Posterior Component ", 1:length(w), ": N(", round(postPars[1, ], 3), ",", round(postPars[2, ], 3), ")"),
                     "Mixture prior density", "Posterior mixture density", "HPD Interval")
  legend("topright", lty = c(rep(2, length(w)), 1, 1,1),
         legend = legend_labels, lwd = 2,
         col = c(colors[1:length(w)], "darkorange", "black", "red"), bty = "n")
  
  hpd_segments <- function(hpd) {
    
    for (i in seq(1, length(hpd), by = 2)) {
      
      segments(hpd[i], 0, hpd[i], fx(hpd[i], priorPars, w, mu, tau, n)$f, col = "red", lwd = 1, lty = 1)
      segments(hpd[i + 1], 0, hpd[i + 1], fx(hpd[i + 1], priorPars, w, mu, tau, n)$f, col = "red", lwd = 1, lty = 1)
      segments(hpd[i], fx(hpd[i], priorPars, w, mu, tau, n)$f, hpd[i + 1], fx(hpd[i + 1], priorPars, w, mu, tau, n)$f, col = "red", lwd = 1, lty = 1)
    }
  }
  
  if (length(hpd) %in% c(2, 4)) {
    hpd_segments(hpd)
  }
  
  hpd_string <- if (length(hpd) == 2) {
    paste0(alpha_level * 100, "% HPDI is:[", round(hpd[1], 4), ",", round(hpd[2], 4), "]")
  } else if (length(hpd) == 4) {
    paste0(alpha_level * 100, "% HPDI is:[", round(hpd[1], 4), ",", round(hpd[2], 4), "] and [", round(hpd[3], 4), ",", round(hpd[4], 4), "]")
  } else {
    paste("No HPDI for", alpha_level * 100, "%")
  }
  return(hpd_string)
}


ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "united"),
  titlePanel("Computation of highest posterior density (HPD) intervals for robust mixture priors"),
  tags$style(type = "text/css", "#hpd_text { 
             color: blue; 
             text-align: center; 
             font-size: 24px; 
             font-weight: bold;
             }"),
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        column(
          width = 12,
          helpText(
            "Instructions: please provide values for the following parameters;",
            HTML("the mean (&#956;) and standard deviation (&#964;) of the data (likelihood),"),
            "the prior parameters (mu, sd) and weights (w) for the mixture priors,and finally the alpha level for which the highest posterior density interval would be computed."
          )
        )
      ),
      h4("Likelihood parameters"),
      fluidRow(
        column(
          width = 12,
          numericInput("mu", HTML("mean (&#956;)"), 0.78)
        )
      ),
      fluidRow(
        column(
          width = 12,
          numericInput("tau", HTML("standard deviation (&#964;)"), 1)
        )
      ),
      fluidRow(
        column(
          width = 12,
          numericInput("n", "sample size:", 20)
        )
      ),
      fluidRow(
        column(
          width = 12,
          numericInput("alpha_level", "alpha level:", 0.80, min = 0, max = 1, step = 0.05)
        )
      ),
      h4("Prior parameters"),
      fluidRow(
        column(
          width = 6,
          numericInput("prior11", HTML("Prior &#956;<sub>1</sub>:"), 0, step = 0.1)
        ),
        column(
          width = 6,
          numericInput("prior12", HTML("Prior &#956;<sub>2</sub>:"), 0, step = 0.1)
        )
      ),
      fluidRow(
        column(
          width = 6,
          numericInput("prior21",  HTML("Prior &#963;<sub>1</sub>:"), 0.1414214)
        ),
        column(
          width = 6,
          numericInput("prior22",  HTML("Prior &#963;<sub>2</sub>:"), 10)
        )
      ),
      fluidRow(
        column(
          width = 12,
          tags$div(
            style = "display: inline-block; width: 100%;",
            tags$label(style = "display: inline-block; width: 30%;"),
            tags$span("prior weight 2 is computed as 1 - w1", style = "display: inline-block; width: 70%;")
          )
        )
      ),

      fluidRow(
        column(
          width = 12,
          sliderInput("w1", "Prior weight 1:", value = 0.5, min = 0, max = 1)
        )
      ),
      fluidRow(
        column(
          width = 12,
          sliderInput("x_range", "Range of X-values:", min = -10, max = 10, value = c(-1, 2),step = .5)
        )
      ),
      fluidRow(
        column(
          width = 12,
          actionButton("compute", "Compute", class = "btn-lg btn-success")
        )
      )
    ),
    mainPanel(
      plotOutput("plot", width = "100%", height = "900px"),
      verbatimTextOutput("hpd_text")
    )
  )
)

server <- function(input, output, session) {
  
  
  w1 <- reactive({
    input$w1
  })
  
  w2 <- reactive({
    1 - w1()
  })
  

  observeEvent(input$compute, {
    
    priorPars <- matrix(c(input$prior11, input$prior12,
                          input$prior21, input$prior22), 
                        nrow = 2, byrow = TRUE)
    mu <- input$mu
    tau <- input$tau
    n <- input$n
    alpha_level <- input$alpha_level
    x_range <- input$x_range
    
    
    hpd_string <- reactive({hpdi(priorPars, c(w1(), w2()), mu, tau, n, alpha_level, x_range)
    })
    
    
    output$plot <- renderPlot({
      hpdi(priorPars, c(w1(), w2()), mu, tau, n, alpha_level, x_range)
    })
    
    output$hpd_text <- renderText({
      hpd_string()
    })
  })
}


shinyApp(ui = ui, server = server)

