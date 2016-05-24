# Define UI for bd trees
shinyUI(pageWithSidebar(
  headerPanel("Probability Density of Distributions"),    # Application title
  sidebarPanel(
    selectInput("select", label = "Select Statistical Distribution",
                choices = list("Normal" = 1, "Exponential" = 2,
                               "Gamma" = 3, "Logistic" = 4, "Poisson" = 5,
                               "Beta" = 6), selected = 1), 
    # normal
    conditionalPanel(
      condition = "input.select == 1",
      sliderInput("mu", "mu:", 
                  min = -100, max = 100, value=0, step = 5),
      sliderInput("sigma", "sigma:", 
                  min = 0.1, max = 10, value=1, step = .1),
      sliderInput("n", "sample size:", 
                  min = 2, max = 1000, value=100, step = 1)
    ),
    # exponential    
    conditionalPanel(
      condition = "input.select == 2",
      sliderInput("lambda", "lambda:", 
                  min = 0.1, max = 10, value=1, step = .1)
    ),
    # gamma
    conditionalPanel(
      condition = "input.select == 3",
      sliderInput("kappa", "kappa:", 
                  min = 0.1, max = 10, value=1, step = .1),
      sliderInput("theta", "theta:", 
                  min = 0.3, max = 10, value=1, step = .1)
    ),
    # logistic
    conditionalPanel(
      condition = "input.select == 4",
      sliderInput("mu2", "mu:", 
                  min = -100, max = 100, value=0, step = 5),
      sliderInput("sigma2", "sigma:", 
                  min = 0.1, max = 10, value=1, step = .1)
    ),
    # poisson
    conditionalPanel(
      condition = "input.select == 5",
      sliderInput("lambda2", "lambda:", 
                  min = 1, max = 32, value=1, step = 1)
    ),   
    # beta
    conditionalPanel(
      condition = "input.select == 6",
      sliderInput("alpha", "alpha:", 
                  min = 0.1, max = 10, value=1, step = .1),
      sliderInput("beta", "beta:", 
                  min = 0.1, max = 10, value=1, step = .1)
    )   
  ),  
  mainPanel(
    plotOutput("treePlot", width='100%', height='600px')
  )
))