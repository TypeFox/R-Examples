# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Change in Genotype Frequency"),
  
  sidebarPanel(
    sliderInput("initial.A", "Initial Frequency of A allele:", 
                min = 0, max = 1, value=.5, step =.025),
    sliderInput("pop", "Population Size:", 
                min = 10, max = 1000, value = 100, step = 10),
    sliderInput("gen", "Generations to simualte:", 
                min = 10, max = 500, value = 100, step = 10),
    sliderInput("fit.AA", "Fitness of AA:", 
                min = 0, max = 1, value = 1, step = .05),
    sliderInput("fit.Aa", "Fitness of Aa:", 
                min = 0, max = 1, value = 1, step = .05),
    sliderInput("fit.aa", "Fitness of aa:", 
                min = 0, max = 1, value = 1, step = .05),
    sliderInput("iter", "Iterations:", 
                min = 1, max = 50, value = 10, step = 1),
    selectInput("var.plot", "Plot:",
                list("A" = 4, 
                     "a" = 5, 
                     "AA" = 1,
                     "Aa"  = 2,
                     "aa"  = 3)),
    checkboxInput(inputId = "traj",
                  label = "Show expected outcome",
                  value = FALSE),
    actionButton("seed.val", 'Refresh'),
    #selectInput("heath", "Benchmarking Tests:",
    #                list("option 1" = "preset", 
    #                     "option 2" = "fly")),\
    sliderInput("qAa", "Mutation Rate from A to a:", 
                min = 0, max = .2, value = 0, step = .02),
    sliderInput("qaA", "Mutation Rate from a to A:", 
                min = 0, max = .2, value = 0, step = .02),
    sliderInput("width", "Line width:", 
                min = .2, max = 6, value = 3, step = .1)
  ),  
  mainPanel(
       h4(textOutput("caption1")),
       h4(textOutput("caption2")),
       plotOutput("genePlot", width='100%', height='600px')
  )
))