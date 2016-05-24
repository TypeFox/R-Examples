library(shiny)

shinyUI(pageWithSidebar(
  headerPanel("dPCR simulation"),
  sidebarPanel(
    numericInput("k", "Positive partitions:", 45, min = 0, step = 1),
    numericInput("n", "Total number of partitions:", 65, min = 0, step = 1),
    numericInput("vol", "Volume of a single partition (nL):", 5, min = 0, 
                 step = 1),
    selectInput("methods", "Choose method of calculating CI:", choices  = 
                  list("exact", "agresti-coull", "asymptotic", "wilson", "prop.test", 
                         "bayes", "logit", "cloglog", "probit", "lrt"), "wilson"),
    numericInput("conf.level", "Confidence interval level:", 0.95, min = 0, max = 1,
                 step = 0.01)
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("Summary (theoretical distribution)", verbatimTextOutput("summary")),
      tabPanel("Theoretical distribution", plotOutput("distr"), 
               tableOutput("table"), plotOutput("pois1"), verbatimTextOutput("volume"))
    )
  )
))

