library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("Bayesian priors: beta and binomial"),

  # Sidebar with a slider input for number of bins
fluidRow(
  p("In Bayesian Statistics, we use prior information (gathered
     from previous research, intuition, or well-informed hypotheses)
    as part of our statistical formulation. This document focuses on 
    working with a prior for the parameter p in a Binomial(n,p) distribution. 
    The beta distribution is the conjugate prior for the Binomial, 
    because we start with a beta prior distribution, utilize the likelihood from our experiment,
    then end up with a beta posterior distribution that uses both the prior 
    and the experimental likelihood.")),
fluidRow(
  column(4,plotOutput("betapriorPlot")),
  column(4,plotOutput("binlikPlot")),
  column(4,plotOutput("betapostPlot"))),
fluidRow(
      column(4,sliderInput("alpha",
                  "Alpha Parameter:",
                  min = 0.1,
                  max = 5,
                  value = 1,
                  step=0.1),
      sliderInput("beta",
                  "Beta Parameter:",
                  min = 0.1,
                  max = 5,
                  value = 1,
                  step=0.1)),
      column(4,sliderInput("n",
                  "Number of trials for the binomial:",
                  min = 1,
                  max = 100,
                  value = 1,
                  step=1),
      htmlOutput("sUI")),
      column(4)
    ),
fluidRow(
  p("As you can see, you can choose a Beta distribution that looks like what you expect
    the distribution of p to look like. We start with a default beta of (1,1), which is
    equivalent to a Uniform(0,1) prior on p, but you set the prior distribution to whatever you'd like."),
  p("Note how the Posterior likelihood follows a Beta (alpha + s,beta + n - s) distribution. "))
  )
)
