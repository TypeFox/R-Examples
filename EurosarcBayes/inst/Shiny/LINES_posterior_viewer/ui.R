library(shiny)



# ui
shinyUI(fluidPage(
  theme="bootstrap_cerulean.css",
  sidebarLayout(
    sidebarPanel(width=3,
      tags$head(
        tags$style(type="text/css", "select { max-width: 360px; }"),
        tags$style(type="text/css", ".span4 { max-width: 360px; }"),
        tags$style(type="text/css",  ".well { max-width: 360px; }")
      ),

      h4("Values of interest"),

      h5("Response"),

      div(id="rlower1",numericInput(inputId="rlower", label="Unacceptable", value = 0.2,min=0,max=1)),
      tags$head(tags$style(type="text/css", "#rlower1 {display: inline-block}")),
      tags$head(tags$style(type="text/css", "#rlower {max-width: 100px}")),

      div(id="rupper1",numericInput(inputId="rupper", label="Desired", value = 0.35,min=0,max=1)),
      tags$head(tags$style(type="text/css", "#rupper1 {display: inline-block}")),
      tags$head(tags$style(type="text/css", "#rupper {max-width: 100px}")),

      h5("Toxicity"),

      div(id="tlower1",numericInput(inputId="tlower", label="Desired", value = 0.1,min=0,max=1)),
      tags$head(tags$style(type="text/css", "#tlower1 {display: inline-block}")),
      tags$head(tags$style(type="text/css", "#tlower {max-width: 100px}")),

      div(id="tupper1",numericInput(inputId="tupper", label="Unacceptable", value = 0.3,min=0,max=1)),
      tags$head(tags$style(type="text/css", "#tupper1 {display: inline-block}")),
      tags$head(tags$style(type="text/css", "#tupper {max-width: 100px}")),

      selectInput(inputId="resp.endpoint", label="Response endpoint", choices=c("None","Futility","Efficacy")),
      selectInput(inputId="tox.endpoint", label="Toxicity endpoint", choices=c("None","Toxicity","No Toxicity")),

      h4("Prior information"),

      h5("Response"),
      div(id="rprior1",numericInput(inputId="rpriora", label="Success", value = 1)),
      tags$head(tags$style(type="text/css", "#rprior1 {display: inline-block}")),
      tags$head(tags$style(type="text/css", "#rpriora {max-width: 100px}")),

      div(id="rprior2",numericInput(inputId="rpriorb", label="Failure", value = 1)),
      tags$head(tags$style(type="text/css", "#rprior2 {display: inline-block}")),
      tags$head(tags$style(type="text/css", "#rpriorb {max-width: 100px}")),

      h5("Toxicity"),
      div(id="tprior1",numericInput(inputId="tpriora", label="Toxicity", value = 1)),
      tags$head(tags$style(type="text/css", "#tprior1 {display: inline-block}")),
      tags$head(tags$style(type="text/css", "#tpriora {max-width: 100px}")),

      div(id="tprior2",numericInput(inputId="tpriorb", label="No toxicity", value = 1)),
      tags$head(tags$style(type="text/css", "#tprior2 {display: inline-block}")),
      tags$head(tags$style(type="text/css", "#tpriorb {max-width: 100px}")),

      h4("Data collected"),

      h5("Response"),
      div(id="rdata1",numericInput(inputId="rdataa", label="Success", value = 0)),
      tags$head(tags$style(type="text/css", "#rdata1 {display: inline-block}")),
      tags$head(tags$style(type="text/css", "#rdataa {max-width: 100px}")),

      div(id="rdata2",numericInput(inputId="rdatab", label="Failure", value = 0)),
      tags$head(tags$style(type="text/css", "#rdata2 {display: inline-block}")),
      tags$head(tags$style(type="text/css", "#rdatab {max-width: 100px}")),

      h5("Toxicity"),
      div(id="tdata1",numericInput(inputId="tdataa", label="Toxicity", value = 0)),
      tags$head(tags$style(type="text/css", "#tdata1 {display: inline-block}")),
      tags$head(tags$style(type="text/css", "#tdataa {max-width: 100px}")),

      div(id="tdata2",numericInput(inputId="tdatab", label="No toxicity", value = 0)),
      tags$head(tags$style(type="text/css", "#tdata2 {display: inline-block}")),
      tags$head(tags$style(type="text/css", "#tdatab {max-width: 100px}"))



    ),

    mainPanel(width=9,
      plotOutput("prior.graph",  width = "auto"),
      plotOutput("posterior.graph",  width = "auto")
    )
  ))
)
