library(shiny)

shinyUI(fluidPage(
  titlePanel(HTML("OncoType DX Recurrence Scores and Outcomes<br>Number Needed to Treat (NNT), with Adverse Events (AE) Details")),

  fluidRow(
    column(1, ""),
    column(11,
           div(class="well container-fluid",
               #               fluidRow(style="horizontal-align:center; text-align:center;",
               column(width=6,
                      numericInput(inputId = "RSchosen", label = "Change the Recurrence Score (RS) here:", value = 30, min = 5, max = 50)
               ),
               #               style="horizontal-align:center; text-align:center;"),
               column(width=6,
                      p(style = "font-weight:bold", "Or set the RS to a standard cutoff here."),
                      actionButton("ODXlow", "ODXlow RS = 17", style = "background-color:#FFE8E8"),
                      actionButton("ODXhigh", "ODXhigh RS = 30", style = "background-color:#FFE8E8"),
                      actionButton("TailorX_low", "TailorX_low = 11", style = "background-color:#E8E8FF"),
                      actionButton("TailorX_high", "TailorX_high = 25", style = "background-color:#E8E8FF")
               )
           ),
           fluidRow(
             column(7,
                    plotOutput(outputId = "skislope",
                               click="skislope_click",
                               hover="skislope_hover"),
                    numericInput(inputId = "ytop", label = "Upper limit of NNT scale", value = 100, min = 10, max = 200, step = 1)),

             column(5,
                    plotOutput(outputId = "AEstack")
             )
           )
    ))
)
)

