#####  shinyROC ui
appName = "shinyROC"
cat("Launching ", appName, "\n")

ui = shinyUI(fluidPage(
  tags$head(tags$script(HTML('
      Shiny.addCustomMessageHandler("jsCode",
        function(message) {
          console.log(message)
          eval(message.code);
        }
      );
    '
  ))),
  uiOutput("debugTools"),
  fluidRow(
    column(7,
           h1("ROC-type curves")),
    column(5, a(
      href="Using_the_NNTbiomarker_package.htm", rel="help", target="_blank",
      #href="../doc/Using_the_NNTbiomarker_package.html", rel="help", target="_blank",
      #href="http://www.github.org/professorbeautiful/NNTbiomarkerHome/man/Using_the_NNTbiomarker_package.html",
      fluidRow(
        column(3,
               style="background:yellow",
               strong(em("Click for information:",
                         style="color:darkgreen")))
        ,
        column(1, style="background:yellow",
               actionButton(inputId = "Info", label="",
                            style="background:lightgreen",
                            icon=icon("info-sign", lib="glyphicon")))
      )
    )
    )
  ),
  hr(),
  fileInput(inputId = "data", label = "Choose data file"),
  numericInput(inputId = "N", "N", value=1000, min = 4, step = 1),
  plotOutput(outputId = "PVplotID",
             height = 350, width=350,
             # Equivalent to: click = clickOpts(id = "plot_click")
             click = "plot_click",
             dblclick = dblclickOpts(id = "plot_dblclick"),
             hover = hoverOpts(id = "plot_hover"),
             brush = brushOpts(id = "plot_brush")
  )
)
)

server = shinyServer(function(input, output, session) {
  rValues = reactiveValues()
  output$PVplotID = renderPlot({
    rValues$plotReturnValue = ROCplots(whichPlots = "pv")
  })
  output$hover_info <- renderUI({
    cat("input$plot_hover:\n")
    hoverData = input$plot_hover
    ppv = hoverData$x
    npv = hoverData$y
    argmin()

  })
})

shinyApp(ui = ui, server = server)
