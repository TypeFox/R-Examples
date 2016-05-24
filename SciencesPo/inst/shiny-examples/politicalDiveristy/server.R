library(shiny)
library(ggplot2)

datasets <- list(
  "Random distribution" =
    data.frame("Normal" =
                 rnorm(1000, 150, 20),
               "Bimodal" =
                 c(rnorm(500, 200, 50), rnorm(500, 400, 50)),
               "Uniform" =
                 runif(1000, 100, 200),
               check.names = FALSE),
  "faithful" = faithful,
  "rock" = rock,
  "cars" = cars)

shinyServer(function(input, output, session) {
  # show/hide the marginal plots settings
  observe({
    toggle(id = "marginal-settings", anim = TRUE,
           time = 0.3, condition = input$show_marginal)
  })

  output$dataset_select <- renderUI({
    selectInput("dataset", "Choose a dataset:", names(datasets))
  })

  datasetInput <- reactive({
    if (is.null(input$dataset)) {
      return(NULL)
    }
    datasets[[input$dataset]]
  })

  output$x_var_select <- renderUI({
    dataset <- datasetInput()
    selectInput("x_var", "X variable",
                colnames(dataset), colnames(dataset)[1])
  })

  output$y_var_select <- renderUI({
    dataset <- datasetInput()
    selectInput("y_var", "Y variable",
                colnames(dataset), colnames(dataset)[2])
  })

  # the code to reproduce the plot
  code <- reactive({
    code <- sprintf(paste0(
      "p <- ggplot(`%s`, aes_string('%s', '%s')) +\n",
      "  geom_point() + theme_bw(%s)"),
      input$dataset, input$x_var, input$y_var, fontSize()
    )

    if (input$xtrans == "log") {
      code <- paste0(code, " + scale_x_log10()")
    } else if (input$xtrans == "reverse") {
      code <- paste0(code, " + scale_x_reverse()")
    }
    if (input$ytrans == "log") {
      code <- paste0(code, " + scale_y_log10()")
    } else if (input$ytrans == "reverse") {
      code <- paste0(code, " + scale_y_reverse()")
    }

    code <- paste0(code, "\n\n")

    if (input$show_marginal) {
      code <- paste0(code, sprintf(paste0(
        "ggExtra::ggMarginal(\n",
        "  p,\n",
        "  type = '%s',\n",
        "  margins = '%s',\n",
        "  size = %s,\n",
        "  col = '%s',\n",
        "  fill = '%s'\n",
        ")"),
        input$type, input$margins, size(), input$col,
        input$fill))
    } else {
      code <- paste0(code, "p")
    }
  })

  # hide the loading message
  hide("loading-content", TRUE, "fade")
})
