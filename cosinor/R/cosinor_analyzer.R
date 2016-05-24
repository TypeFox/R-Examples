#' @import shiny
NULL

#' Shiny application to demonstrate cosinor fit
#'
#' Given a dataset, specify the outcome, time variable, and optional covariates. The app will then perform a cosinor analysis and plot the results.
#'
#' @param data Data frame to analyze
#'
#'
#' @examples
#' \dontrun{
#' library(shiny)
#' cosinor_analyzer(vitamind)
#' }
#'
#' @export
#'
#'
cosinor_analyzer <- function(data = vitamind){

  shinyApp(

    ui = fluidPage( title = paste("Cosinor analysis of", as.character(bquote(data))),
      sidebarLayout(
        sidebarPanel(
          selectInput("Y", "Outcome variable", choices = colnames(data)),
          selectInput("time", "Time variable", choices = colnames(data)),
          selectInput("X", "Covariates affecting the mean", choices = colnames(data), multiple = TRUE),
          selectInput("X2", "Covariates affecting the amplitude and acrophase", choices = colnames(data), multiple = TRUE),
          textInput("period", "Period length", value = "12"),
          checkboxInput("points", "Show points", value = TRUE),
          actionButton("go", label = "Apply Changes", icon = icon("play")),
          helpText("Testing options"),
          selectInput("tX", "Variable(s) to test", choices = colnames(data), multiple = TRUE),
          radioButtons("param", "Parameter to test", choices = list(Amplitude = "amp", Acrophase = "acr"), inline = TRUE),
          actionButton("test", label = "Run Test", icon = icon("play"))
          ),



        mainPanel(
          plotOutput("mainPlot"),
          verbatimTextOutput("summary"),
          verbatimTextOutput("test")


          )
        )
      ),

    server = function(input, output){

      myFit <- reactive({


        ## make formula
        tv <- paste0("time(", input$time, ")")

        if(is.null(input$X)) { covs <- NULL } else {
          covs <- paste(input$X, collapse = " + ")
        }
        if(is.null(input$X2)) { amps <- NULL } else {
          amps <- paste(paste0("amp.acro(", input$X2, ")"), collapse = " + ")
        }

        myform <- as.formula(paste(input$Y, "~", paste(c(tv, covs, amps), collapse = " + ")))

        cosinor.lm(myform, data, period = as.numeric(input$period))



      })

      output$mainPlot <- renderPlot({

        if(input$go == 0) return()



         p1 <-  ggplot.cosinor.lm(myFit(), input$X)

         if(input$points){

           p1 + geom_point(data = data, aes_string(y = input$Y, x = input$time, col = NULL))

         } else p1

        })


      output$summary <- renderPrint({

        if(input$go == 0) return()
        summary(myFit())

        })

      output$test <- renderPrint({

        if(input$test == 0) return()
        test_cosinor(myFit(), input$tX, input$param)

      })

    }

    )

}


