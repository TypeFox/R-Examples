
#' @title Interactive Model Performance Evaluation & Comparison
#' @description Interactive version of the staticPerfMeasures function
#' @param list_models A list of one (or more) dataframes for each model whose performance is to be evaluated. Each dataframe should comprise of 2 columns with the first column indicating the class labels (0 or 1)
#' and the second column providing the raw predicted probabilities
#' @param sample_size_concord For computing concordance-discordance measures (and c-statistic) a random sample
#' is drawn from each dataset (if nrow(dataset) > 5000). Default sample size of 5000 can be adjusted by changing the value of this
#' argument
#' @param model_function Models can be created interactively, if required. For this option to work, a model function
#' should be passed as an argument. The model function should take a formula as an argument, and return a
#' a dataframe as output (dataframe should comprise of 2 columns with the first column indicating the class labels (0 or 1)
#' and the second column provding the raw predicted probabilities) Refer to the example section for more details
#' @param data The name of the data-set. The Independent Variable (IV) names, for interactive model building, is picked up from this data set
#' @param y The column name of the Dependent Variable (DV), for interactive model building
#'
#' @return
#' This function will launch a ShinyApp.Input parameters (such as the number of bins,
#' the "g" argument in the static version of this function) can be adjusted
#' through app widgets. The 'Run-Analysis' button in the app, will generate model performance
#' output basis selected input parameters
#'
#' For interactive Model building, a model function, data set & the dependent variable name
#' should be passed as arguments. Interactive model building option
#' creates additional input widgets in the app. This includes -
#'
#' A drop down to select independent variables (the names of the variables will be picked up from the data argument)
#'
#' An input slider to include additional models (upto 4 additional models can be created). Each additional model
#' updates the original model created. For e.g. consider the dataset has 10 IVs: x1-x10. Original model
#' was created by selecting x1-x4 from the drop down list. If we need to create a second model, by including x5 and excluding x3 simply type,
#' "+ x5 - x3" in the input text box
#'
#' @export
#' @import dplyr ggplot2 shiny
#' @importFrom tidyr gather
#' @importFrom stats as.formula pchisq quantile
#'
#' @examples
#' # Without interactive model development
#' model_1 <- glm(Species ~ Sepal.Length,data=iris,family=binomial)
#' model_2 <- glm(Species ~ Sepal.Width, data=iris, family = binomial)
#' df1 <- data.frame(model_1$y,fitted(model_1))
#' df2 <- data.frame(model_2$y,fitted(model_2))
#'
#' \dontrun{
#' #This will launch a Shiny App
#' interPerfMeasures(list_models = list(df1,df2))}
#'
#'# With interactive model development
#' glm_model <- function(formula) {
#'    glm_model <- glm(formula, data = iris, family = "binomial")
#'    out <- data.frame(glm_model$y, fitted(glm_model))
#'    out }
#'  \dontrun{
#'  #This will launch a Shiny App
#'  interPerfMeasures (model_function = glm_model,data=iris,y="Species")}
#'

interPerfMeasures <- function(list_models,sample_size_concord = 5000,model_function = NULL,data = NULL,y = NULL) {

  if (is.null(model_function)) {
    error_handler(list_models,arg=0,method = "Performance_Measure")
  }

  return_fun_value <- function(model_function,form) {

    model_function(form)
  }



  shinyApp(
    ui = fluidPage(
      sidebarLayout(
        sidebarPanel(uiOutput("build_model_interactively"),
                     sliderInput("bins", "Select no of Groups", 1, 30, 10),
                     selectInput("perf_measures",
                          label = strong("Select Performance Measure"),
                          choices = list("Hosemr Lemeshow" = "hosmer", "Calibration Plot" = "calibration",
                                          "Lift Index"="lift","Conordance-Discordance"="concord")),

                     conditionalPanel(
                         condition = " input.build_model == 'Y' ",
                         uiOutput("variables"),
                         sliderInput("numModels","No. of additional models",
                          min=0,max=4,value=0),
                         uiOutput("addModels")),

                     actionButton("perf_measure_button","Run Analysis")),

        mainPanel(conditionalPanel(
                    condition = " input.perf_measures == 'hosmer' ",
                    dataTableOutput("hosmer_results")),

                  h4(strong("Goodness of fit tests"),style="color:brown"),
                  plotOutput("plot_perf_measure"),
                  dataTableOutput("df_perf_measure")
        )
      )
    ),

    server = function(input, output, session) {

      perf_func_reactive <- reactive({

        if (input$build_model == 'Y')   {
          var_list_combined <- list()
          var_list_combined[[1]] <- paste(input$Var1,collapse="+")

          ModelsCount <- as.numeric(input$numModels)

          if (ModelsCount > 0) {

            for (i in 1:ModelsCount){
                if ((input[[paste("Mod",i+1)]]) != "") {
                  var_list_i <- paste(var_list_combined[[1]],input[[paste("Mod",i+1)]],sep="")
                  var_list_combined[[i+1]] <- var_list_i
                }
            }
          }

          run_formula <- function(x) {
            formula_val <- as.formula(paste(y," ~ ",x,sep=""))
            df=return_fun_value(model_function,formula_val)
            df}


          out_df <- lapply(seq_along(var_list_combined),function(i) run_formula(var_list_combined[[i]]))
          out <- staticPerfMeasures(out_df,g = input$bins,perf_measures = input$perf_measures, sample_size_concord)

        } else {
          g <- as.numeric(input$bins)
          out <- staticPerfMeasures(list_models,g,perf_measures = input$perf_measures,sample_size_concord)
        }


        if (input$perf_measures == 'hosmer'){
          perf_df <- out$data$hosmer_df
          hosmer_test <- out$data$hosmer_results
          perf_plot <- out$plots$hosmer
          return(list(df=perf_df,plot=perf_plot,hosmer_test=hosmer_test))
        }

        else{
          perf_df <- out[[1]]
          perf_plot <- out[[2]]
          return(list(df=perf_df,plot=perf_plot))
        }


    })

      button_click <- eventReactive(input$perf_measure_button,{

        perf_func_reactive()

      })

      output$plot_perf_measure <- renderPlot({

        button_click()$plot

      })

      output$df_perf_measure <- renderDataTable({

        df <- as.data.frame(button_click()$df)
        df

      })

      output$hosmer_results <- renderDataTable({

        df <- as.data.frame(button_click()$hosmer_test)
        df

      })

      output$variables <- renderUI({

        selectizeInput("Var1","Select Independent Variables",colnames(data),multiple = T)

      })

      output$addModels <- renderUI({

        ModelsCount <- as.integer(input$numModels)
        if (ModelsCount > 0){
          lapply(1:ModelsCount,function(i) {

            textInput(inputId = paste("Mod",i+1),label="Update Original Model",value="")
          })
        }

      })

      output$build_model_interactively <- renderUI({

        if (!is.null(model_function)) {
          radioButtons("build_model","Build Model Interactively?",
                     choices=c("Yes"="Y","No"="N"),selected="Y")
        } else {

          radioButtons("build_model","Build Model Interactively?",
                       choices=c("Yes"="Y","No"="N"),selected="N")
        }

      })

      session$onSessionEnded(function() {

        stopApp()

      })
    }
  )
}


###################### Interactive Confusion ##################################


#' @title Interactive confusion matrix
#' @description Interactive version of the staticConfMat function
#' @param list_models A list of one (or more) dataframes for each model whose performance is to be evaluated. Each dataframe should comprise of 2 columns with the first column indicating the class labels (0 or 1)
#' and the second column provding the raw predicted probabilities
#' @param model_function Models can be created interactively, if required. For this option to work, a model function
#' should be passed as an argument. The model function should take a formula as an argument, and return a
#' a dataframe as output (dataframe should comprise of 2 columns with the first column indicating the class labels (0 or 1)
#' and the second column provding the raw predicted probabilities) Please refer to the example section for more details
#' @param data The name of the data-set. The Independent Variable (IV) names, for interactive model building, is picked up from this data set
#' @param y The column name of the Dependent Variable (DV), for interactive model building
#'
#' @return
#' This function will launch a ShinyApp.Input parameters (such as the probability threshold,
#' the "t" argument in the static version of this function) can be adjusted
#' through app widgets. The 'Run-Analysis' button in the app, will generate model performance
#' output basis selected input parameters
#'
#' For interactive Model building, a model function, data set & the dependent variable name
#' should be passed as arguments. Interactive model building option
#' creates additional input widgets in the app. This includes -
#'
#' A drop down to select independent variables (the names of the variables will be picked up from the data argument)
#'
#' An input slider to include additional models (upto 4 additional models can be created). Each additional model
#' updates the original model created. For e.g. consider the dataset has 10 IVs: x1-x10. Original model
#' was created by selecting x1-x4 from the drop down list. If we need to create a second model, by including x5 and excluding x3 simply type,
#' "+ x5 - x3" in the input text box
#'
#'
#' @export
#' @import dplyr ggplot2 shiny
#' @importFrom tidyr gather
#' @importFrom stats as.formula pchisq quantile
#' @examples
#'# Without interactive model development
#' model_1 <- glm(Species ~ Sepal.Length,data=iris,family=binomial)
#' model_2 <- glm(Species ~ Sepal.Width, data=iris, family = binomial)
#' df1 <- data.frame(model_1$y,fitted(model_1))
#' df2 <- data.frame(model_2$y,fitted(model_2))
#' \dontrun{
#' #This will launch a Shiny App
#' interConfMatrix(list_models = list(df1,df2))}
#'
#' # With interactive model development
#' glm_model <- function(formula) {
#'    glm_model <- glm(formula, data = iris, family = "binomial")
#'    out <- data.frame(glm_model$y, fitted(glm_model))
#'    out }
#'  \dontrun{
#'  #This will launch a Shiny App
#'  interConfMatrix(model_function=glm_model,data=iris,y="Species")}
#'


interConfMatrix <- function(list_models,model_function=NULL,data=NULL,y=NULL) {

  if (is.null(model_function)) {
    error_handler(list_models,arg=0,method = "Confusion")
  }

  return_fun_value <- function(model_function,form) {

    model_function(form)
  }

  shinyApp(
    ui = fluidPage(
      sidebarLayout(
        sidebarPanel(uiOutput("build_model_interactively"),

                     conditionalPanel(
                       condition = " input.build_model == 'Y' ",
                       uiOutput("variables"),
                       sliderInput("numModels","No. of additional models",
                                   min=0,max=4,value=0),
                       uiOutput("addModels")),

                     radioButtons("Option",'Select Option',
                        choices=list("Confusion Matrix"="matrix","Plots & Metrics"="plot"),selected = "matrix"),

                     conditionalPanel(
                       condition = " input.Option == 'matrix' ",
                       sliderInput("Threshold", "Set Threshold", 0, 1, 0.1)),

                     conditionalPanel(
                       condition = " input.Option == 'plot' ",
                       selectInput("Metric","Select Metric",
                         choices = list("Accuracy"="Acc","TPR"="TPR","FPR"="FPR","Precision"="Prec")),
                       sliderInput("Reps", "No of repetitions", 10, 100, 10)),

                     actionButton("Confusion_Button","Run Analysis")),
        mainPanel(
          conditionalPanel(
            condition = " input.Option == 'matrix' ",
            h4(strong("Confusion Matrix"),style="color:brown"),
            dataTableOutput('confusion_df'),
            dataTableOutput('metrics_df')),

        conditionalPanel(
            condition = "input.Option=='plot'",
            h4(strong("Selected Metric vs Prob Threshold"),style="color:brown"),
            plotOutput('metrics_range_plot'),
            dataTableOutput('metrics_range_df'))

      )
    )
  ),
    server = function(input, output,session) {

      confusion_func_shiny <- reactive({

        if (input$build_model == 'Y')   {
          var_list_combined <- list()
          var_list_combined[[1]] <- paste(input$Var1,collapse="+")
          ModelsCount <- as.integer(input$numModels)

          if (ModelsCount > 0) {

            for (i in 1:ModelsCount){
              if ((input[[paste("Mod",i+1)]]) != "") {
                var_list_i <- paste(var_list_combined[[1]],input[[paste("Mod",i+1)]],sep="")
                var_list_combined[[i+1]] <- var_list_i
              }
            }
          }


          build_model <- function(x) {
            formula_val <- as.formula(paste(y," ~ ",x,sep=""))
            df=return_fun_value(model_function,formula_val)
            df}


          out_df <- lapply(seq_along(var_list_combined),function(i) build_model(var_list_combined[[i]]))
          out <- staticConfMatrix(out_df,t= as.numeric(input$Threshold),reps=input$Reps)
        } else {

          t <- as.numeric(input$Threshold)
          out <- staticConfMatrix(list_models,t=t,reps=input$Reps)
        }

        out_list <- list(metrics_range=out[[1]],confusion=out[[2]],metrics=out[[3]])

      })

      confusion_event <- eventReactive(input$Confusion_Button,{

        confusion_func_shiny()

      })

      output$confusion_df <- renderDataTable({

        as.data.frame(confusion_event()$confusion)

      })
      output$metrics_df <- renderDataTable({

        as.data.frame(confusion_event()$metrics)

      })

      output$metrics_range_plot <- renderPlot({

        out_metrics_range <- as.data.frame(confusion_event()$metrics_range)

        Metric <- as.character(input$Metric)

        plot_val <- ggplot(out_metrics_range,aes_string(x="Threshold",y=Metric,color="Model")) + geom_line()
        plot_val

      })

      output$metrics_range_df <- renderDataTable({

        as.data.frame(confusion_event()$metrics_range)

      })
      output$variables <- renderUI({

        selectizeInput("Var1","Select Independent Variables",colnames(data),multiple = T)
      })

      output$addModels <- renderUI({

        ModelsCount <- as.integer(input$numModels)
        if (ModelsCount > 0){
          lapply(1:ModelsCount,function(i) {

            textInput(inputId = paste("Mod",i+1),label="Update Original Model",value="")
          })
        }

      })

      output$build_model_interactively <- renderUI({

        if (!is.null(model_function)) {
          radioButtons("build_model","Build Model Interactively?",
                       choices=c("Yes"="Y","No"="N"),selected="Y")
        } else {

          radioButtons("build_model","Build Model Interactively?",
                       choices=c("Yes"="Y","No"="N"),selected="N")
        }

      })


      session$onSessionEnded(function() {

        stopApp()

      })

    }
  )
}
