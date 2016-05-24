
#' @title Launch Tutorial App
#' @description Invoking this function will launch the Tutorial App;
#' The app includes a ReadMe introduction which provides a quick overview on how to use the app
#' @return Launches app
#' @import dplyr shinydashboard ggplot2 rmarkdown shiny
#' @importFrom stats runif
#' @examples
#' # Simply invoke the function to launch the app
#' \dontrun{
#' launch_tutorial()}
#' @export


launch_tutorial <- function() {

  tutorial_set <- list(basic_operations = basic_operations,dplyr_tutorial = dplyr_tutorial,
                       loops_tutorial = loops_tutorial,model_tutorial = model_tutorial)
  diamonds <- ggplot2::diamonds
  diamonds <- diamonds[sample(1:nrow(diamonds),size=1000),]

  marks <- data.frame(replicate(5,round(runif(100,0,100))))
  colnames(marks) <- paste("Subject", seq(1,5,1), sep="")

  Task_List <- Use.Case <- NULL


  shinyApp(ui <- dashboardPage(
    dashboardHeader(title = "R Tutorial"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Read Me", tabName = "readme", icon = icon("dashboard")),
        menuItem("Select a topic", tabName = "topics", icon = icon("dashboard"),
                 menuSubItem(icon=NULL, tabName="topics",
                             selectInput("topic_select","Choose a topic", choices = c("Basic operations on dataset" = "basic_operations", "Data manipulation"="dplyr_tutorial",
                                                                                      "Loops and functions"="loops_tutorial", "Basic model development"="model_tutorial"), selected = "basic_operations"))))
    ),
    dashboardBody(
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"
      ),

      tabItems(
        # First tab content
        tabItem(tabName = "topics",

                fluidRow(
                  box(title = "Select a Use Case", status = "primary", solidHeader = T,
                      DT::dataTableOutput("use_cases")),

                  box(
                    title = "Select a Task", status = "primary", solidHeader = T,
                    uiOutput("func_list")
                  )
                ),

                fluidRow(
                  box(title = "Code Output", status = "primary", solidHeader = T,
                      DT::dataTableOutput("code_out")),

                  box(
                    title = "Code and Comments", status = "primary", solidHeader = T,
                    htmlOutput("comment"))

                )
        ),
        tabItem(tabName = "readme",
                fluidRow(
                  box(title = "Read Me", status = "primary", solidHeader = T, width = 12,
                      htmlOutput("read_me"))))
      ))),


    server <- function(input, output, session) {

      diamonds_sample <- diamonds[sample(nrow(diamonds),size=200),]

      output$func_list <- renderUI({

        # Extract the set of tasks/functions for the selected topic
        selected_topic <- tutorial_set[[input$topic_select]]
        task_list = unique(as.character(selected_topic$Task_List))

        # Dynamically update the radiobutton widget
        radioButtons("func",'', choices = task_list)


      })

      output$use_cases <- DT::renderDataTable({

        extract_usecase <- tutorial_set[[input$topic_select]] %>% filter(Task_List == input$func) %>%
                                                                  select(Use.Case)
        DT::datatable(data = extract_usecase, selection = 'single',
                      options = list(scrollX = TRUE, rownames= FALSE))

      })



      output$code_out <- DT::renderDataTable({

        selected_topic <- tutorial_set[[input$topic_select]] %>% filter(Task_List == input$func)

        index <- as.numeric(input$use_cases_row_last_clicked)
        form <- as.character(selected_topic[index,"Formula"])


        code_output <- eval(parse(text = form))


        DT::datatable(data = code_output, selection = 'single',
                      options = list(scrollX = TRUE, rownames= FALSE))

      })

      output$comment <- reactive({

        selected_topic <- tutorial_set[[input$topic_select]] %>% filter(Task_List == input$func)
        index <- as.numeric(input$use_cases_row_last_clicked)
        Comments <- as.character(selected_topic[index,"Comment"])

        gen_markdown(Comments)

      })

      output$read_me <- reactive({


        val <- as.character(readme_file[1,1])
        gen_markdown(val)

      })

      session$onSessionEnded(function() {
        stopApp()
      })


    }

  )
}


#' @title Extract Tutorial datasets
#' @description Extract the datasets used by the tutorial app
#' @return A list object containing the tutorial datasets for each topic
#' @export
#'
#' @examples
#' # If you would like to extract the underlying datasets used for building the app
#' # you can invoke this function
#' show_tutorial_datasets()

show_tutorial_datasets <- function() {


  tutorial_set <- list(basic_operations = basic_operations,dplyr_tutorial = dplyr_tutorial,
                       loops_tutorial = loops_tutorial,model_tutorial = model_tutorial)
  tutorial_set

}
