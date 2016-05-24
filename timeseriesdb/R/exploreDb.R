#' Start a GUI to Explore Data
#' 
#' Start a graphical user interface in the user's standard web browser to search 
#' and explore time series data. Data can be searched using regular expressions for
#' keys. Hits can subsequently selected from a select box and are plotted in a joint
#' time series plot. 
#' 
#' @param con PostgreSQL Connection object
#' @param browser logical should app be fired up in the web browser? Defaults to TRUE.
#' @export 
exploreDb <- function(con, browser = T){
  
  if(!dbIsValid(con)) stop("Database connection is not valid. Can't start exploring data.")
  
  shiny::shinyApp(ui = shiny::navbarPage("timeseriesdb Data Explorer",
                           shiny::tabPanel("Build Query",
                                    shiny::fluidRow(
                                      shiny::selectInput("query_type","Select Query Type",
                                                  c("Key Based Query" = "key",
                                                    "Load Pre-Defined Set" = "set",
                                                    "Search Localized Meta Information" = "md")),
                                      shiny::uiOutput("query_builder")  
                                    )
                           ),
                           shiny::tabPanel("Plot and Export",
                                    shiny::fluidRow(
                                      shiny::column(6,shiny::tags$h2("Variable Selection"),
                                             shiny::uiOutput("choices")
                                      ),
                                      shiny::column(4,
                                             shiny::tags$h2("Store As Set"),
                                             shiny::tags$form(
                                               shiny::textInput("set_name", "Name", "")
                                               , shiny::br()
                                               , shiny::actionButton("button2", "Store the time series set"),
                                               shiny::textOutput("store_set") 
                                             )
                                             
                                      ),
                                      shiny::column(2,shiny::tags$h2("Export"),
                                             shiny::radioButtons("wide", "Use wide format?",
                                                          c("Yes" = "T",
                                                            "No" = "F")),
                                             shiny::downloadButton('download', 'Download'))
                                    ),
                                    shiny::fluidRow(
                                      shiny::column(10,shiny::plotOutput("plot")),
                                      shiny::column(2,shiny::uiOutput("legend_control"))
                                    )
                           ),
                           header = 
                             shiny::tags$style(shiny::HTML("
                                           @import url('//fonts.googleapis.com/css?family=Lato|Cabin:400,700');
                                           
                                           h2 {
                                           font-family: 'Lato';
                                           font-weight: 500;
                                           line-height: 1.1;
                                           color: #A2C3C9;
                                           }
                                           
                                           select {
                                           width:400px !important;
                                           height:150px !important;
                                           }
                                           
                                           input[type='text']{
                                           width:300px !important;
                                           }

                                           .row{
                                            margin-left:15px !important;
                                           }                                           
                                           
                                           
                                           "))
  ),
  server = function(input,output){
    # reactive stuff ----------------
    query_type <- shiny::reactive({
      out <- input$query_type
      if(out == "key"){
        class(out) <- append("key",class(out))
      } else if(out == "set") {
        class(out) <- append("set",class(out))
      } else {
        class(out) <- append(c("md","key"),class(out))
      }
      out
    })
    
    
    keys <- shiny::reactive({
      searchKeys(query_type(),input = input,con = con)
    })
    
    
    
    
    
    
    # outputs ----------------
    # flexible query builder 
    output$query_builder <- shiny::renderUI({
      createUI(query_type(),con)        
    })
    
    # display the hits 
    output$hits <- shiny::renderText({
      input$button1
      paste0(length(shiny::isolate(keys())),
             " series found. Switch to the next tab to proceed.")
    })
    
    
    # flexible choices boxes
    output$choices <- shiny::renderUI({
      createChoices(query_type(),input = input,con = con,
                    keys = keys())
      
    })
    
    # store set 
    output$store_set <- shiny::renderText({
      
      input$button2
      
      otext <- ""
      
      set_list <- shiny::isolate({
        li <- as.list(rep(input$search_type, length(input$in5)))
        names(li) <- input$in5
        li
      })
      
      if(length(set_list) > 0 && shiny::isolate(input$set_name) != "") {
        storeTsSet(con, shiny::isolate(input$set_name), set_list)
        
        otext <- paste('You have stored the set ', shiny::isolate(input$set_name), '.')
      }
      
      otext
    })
    
    
    
    
    # download handler for export
    output$download <- shiny::downloadHandler(
      filename = function(){paste0("time_series_export_",
                                   gsub(" |:|-","_",Sys.time()),".csv")}, #input$fname,
      content = function(file){
        # write.table(isolate(keys())[[1]],file)
        # don't forget to change separator
        exportTsList(shiny::isolate(keys())[input$in5],fname = file,cast = input$wide) 
      }
      
      
    )
    
    # plot that reactive so changes in selection
    output$plot <- shiny::renderPlot({
      if(is.null(input$in5)) return(NULL)
      
      li <- shiny::isolate(keys())
      li <- li[input$in5]
      class(li) <- append(class(li),"tslist")
      plot(li,use_legend = ifelse(input$legend == "yes",T,F),
           shiny_legend = T)  
    })
    
    # switch legend on/off 
    # legends are not suitable when 
    # there are to many series selected
    output$legend_control <- shiny::renderUI({
      if(is.null(input$in5)) return(NULL)
      shiny::column(2,shiny::radioButtons("legend", "Use legend?",
                            c("Yes" = "yes",
                              "No" = "no")))
      
      
    })
  }, options=list(launch.browser = browser))
}
