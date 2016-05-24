#' Polymorphic Helpers for the Shiny Based GUI
#' 
#' Though R is not strictly OO it is sometimes useful to add some object orientation. 
#' In the case of a shiny app with an dynamically created user interface a bit
#' of pseudo polymorphism helps to avoid too many but,ifs and elses. 
#' Depending on the query type that is selected by the user a UI creation
#' method is called. 
#' 
#' @param x character query type.
#' @param con PostgreSQL Connection object.
#' @param ... additional arguments to be passed to the respective methods.
#' @rdname createUI
createUI <- function(x, con, ...) UseMethod("createUI")

#' @rdname createUI
createUI.key <- function(x, con, ...){
  
  md_keys <- DBI::dbGetQuery(con,"SELECT DISTINCT k FROM 
                        (SELECT skeys(meta_data) as k 
                        FROM meta_data_unlocalized) as dt;")$k
  
  names(md_keys) <- md_keys
  st_keys <- c(c("Main ts_key" = "ts_key"),md_keys)
  
  shiny::column(6,
         shiny::radioButtons("search_type", "Key Based Query",
                      st_keys),
         shiny::tags$form(
           shiny::textInput("key", "Search for Key", "")
           , shiny::br()
           , shiny::actionButton("button1", "Search timeseriesdb")
         ),
         shiny::textOutput("hits")
  )    
}

#' @rdname createUI
createUI.set <- function(x, con, ...){
  sets <- listTsSets(con)
  shiny::column(6,
         shiny::selectInput('button1',
                     "Select a Set of Time Series",
                     sets,
                     multiple = T, selectize=FALSE)
  )
  
  
}

#' @rdname createUI
createUI.md <- function(x,con ,...){
  md_keys <- DBI::dbGetQuery(con,"SELECT DISTINCT k FROM 
                        (SELECT skeys(meta_data) as k 
                        FROM meta_data_localized) as dt;")$k
  
  names(md_keys) <- md_keys
  st_keys <- md_keys
  
  shiny::column(6,
         shiny::radioButtons("search_type", "Key Based Query",
                      st_keys),
         shiny::tags$form(
           shiny::textInput("key", "Search for Key", "")
           , shiny::br()
           , shiny::actionButton("button1", "Search timeseriesdb")
         ),
         shiny::textOutput("hits")
  )
  
  
}





#' Create Choices Boxes Dynamically
#' 
#' Dynamically create choices boxes depending on interactive choices by the user. 
#' The UI reacts to user input. We use a somewhat polymorphic approach here to avoid
#' too much if/else and improve readabitity. 
#'
#' @param x query type as character. 
#' @param input input that triggers action. Defaults to NULL. 
#' @param keys results to rendered to the choice box. Defaults to NULL. 
#' @param con PostgreSQL Connection object.
#' @param ... addtional arguments to be passed to the methods.
#' @rdname createChoices
createChoices <- function(x, input = NULL, keys = NULL, con,  ...) UseMethod("createChoices")

#' @rdname createChoices
createChoices.key <- function(x, input = NULL, keys = NULL,...){
  
  input$button1
  if(length(keys) != 0){
    shiny::selectInput('in5', paste0('Select keys (',
                              length(shiny::isolate(keys)),' hits)'),
                names(shiny::isolate(keys)),
                multiple = T, selectize=FALSE)    
  } else {
    NULL
  }
}

#' @rdname createChoices
createChoices.md <- function(x,input = NULL,keys = NULL, con, ...){
NextMethod()  
  
}





#' @rdname createChoices
createChoices.set <- function(x,input = NULL,keys = NULL, con, ...){
  input$get_series
  if(length(shiny::isolate(input$button1)) != 0){
    ts_keys <- loadTsSet(con,input$button1)$keys
    shiny::selectInput('in5', paste0("Time Series in set '",input$button1,"'"),
                ts_keys,
                multiple = T, selectize=FALSE)    
  } else {
    NULL
  }
}



#' Search Keys From Pre-Defined Or Key Based Queries
#' @param x query type as character. 
#' @param input input that triggers action. Defaults to NULL. 
#' @param con PostgreSQL Connection object.
#' @param ... addtional arguments to be passed to the methods.
#' @rdname searchKeys
searchKeys <- function(x,input = NULL, con ,...) UseMethod("searchKeys")

#' @rdname searchKeys
searchKeys.key <- function(x,input = NULL, con, ...) {
  if(input$key != ""){
    if(input$search_type == "ts_key"){
      keys <- con %k% input$key # double check this
      keys  
    } else{
      "%m%" <- createMetaDataHandle(input$search_type,keep_keys = T)
      keys <- con %m% input$key
      keys
    }
  } else {
    NULL
  }
}

#' @rdname searchKeys
searchKeys.md <- function(x,input = NULL, con, ...) {
  if(input$key != ""){
    "%m%" <- createMetaDataHandle(input$search_type,keep_keys = T,
                                  tbl = "meta_data_localized")
    keys <- con %m% input$key
    keys
  } else {
    NULL
  }
}



#' @rdname searchKeys
searchKeys.set <- function(x,input = NULL, con, ...){
  # it's clear we're looking for a set named input$button1
  # for the current user, could improve this using tryCatch
  # instead of if
  if(!is.null(input$button1)){
    keys_chars <- loadTsSet(con,input$button1)$keys
    keys <- readTimeSeries(keys_chars,con)
    keys  
  } else {
    NULL
  }
  
}
