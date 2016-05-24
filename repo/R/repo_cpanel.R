#' Open a visual interface to the repo
#'
#' Opens a browser window with a Shiny interface to a repo. The
#' interface is preliminary and has some exploration features together
#' with a "Load into workspace" button for a selected item.
#'
#' @param reporoot An object of class repo. Can be NULL like for repo_open.
#' @param env Environment to export variables to. Defaults to globalenv.
#' @return Used for side effects.

repo_cpanel <- function(reporoot=NULL, env=globalenv()){

    if(is.null(reporoot))
        repo <- repo_open() else repo <- repo_open(reporoot)
    
    items <- sapply(repo$entries(), get, x="name")    
    
    shiny::shinyApp(
        server = function(input, output)
        {
            checkItemSelected <- function() {
                if(is.null(input$items))
                    return(FALSE)
                return(TRUE)
            }
                output$items <- shiny::renderUI({
                    shiny::selectizeInput("items", "Select item", choices=items, selected=1)
                })
                output$itemDetails_name <- shiny::renderText({
                    if(checkItemSelected())
                        repo$entries()[[which(items==input$items)]]$name
                })
                output$itemDetails_description <- shiny::renderText({
                    if(checkItemSelected())
                        repo$entries()[[which(items==input$items)]]$description
                })
                output$itemDetails_tags <- shiny::renderText({
                    if(checkItemSelected())
                    paste(repo$entries()[[which(items==input$items)]]$tags,
                          collapse=", ")
                })
                output$itemDetails_dimensions <- shiny::renderText({
                    if(checkItemSelected())
                        paste(repo$entries()[[which(items==input$items)]]$dims,
                              collapse="x")
                })
                output$itemDetails_size <- shiny::renderText({
                    if(checkItemSelected())
                        paste(repo$entries()[[which(items==input$items)]]$size)
                })
                output$itemDetails_time <- shiny::renderText({
                    if(checkItemSelected())
                        paste(repo$entries()[[which(items==input$items)]]$timestamp)
                })
                output$itemDetails_source <- shiny::renderText({
                    if(checkItemSelected())
                        paste(repo$entries()[[which(items==input$items)]]$source,
                          collapse=", ")
                })
                output$itemDetails_attached <- shiny::renderText({
                    if(checkItemSelected())
                        paste(repo$entries()[[which(items==input$items)]]$attachedto,
                              collapse=", ")
                })
                output$itemDetails_depends <- shiny::renderText({
                    if(checkItemSelected())
                        paste(repo$entries()[[which(items==input$items)]]$depends,
                              collapse=", ")
                })
                output$itemDetails_md5 <- shiny::renderText({
                    if(checkItemSelected())
                        paste(repo$entries()[[which(items==input$items)]]$checksum,
                              collapse=", ")
                })
                output$itemDetails_path <- shiny::renderText({
                    if(checkItemSelected())
                        paste(repo$entries()[[which(items==input$items)]]$dump)
                })

                output$pies <- shiny::renderPlot({
                    repo$pies()
                })

                shiny::observeEvent(input$loadButton, {
                    if(checkItemSelected())
                        assign(input$items, repo$get(input$items), envir=env)
                })

            ## observeEvent(input$updateItemButton, {
            ##     if(checkItemSelected())
            ##     {
            ##         selitem <- which(items==input$items)
            ##         repo$set(input$items, newname=input$items,
            ##                  description = )
            ##     }
            ## })
            },

        
        ui=shiny::fluidPage(

            shiny::titlePanel("Repo Control Panel"),

            shiny::sidebarLayout(
                shiny::sidebarPanel(
                    shiny::h6(paste0("Current repo: ", repo$root())),
                    shiny::uiOutput("items"),
                    shiny::actionButton("loadButton", "Load into workspace"),
                    shiny::h4("Items size overview"),
                    shiny::plotOutput("pies")
                    ),

                shiny::mainPanel(
                    shiny::h2("Item details"),
                    shiny::verbatimTextOutput("itemDetails_name"),
                    shiny::tags$p("Description:"),
                    shiny::verbatimTextOutput("itemDetails_description"),
                    shiny::tags$p("Tags:"),
                    shiny::verbatimTextOutput("itemDetails_tags"),
                    shiny::tags$p("Dimensions:"),
                    shiny::verbatimTextOutput("itemDetails_dimensions"),
                    shiny::tags$p("Size on disk:"),
                    shiny::verbatimTextOutput("itemDetails_size"),
                    shiny::tags$p("Timestamp:"),
                    shiny::verbatimTextOutput("itemDetails_time"),
                    shiny::tags$p("Source:"),
                    shiny::verbatimTextOutput("itemDetails_source"),                                        
                    shiny::tags$p("Attached to:"),
                    shiny::verbatimTextOutput("itemDetails_attached"),
                    shiny::tags$p("Depends on:"),
                    shiny::verbatimTextOutput("itemDetails_depends"),
                    shiny::tags$p("MD5 checksum:"),
                    shiny::verbatimTextOutput("itemDetails_md5"),
                    shiny::tags$p("Stored in:"),
                    shiny::verbatimTextOutput("itemDetails_path")
                    )
                )
            )
        )
}
