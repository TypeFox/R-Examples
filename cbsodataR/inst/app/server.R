library(shiny)
library(cbsodata)

shinyServer(function(input, output, session) {
  tables <- get_tables(Language="nl", select=c("Title", "Summary", "Identifier"))
  tables$ShortDescription <- NULL
  choices = tables$Identifier
  names(choices) = tables$Title
  options = data.frame(label=tables$Title, value=tables$Identifier)
  updateSelectizeInput(session, "table", choices=options, server = T,
                       options=list(placeholder="Select table"))
  
  output$table_list <- renderDataTable(tables)
})