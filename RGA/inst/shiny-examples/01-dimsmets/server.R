library(shiny)

shinyServer(function(input, output) {
    output$table <- DT::renderDataTable({
        if (input$group != "All")
            ga <- ga[ga$group == input$group,]
        if (input$type != "All")
            ga <- ga[ga$type == input$type,]
        if (!input$status)
            ga <- ga[ga$status != "DEPRECATED",]
        if (input$calc)
            ga <- ga[!is.na(ga$calculation),]
        if (input$segments)
            ga <- ga[ga$allowedInSegments,]
        if (!is.null(input$columns))
            selected <- c(selected, input$columns)
        ga[, selected, drop = FALSE]
    },
        rownames = FALSE,
        selection = "none", filter = "top", style = "bootstrap",
        extensions = c("Responsive"),
        options = list(
            autoWidth = TRUE, paging = FALSE, searchHighlight = TRUE,
            dom = 'irt'
        )
    )
})
