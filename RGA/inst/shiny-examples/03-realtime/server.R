library(shiny)

shinyServer(function(input, output, session) {

    auth_token <- reactive({
        query <- parseQueryString(session$clientData$url_search)
        if (is.null(query$code))
            return(NULL)
        get_auth_token(creds, query$code)
    })

    output$auth_button <- renderUI({
        if (!is.null(auth_token()))
            return(NULL)
        tags$a(id = "auth_link", class = "btn btn-default", href = auth_url, "Authorize app", class = "text-center")
    })

    output$active_users <- renderText({
        if (is.null(auth_token()))
            return(invisible(NULL))
        invalidateLater(delay, session)
        paste(RGA::get_realtime(id, metrics = "rt:activeUsers", token = auth_token())$active.users)
    })
})
