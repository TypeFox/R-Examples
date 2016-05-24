library(shiny)

shinyServer(function(input, output, session) {

    auth_token <- reactive({
        query <- parseQueryString(session$clientData$url_search)
        if (is.null(query$code))
            return(NULL)
        get_auth_token(creds, query$code)
    })

    output$auth_links <- renderUI({
        if (!is.null(auth_token()))
            return(NULL)
        list(
            p("Auth ", a("link", href = auth_url), "or button ", tags$a(id = "auth_link", class = "btn btn-default", href = auth_url, "Authorize app", class = "text-center"))
        )
    })

    output$accounts <- DT::renderDataTable({
        if (is.null(auth_token()))
            return(invisible(NULL))
        RGA::list_accounts(token = auth_token())
    })
})
