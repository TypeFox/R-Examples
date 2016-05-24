library(shiny)
library(sglr)

## Define server logic required to summarize and view the selected dataset
shinyServer(function(input, output) {
    glrResult <- reactive(({
        glrSearch(p = c(input$p0, input$p1),
                  alpha = input$alpha,
                  beta = input$beta,
                  stepSize = input$stepSize,
                  tol = input$tol,
                  maxIter = input$maxIter,
                  gridIt = TRUE,
                  nGrid = input$nGrid,
                  verbose = FALSE)
    }))

    ## Show the alphaTable

    output$alphaTable <- renderTable(({
        table <- glrResult()$alphaTable
        rownames(table) <- sub("b1", "$b_1$", rownames(table))
        colnames(table) <- sub("b0", "$b_0$", colnames(table))
        table
    }))
    ## Show the betaTable

    output$betaTable <- renderTable(({
        table <- glrResult()$betaTable
        rownames(table) <- sub("b1", "$b_1$", rownames(table))
        colnames(table) <- sub("b0", "$b_0$", colnames(table))
        table
    }))

    output$plot <- renderPlot(({
        print(plotBoundary(b0 = input$bb0,
                     b1 = input$bb1,
                     p = c(input$bp0, input$bp1)))
    }))

    output$boundary <- renderTable(({
        boundary <- computeBoundary(b0 = input$bb0,
                                    b1 = input$bb1,
                                    p = c(input$bp0, input$bp1),
                                    tol = input$btol)
        d <- cbind(total=seq(length(boundary$lower)), lower=boundary$lower, upper=boundary$upper)
        storage.mode(d) <- "integer"
        colnames(d) <- c("Total No. of AEs", "Lower Boundary (Vaccine AEs)", "Upper Boundary (Vaccine AEs)")
        d
    }))

    output$estimate <- renderText(({
        boundary <- computeBoundary(b0 = input$bb0,
                                    b1 = input$bb1,
                                    p = c(input$bp0, input$bp1),
                                    tol = input$btol)

        paste("Type I Error:", format(boundary$estimate[1], digits=4),
              ", Type II Error:", format(boundary$estimate[2], digits=4))
    }))

    output$maxAE <- renderText(({
        boundary <- computeBoundary(b0 = input$bb0,
                                    b1 = input$bb1,
                                    p = c(input$bp0, input$bp1),
                                    tol = input$btol)
        paste("Maximum Total Number of AEs:", length(boundary$lower))
    }))

})
