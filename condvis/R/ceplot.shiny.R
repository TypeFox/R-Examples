ceplot.shiny <-
## this code really needs a rewrite!
function(data, model, response = NULL, S = NULL, C = NULL, cex.axis = NULL, 
    cex.lab = NULL, tck = NULL, Corder = "default")
{
    ui <- NULL
    server <- NULL
    data <- na.omit(data)
    if(!requireNamespace("shiny", quietly = TRUE))
        stop("requires package 'shiny'")
    else if (!exists("shinyApp")) attachNamespace("shiny")    
    model <- if (!identical(class(model), "list"))
        list(model)
    else model
    model.name <- if (!is.null(names(model)))
        names(model)
    else NULL
    varnamestry <- try(getvarnames(model[[1]]), silent = TRUE)
    response <- if (is.null(response))
        if (class(varnamestry) != "try-error")
           which(colnames(data) == varnamestry$response[1])
        else stop("could not extract response from 'model'.")
    else if (is.character(response))
            which(colnames(data) == response)
        else response
    S <- if(is.null(S)){
         (1:ncol(data))[-response][1L]
        } else if (is.character(S))
            vapply(S, function(x) which(colnames(data) == x), numeric(1))
            else S
    C <- if (is.null(C))
        arrangeC(data[, -c(response, S)])
    else C
    try(
        if (class(varnamestry) != "try-error"){
            possibleC <- unique(unlist(lapply(lapply(model, getvarnames), `[[`, 
                2)))
            possibleC <- possibleC[possibleC %in% colnames(data)]
            C <- arrangeC(data[, possibleC[!(possibleC %in% colnames(data)[S])], 
                drop = FALSE], method = Corder)
        }     
    , silent = TRUE)
    C <- if (all(vapply(C, is.numeric, logical(1))))
        as.list(C)
    else if (all(vapply(C, is.character, logical(1))))
            lapply(C, match, table = colnames(data))
        else
            stop("'C' should be a vector or list (containing vectors of length",
                 " 1 or 2) with integer column indices or character variable",
                 " names from 'data'.")
    uniqC <- unique(unlist(C))
    n.selector.cols <- ceiling(length(C) / 4L)
    if (any(response %in% uniqC))
        stop("cannot have 'response' variable in 'C'")
    if (any(response %in% S))
        stop("cannot have 'response' variable in 'S'")
    if (!identical(length(unique(vapply(lapply(model, getvarnames),
        `[[`, character(1), 1))), 1L))
        stop("cannot compare models with different response variables")
    if (!identical(length(intersect(S, uniqC)), 0L))
        stop("cannot have variables common to both 'S' and 'C'")
    Xc.cond <- data[1, uniqC, drop = FALSE]
    rownames(Xc.cond) <- "Xc.cond"
    tmp <- new.env()
    assign("Xc.cond", Xc.cond, tmp)
    Xc <- data[, uniqC, drop = FALSE]
    xcplotsize <- 190
    contplot <- !is.factor(data[, response]) & !any(vapply(data[, S, drop = FALSE], is.factor, logical(1))) & identical(length(S), 2L)
    eval(parse(text = paste("
        ui <- fluidPage(
            fluidRow(
                column(5
                    , fluidRow(
                        column(8, 
                            if (contplot) {
                                tabsetPanel(
                                    tabPanel('Contour', plotOutput('plotS', height = '100%', width = '80%'), value = 1),
                                    tabPanel('Perspective', plotOutput('plotS2', height = '100%', width = '80%'), value = 2)
                                , id = 'tab')
                            } else plotOutput('plotS', height = '100%', width = '80%')
                        ),
                        column(3, if (identical(length(S), 2L)) plotOutput('legend', height = 400, width = '20%'))
                    )
                    , actionButton('saveButton', 'Take snapshot (pdf)')
                    , conditionalPanel(condition = 'input.tab == 2', numericInput('phi', 'Vertical rotation: ', 20, -180, 180))
                    , conditionalPanel(condition = 'input.tab == 2', numericInput('theta', 'Horizontal rotation: ', 45, -180, 180))
                    , sliderInput('sigma', 'Weighting function parameter: ', 0.01, 5, step = 0.01, value = 1)
                    , radioButtons('type', 'Weighting function type:', c('euclidean', 'maxnorm'))
                ),
                column(7,
                    fluidRow( helpText(strong('Condition selector plots')) ),
                    fluidRow(
                    column(3,",
                        paste(" if (length(C) >= ", 1:4, ") {plotOutput('plotC", 1:4, "', height = ", xcplotsize,", width = ", xcplotsize,", click = 'plotC", 1:4, "click') }", sep = "", collapse = ",\n")
                    ,"
                    )
                    , column(3,",
                        paste("if (length(C) >= ", 5:8, ") {plotOutput('plotC", 5:8, "', height = ", xcplotsize,", width = ", xcplotsize,", click = 'plotC", 5:8, "click') }", sep = "", collapse = ",\n")
                    ,"
                    )
                    , column(3,",
                        paste("if (length(C) >= ", 9:12, ") {plotOutput('plotC", 9:12, "', height = ", xcplotsize,", width = ", xcplotsize,", click = 'plotC", 9:12, "click') }", sep = "", collapse = ",\n")
                    ,"
                    )  
                    , column(3,",
                        paste("if (length(C) >= ", 13:16, ") {plotOutput('plotC", 13:16, "', height = ", xcplotsize,", width = ", xcplotsize,", click = 'plotC", 13:16, "click') }", sep = "", collapse = ",\n")
                    ,"
                    ) 
                    )
                    , fluidRow( helpText(strong('Current condition/section')) )
                    , fluidRow( tableOutput('text') )
                )               
            )
        )
    ")))
    
    eval(parse(text = paste("
        server <-
        function (input, output){
            output$text <- renderTable({
                Xc.cond <- get('Xc.cond', envir = tmp)", paste("
            if (!is.null(input$plotC", 1:length(C), "click$x)){
                arefactors <- unlist(lapply(data[, C[[", 1:length(C), "]], drop = FALSE], is.factor))
                if (identical(length(arefactors), 1L)){
                    if(arefactors)
                        Xc.cond[, names(data)[C[[", 1:length(C), "]]]] <- factor(levels(data[, names(data)[C[[", 1:length(C), "]]]])[which.min(abs(input$plotC", 1:length(C), "click$x - (1:length(levels(data[, names(data)[C[[", 1:length(C), "]]]])))))], levels = levels(data[, names(data)[C[[", 1:length(C), "]]]]))
                    else Xc.cond[, names(data)[C[[", 1:length(C), "]]]] <- input$plotC", 1:length(C), "click$x
                }
                if (identical(length(arefactors), 2L)){
                    if (all(arefactors)){
            
                    } else if (any(arefactors)){
                       Xc.cond[, names(data)[C[[", 1:length(C), "]]][which(arefactors)]] <- factor(levels(data[, names(data)[C[[", 1:length(C), "]]][which(arefactors)]])[which.min(abs(input$plotC", 1:length(C), "click$x - (1:length(levels(data[, names(data)[C[[", 1:length(C), "]]][which(arefactors)]])))))], levels = levels(data[, names(data)[C[[", 1:length(C), "]]][which(arefactors)]]))
                        Xc.cond[, names(data)[C[[", 1:length(C), "]]][which(!arefactors)]] <- input$plotC", 1:length(C), "click$y
                        } else {
                            Xc.cond[, names(data)[C[[", 1:length(C), "]]][1]] <- input$plotC", 1:length(C), "click$x
                            Xc.cond[, names(data)[C[[", 1:length(C), "]]][2]] <- input$plotC", 1:length(C), "click$y
                        }
            
                }
                assign('Xc.cond', Xc.cond, envir = tmp) 
            }
            ", sep = "", collapse = ""),"
            colnames(Xc.cond) <- vapply(colnames(Xc.cond), substr, character(1), 1, 3)
            Xc.cond
            })
            output$legend <- renderPlot({
                xslegend(y = data[, response], name = colnames(data)[response])
            }, width = 100, height = 400)
            output$plotS <- renderPlot({
            Xc.cond <- get('Xc.cond', envir = tmp)", paste("
            if (!is.null(input$plotC", 1:length(C), "click$x)){
                arefactors <- unlist(lapply(data[, C[[", 1:length(C), "]], drop = FALSE], is.factor))
                if (identical(length(arefactors), 1L)){
                    if(arefactors)
                        Xc.cond[, names(data)[C[[", 1:length(C), "]]]] <- factor(levels(data[, names(data)[C[[", 1:length(C), "]]]])[which.min(abs(input$plotC", 1:length(C), "click$x - (1:length(levels(data[, names(data)[C[[", 1:length(C), "]]]])))))], levels = levels(data[, names(data)[C[[", 1:length(C), "]]]]))
                    else Xc.cond[, names(data)[C[[", 1:length(C), "]]]] <- input$plotC", 1:length(C), "click$x
                }
                if (identical(length(arefactors), 2L)){
                    if (all(arefactors)){
            
                    } else if (any(arefactors)){
                       Xc.cond[, names(data)[C[[", 1:length(C), "]]][which(arefactors)]] <- factor(levels(data[, names(data)[C[[", 1:length(C), "]]][which(arefactors)]])[which.min(abs(input$plotC", 1:length(C), "click$x - (1:length(levels(data[, names(data)[C[[", 1:length(C), "]]][which(arefactors)]])))))], levels = levels(data[, names(data)[C[[", 1:length(C), "]]][which(arefactors)]]))
                        Xc.cond[, names(data)[C[[", 1:length(C), "]]][which(!arefactors)]] <- input$plotC", 1:length(C), "click$y
                        } else {
                            Xc.cond[, names(data)[C[[", 1:length(C), "]]][1]] <- input$plotC", 1:length(C), "click$x
                            Xc.cond[, names(data)[C[[", 1:length(C), "]]][2]] <- input$plotC", 1:length(C), "click$y
                        }
            
                }
                assign('Xc.cond', Xc.cond, envir = tmp) 
            }
            ", sep = "", collapse = ""),"
                    vw <- visualweight(xc = Xc, xc.cond = get('Xc.cond', envir = tmp), sigma = input$sigma, distance = input$type)
                    k <- vw$k
                    data.colour <- rgb(1 - k, 1 - k, 1 - k)
                    data.order <- vw$order
                    plotxsobject <- plotxs1(xs = data[, S, drop = FALSE],
                        y = data[, response, drop = FALSE], xc.cond = get('Xc.cond', envir = tmp), model = model,
                        model.colour = NULL, model.lwd = NULL, model.lty = NULL,
                        model.name = model.name, yhat = NULL, mar = NULL,
                        data.colour = data.colour, data.order = data.order, view3d = FALSE)
            }, width = 400, height = 400)
            output$plotS2 <- renderPlot({
            Xc.cond <- get('Xc.cond', envir = tmp)", paste("
            if (!is.null(input$plotC", 1:length(C), "click$x)){
                arefactors <- unlist(lapply(data[, C[[", 1:length(C), "]], drop = FALSE], is.factor))
                if (identical(length(arefactors), 1L)){
                    if(arefactors)
                        Xc.cond[, names(data)[C[[", 1:length(C), "]]]] <- factor(levels(data[, names(data)[C[[", 1:length(C), "]]]])[which.min(abs(input$plotC", 1:length(C), "click$x - (1:length(levels(data[, names(data)[C[[", 1:length(C), "]]]])))))], levels = levels(data[, names(data)[C[[", 1:length(C), "]]]]))
                    else Xc.cond[, names(data)[C[[", 1:length(C), "]]]] <- input$plotC", 1:length(C), "click$x
                }
                if (identical(length(arefactors), 2L)){
                    if (all(arefactors)){
            
                    } else if (any(arefactors)){
                       Xc.cond[, names(data)[C[[", 1:length(C), "]]][which(arefactors)]] <- factor(levels(data[, names(data)[C[[", 1:length(C), "]]][which(arefactors)]])[which.min(abs(input$plotC", 1:length(C), "click$x - (1:length(levels(data[, names(data)[C[[", 1:length(C), "]]][which(arefactors)]])))))], levels = levels(data[, names(data)[C[[", 1:length(C), "]]][which(arefactors)]]))
                        Xc.cond[, names(data)[C[[", 1:length(C), "]]][which(!arefactors)]] <- input$plotC", 1:length(C), "click$y
                        } else {
                            Xc.cond[, names(data)[C[[", 1:length(C), "]]][1]] <- input$plotC", 1:length(C), "click$x
                            Xc.cond[, names(data)[C[[", 1:length(C), "]]][2]] <- input$plotC", 1:length(C), "click$y
                        }
            
                }
                assign('Xc.cond', Xc.cond, envir = tmp) 
            }
            ", sep = "", collapse = ""),"
                    vw <- visualweight(xc = Xc, xc.cond = get('Xc.cond', envir = tmp), sigma = input$sigma, distance = input$type)
                    k <- vw$k
                    data.colour <- rgb(1 - k, 1 - k, 1 - k)
                    data.order <- vw$order
                    plotxsobject <- plotxs.shiny(xs = data[, S, drop = FALSE],
                        y = data[, response, drop = FALSE], xc.cond = get('Xc.cond', envir = tmp), model = model,
                        model.colour = NULL, model.lwd = NULL, model.lty = NULL,
                        model.name = model.name, yhat = NULL, mar = NULL,
                        data.colour = data.colour, data.order = data.order, view3d = TRUE, phi3d = input$phi, theta3d = input$theta)
            }, width = 400, height = 400)
            ", paste("
            output$plotC", 1:length(C), " <- renderPlot({
                o <- plotxc(xc = data[, C[[", 1:length(C), "]]], 
                    xc.cond = get('Xc.cond', envir = tmp)[, names(data)[C[[", 1:length(C), "]]]],
                    name = colnames(data)[C[[", 1:length(C), "]]],
                    select.colour = 'blue', select.lwd = 2, cex.axis = cex.axis,
                    cex.lab = cex.lab, tck = tck, shiny = TRUE)            
                Xc.cond <- get('Xc.cond', envir = tmp)
            if (!is.null(input$plotC", 1:length(C), "click$x)){
                arefactors <- unlist(lapply(data[, C[[", 1:length(C), "]], drop = FALSE], is.factor))
                if (identical(length(arefactors), 1L)){
                    if(arefactors)
                        Xc.cond[, names(data)[C[[", 1:length(C), "]]]] <- factor(levels(data[, names(data)[C[[", 1:length(C), "]]]])[which.min(abs(input$plotC", 1:length(C), "click$x - (1:length(levels(data[, names(data)[C[[", 1:length(C), "]]]])))))], levels = levels(data[, names(data)[C[[", 1:length(C), "]]]]))
                    else Xc.cond[, names(data)[C[[", 1:length(C), "]]]] <- input$plotC", 1:length(C), "click$x
                }
                if (identical(length(arefactors), 2L)){
                    if (all(arefactors)){
                        varnames <- o$name
                        sptmp <- o$sptmp
                        rectcoords <- data.frame(sptmp$xleft, sptmp$xright, 
                            sptmp$ybottom, sptmp$ytop)
                        if (c(input$plotC", 1:length(C), "click$x, input$plotC", 1:length(C), "click$y) %inrectangle% 
                            c(min(sptmp$xleft), max(sptmp$xright) ,
                            min(sptmp$ybottom), max(sptmp$ytop)) ){
                                comb.index <- apply(rectcoords, 1L, 
                                    `%inrectangle%`, point = c(input$plotC", 1:length(C), "click$x, input$plotC", 1:length(C), "click$y))
                                if (any(comb.index)){
                                    Xc.cond[, names(data)[C[[", 1:length(C), "]]][1]] <- as.factor(sptmp$xnames)[comb.index]
                                    Xc.cond[, names(data)[C[[", 1:length(C), "]]][2]] <- as.factor(sptmp$ynames)[comb.index]
                                }     
                            }
                    } else if (any(arefactors)){
                       Xc.cond[, names(data)[C[[", 1:length(C), "]]][which(arefactors)]] <- factor(levels(data[, names(data)[C[[", 1:length(C), "]]][which(arefactors)]])[which.min(abs(input$plotC", 1:length(C), "click$x - (1:length(levels(data[, names(data)[C[[", 1:length(C), "]]][which(arefactors)]])))))], levels = levels(data[, names(data)[C[[", 1:length(C), "]]][which(arefactors)]]))
                        Xc.cond[, names(data)[C[[", 1:length(C), "]]][which(!arefactors)]] <- input$plotC", 1:length(C), "click$y
                        } else {
                            Xc.cond[, names(data)[C[[", 1:length(C), "]]][1]] <- input$plotC", 1:length(C), "click$x
                            Xc.cond[, names(data)[C[[", 1:length(C), "]]][2]] <- input$plotC", 1:length(C), "click$y
                        }
            
                }
                assign('Xc.cond', Xc.cond, envir = tmp) 
            }            
            
                plotxc(xc = data[, C[[", 1:length(C), "]]], 
                    xc.cond = get('Xc.cond', envir = tmp)[, names(data)[C[[", 1:length(C), "]]]],
                    name = colnames(data)[C[[", 1:length(C), "]]],
                    select.colour = 'blue', select.lwd = 2, cex.axis = cex.axis,
                    cex.lab = cex.lab, tck = tck, shiny = TRUE)
            })
            ", sep = "", collapse = "\n"),"
            observeEvent(input$saveButton, {
                n.selector.cols <- ceiling(length(C) / 4L)
                select.colwidth <- max(min(0.18 * n.selector.cols, 0.45), 0.2)  
                width <- 8.5 + 2 * n.selector.cols 
                filename <- paste('snapshot_', gsub(':', '.', gsub(' ', '_', Sys.time())), '.pdf', sep = '') 
                pdf(file = filename, width = width, height = 8)
                ceplot.static(data = data, model = model, response = response, 
                    S = S, C = C, sigma = input$sigma, distance = input$type, 
                    cex.axis = cex.axis, cex.lab = cex.lab, tck = tck, 
                    view3d = if (!is.null(input$tab)) input$tab == 2 else FALSE, Xc.cond = get('Xc.cond', envir = tmp), 
                    theta3d = if (!is.null(input$theta)) input$theta else NULL, phi3d = if (!is.null(input$phi)) input$phi else NULL)
                dev.off()
            })
        }
    ")))


    shiny::shinyApp(ui, server)
}
