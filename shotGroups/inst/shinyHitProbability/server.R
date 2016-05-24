source("helper.R")

shinyServer(function(input, output) {
    #####-----------------------------------------------------------------------
    ## provide the data - reactive conductor
    #####-----------------------------------------------------------------------
    coords <- reactive({
        ## only change when explicitly applied
        input$applyData

        ## isolate against non-applied changes in data input UI elements
        isolate({
            if(input$datIn == '1') {
                ## built in data
                get(dataBuiltIn[input$builtInData])
            } else if(input$datIn == '2') {
                ## upload files
                if(!is.null(input$fileUpload)) {
                    fPath <- input$fileUpload$datapath

                    if(input$fileType == '1') {          ## OnTarget 1.*
                        readDataOT1(fPath=dirname(fPath),  fNames=basename(fPath))
                    } else if(input$fileType == '2') {   ## OnTarget 2.*, 3.*
                        readDataOT2(fPath=dirname(fPath),  fNames=basename(fPath))
                    } else if(input$fileType == '3') {   ## other
                        readDataMisc(fPath=dirname(fPath), fNames=basename(fPath))
                    }
                } else {
                    NULL
                }
            } else if(input$datIn == '3') {
                ## paste data
                fPath <- tempfile()
                writeLines(input$datPaste, fPath)
                if(input$fileType == '1') {          ## OnTarget 1.*
                    readDataOT1(dirname(fPath),  fNames=basename(fPath))
                } else if(input$fileType == '2') {   ## OnTarget 2.*, 3.*
                    readDataOT2(dirname(fPath),  fNames=basename(fPath))
                } else if(input$fileType == '3') {   ## other
                    readDataMisc(dirname(fPath), fNames=basename(fPath))
                }
            } else {
                NULL
            }
        })
    })

    #####---------------------------------------------------------------------------
    ## provide file information - UI element
    #####---------------------------------------------------------------------------
    output$fileInfo <- renderUI({
        xy <- coords()
        dstTrgt <- if(!is.null(xy$distance) && !all(is.na(xy$distance))) {
            xy$distance[1]
        } else {
            "not available"
        }

        nGroups <- if(!is.null(xy$series) && !all(is.na(xy$series))) {
            nlevels(xy$series)
        } else {
            nlevels(xy$group)
        }

        comment <- attributes(xy)$comment
        ammo <- if(!is.null(xy$ammunition) && !all(is.na(xy$ammunition)) && !all(xy$ammunition == "")) {
            paste(unique(xy$ammunition), collapse=", ")
        } else {
            NULL
        }

        ## isolate against non-applied changes in data input UI elements
        isolate({
            x <- if(input$datIn == "1") {
                paste0("<p>For details on this data set (measurement units etc.), see ",
                       "<a href='http://www.rdocumentation.org/packages/shotGroups/functions/",
                       dataBuiltIn[input$builtInData], "'>", dataBuiltIn[input$builtInData],
                       "</a></p>Name: ", dataBuiltIn[input$builtInData], "<br />")
            } else if(input$datIn == "2") {
                if(!is.null(input$fileUpload)) {
                    paste0("<p>Name: ", paste(basename(input$fileUpload$name), collapse=", "), "<br />")
                } else {
                    "<p>"
                }
            } else {
                paste0("<p>Pasted data<br />")
            }

            y <- paste0(x, "Number of shots: ", nrow(xy),
                   "<br />Distance to target ", dstTrgt,
                   "<br />Number of groups: ", nGroups,
                   ifelse(!is.null(comment), "<br />Additional information: ", ""), comment,
                   ifelse(!is.null(ammo),    "<br />Ammunition: ", ""), ammo, "</p>")
            HTML(y)
        })
    })

    #####---------------------------------------------------------------------------
    ## provide short file information - UI element
    #####---------------------------------------------------------------------------
    output$fileInfoShort <- renderUI({
        xy <- coords()

        ## isolate against non-applied changes in data input UI elements
        isolate({
            nGroups <- if(!is.null(xy$series) && !all(is.na(xy$series))) {
                nlevels(xy$series)
            } else {
                nlevels(xy$group)
            }

            x <- if(input$datIn == "1") {
                paste0("<p>Name: ", dataBuiltIn[input$builtInData], "<br />")
            } else if((input$datIn == "2") && !is.null(input$fileUpload)) {
                paste0("<p>Name: ", paste(basename(input$fileUpload$name), collapse=", "), "<br />")
            } else {
                paste0("<p>Pasted data<br />")
            }

            y <- paste0(x, "Number of shots: ", nrow(xy),
                        "<br />Number of groups: ", nGroups, "</p>")
            HTML(y)
        })
    })

    #####---------------------------------------------------------------------------
    ## provide file name + ammo information - reactive conductor
    #####---------------------------------------------------------------------------
    fileName <- reactive({
        xy <- coords()
        ammo <- if(!is.null(xy$ammunition) && !all(is.na(xy$ammunition)) && !all(xy$ammunition == "")) {
            paste(unique(xy$ammunition), collapse=", ")
        } else {
            NULL
        }

        ## isolate against non-applied changes in data input UI elements
        isolate({
            fName <- if(input$datIn == "1") {
                dataBuiltIn[input$builtInData]
            } else if(input$datIn == "2") {
                if(!is.null(input$fileUpload)) {
                    paste(basename(input$fileUpload$name), collapse=", ")
                } else {
                    NULL
                }
            } else {
                paste0("Pasted data")
            }
            c(fName, ammo, "")
        })
    })

    #####---------------------------------------------------------------------------
    ## distance to target, unit distance, unit xy-coords - UI element
    #####---------------------------------------------------------------------------
    output$unitDstXY <- renderUI({
        xy <- coords()
        dstTarget <- if(!is.null(xy$distance) && !all(is.na(xy$distance)) &&
                        !all(xy$distance == "")) {
            ## distance to target is given in input data
            xy$distance[1]
        } else {
            ## default distance to target
            100
        }

        ## dst to target, unit dst, unit xy
        inputPanel(numericInput("dstTrgt", h5("Distance to target"),
                                min=0, step=1, value=dstTarget),
                   selectInput("unitDst", h5("Measurement unit distance"),
                               choices=unitsDst, selected=2),
                   selectInput("unitXY", h5("Measurement unit coordinates"),
                               choices=unitsXY, selected=3))
    })

    #####---------------------------------------------------------------------------
    ## string for conversion argument - unit distance to unit xy-coords
    ## reactive conductor
    #####---------------------------------------------------------------------------
    conversionStr <- reactive({
        paste0(unitsDstInv[input$unitDst], "2",
               unitsXYInv[input$unitXY], collapse="")
    })

    #####---------------------------------------------------------------------------
    ## hit probability -> radius
    #####---------------------------------------------------------------------------
    ## CEP output list - reactive conductor
    CEPListRadius <- reactive({
        xy <- coords()
        if(!is.null(xy)) {
            ## if no CEP type is selected -> fall back to default CorrNormal
            CEPtype <- if(!is.null(input$hitpCEPtype1)) {
                CEPtypesInv[input$hitpCEPtype1]
            } else {
                "CorrNormal"
            }

            ## hit probability -> radius
            x <- getCEP(xy,
                        CEPlevel=input$hitpLevel,
                        dstTarget=input$dstTrgt,
                        conversion=conversionStr(),
                        accuracy=input$hitpAcc,
                        type=CEPtype,
                        doRob=input$hitpDoRob1)$CEP[[1]]
            setNames(list(x), paste0("CEP_", 100*input$hitpLevel, "%"))
        } else {
            NULL
        }
    })

    ## confidence ellipse output list - reactive conductor
    confEllList <- reactive({
        xy <- coords()
        if(!is.null(xy)) {
            x <- getConfEll(xy,
                            level=input$hitpLevel,
                            dstTarget=input$dstTrgt,
                            conversion=conversionStr(),
                            doRob=input$hitpDoRob1)

            x <- if(input$hitpDoRob1) {
                list(confEllSizeRobust=x$sizeRob, confEllShapeRobust=x$shapeRob)
            } else {
                list(confEllSize=x$size, confEllShape=x$shape)
            }

            outNames <- paste0(names(x), "_", 100*input$hitpLevel, "%")
            setNames(x, outNames)
        } else {
            NULL
        }
    })

    ## extrapolation to different distance output list
    extraListRadius <- reactive({
        xy <- coords()
        if(!is.null(xy)) {
            ## if no CEP type is selected -> fall back to default CorrNormal
            CEPtype <- if(!is.null(input$hitpCEPtype1)) {
                CEPtypesInv[input$hitpCEPtype1]
            } else {
                "CorrNormal"
            }

            ## hit probability -> radius
            res1 <- getCEP(xy,
                           CEPlevel=input$hitpLevel,
                           dstTarget=input$dstTrgt,
                           conversion=conversionStr(),
                           accuracy=input$hitpAcc,
                           type=CEPtype,
                           doRob=input$hitpDoRob1)$CEP[[1]]["MOA", , drop=FALSE]

            res2 <- getConfEll(xy,
                               level=input$hitpLevel,
                               dstTarget=input$dstTrgt,
                               conversion=conversionStr(),
                               doRob=input$hitpDoRob1)

            res3 <- if(input$hitpDoRob1) {
                res2$sizeRob["MOA", , drop=FALSE]
            } else {
                res2$size["MOA", , drop=FALSE]
            }

            CEP <- fromMOA(res1,
                           dst=input$hitpExtraDst1,
                           conversion=paste0(unitsDstInv[input$hitpUnitExtraDst1], "2",
                                             unitsXYInv[input$unitXY], collapse=""))
            ConfEll <- fromMOA(res3,
                               dst=input$hitpExtraDst1,
                               conversion=paste0(unitsDstInv[input$hitpUnitExtraDst1], "2",
                                                 unitsXYInv[input$unitXY], collapse=""))
            rownames(CEP)     <- "unit"
            rownames(ConfEll) <- "unit"
            x <- list(CEP=CEP, ConfEll=ConfEll)
            setNames(x, paste0(names(x), "_", 100*input$hitpLevel, "%", "_@",
                               input$hitpExtraDst1, unitsDstInv[input$hitpUnitExtraDst1]))
        } else {
            NULL
        }
    })

    ## CEP output
    output$CEPRadius <- renderPrint({
        CEPListRadius()
    })

    ## confidence ellipse output
    output$confEll <- renderPrint({
        confEllList()
    })

    ## extrapolation to different distance output
    output$extraRadius <- renderPrint({
        extraListRadius()
    })

    ## save output to text file
    output$saveRadius <- downloadHandler(
        filename=function() { "groupRadius.txt" },
        content=function(file) {
            writeLines(fileName(), con=file)
            out1 <- CEPListRadius()
            out2 <- confEllList()
            out3 <- extraListRadius()
            out  <- c(out1, out2, out3)
            Map(textOut, out, names(out), file)
        },
        contentType='text/plain' # MIME type
    )

    #####---------------------------------------------------------------------------
    ## radius -> hit probability
    #####---------------------------------------------------------------------------
    ## CEP output list - reactive conductor
    CEPListHitProb <- reactive({
        xy <- coords()
        if(!is.null(xy)) {
            ## if no CEP type is selected -> fall back to default CorrNormal
            CEPtype <- if(!is.null(input$hitpCEPtype2)) {
                CEPtypesInv[input$hitpCEPtype2]
            } else {
                "CorrNormal"
            }

            x <- getHitProb(xy,
                            r=input$hitpR,
                            unit=hitpRUnitInv[[input$hitpUnitR]],
                            dstTarget=input$dstTrgt,
                            conversion=conversionStr(),
                            accuracy=input$hitpAcc,
                            type=CEPtype,
                            doRob=input$hitpDoRob2)
            setNames(list(x), paste0("hitProbRadius_", input$hitpR))
        } else {
            NULL
        }
    })

    ## extrapolation to different distance output list
    extraListHitProb <- reactive({
        xy <- coords()
        if(!is.null(xy)) {
            ## if no CEP type is selected -> fall back to default CorrNormal
            CEPtype <- if(!is.null(input$hitpCEPtype2)) {
                CEPtypesInv[input$hitpCEPtype2]
            } else {
                "CorrNormal"
            }

            ## radius -> hit probability
            MOA <- getMOA(input$hitpR,
                          dst=input$hitpExtraDst2,
                          conversion=paste0(unitsDstInv[input$hitpUnitExtraDst2], "2",
                                            unitsXYInv[input$unitXY], collapse=""))
            x <- getHitProb(xy,
                            r=MOA,
                            unit="MOA",
                            dstTarget=input$dstTrgt,
                            conversion=conversionStr(),
                            accuracy=input$hitpAcc,
                            type=CEPtype,
                            doRob=input$hitpDoRob2)
            setNames(list(x), paste0("hitProbRadius_", input$hitpR, "_@",
                                     input$hitpExtraDst2, unitsDstInv[input$hitpUnitExtraDst2]))
        } else {
            NULL
        }
    })

    ## CEP output
    output$CEPHitProb <- renderPrint({
        CEPListHitProb()
    })

    ## extrapolation to different distance output
    output$extraHitProb <- renderPrint({
        extraListHitProb()
    })

    ## save output to text file
    output$saveHitProb <- downloadHandler(
        filename=function() { "groupHitProbability.txt" },
        content=function(file) {
            writeLines(fileName(), con=file)
            out1 <- CEPListHitProb()
            out2 <- extraListHitProb()
            out  <- c(out1, out3)
            Map(textOut, out, names(out), file)
        },
        contentType='text/plain' # MIME type
    )
})
