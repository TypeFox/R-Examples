library(DVHmetrics)
source("helper.R")

## max number of IDs (=DVH files) we can make plots for
maxNID <- 100

shinyServer(function(input, output) {
    ## directory where current DVH data is saved
    DVHdir <- tempdir()

    ## reactive conductor
    DVH <- reactive({
        input$applyData
        isolate({
            DVHdata <- if(input$DVHin == '1') {
                dataMZ
            } else if(input$DVHin == '2') {
                if(!is.null(input$DVHupload)) {
                    types <- c('1'="Eclipse", '2'="Cadplan", '3'="Masterplan",
                               '4'="Pinnacle", '5'="Monaco", '6'="HiArt",
                               '7'="RayStation")
                    plans <- c('1'="none", '2'="doseRx")
                    argL  <- if(("DVHadd" %in% input$DVHreadOpts) && file.exists("DVHprev.rds")) {
                        setwd(DVHdir)
                        add <- readRDS("DVHprev.Rds")
                        list(x=input$DVHupload$datapath,
                             type=types[input$DVHtype],
                             courseAsID=("DVHcourse" %in% input$DVHreadOpts),
                             planInfo=plans[input$DVHplanInfo],
                             add=add)
                    } else {
                        list(x=input$DVHupload$datapath,
                             type=types[input$DVHtype],
                             courseAsID=("DVHcourse" %in% input$DVHreadOpts),
                             planInfo=plans[input$DVHplanInfo])
                    }
                    do.call(readDVH, argL)
                } else {
                    NULL
                }
            }

            ## add representation by structure
            DVHdataByStruct <- if(!is.null(DVHdata)) {
                DVHmetrics:::reorgByPat(DVHdata, byPat=FALSE)
            } else {
                NULL
            }

            DvHcurr <- list(DVH=DVHdata, DVHbyStruct=DVHdataByStruct)
            setwd(DVHdir)
            saveRDS(DvHcurr$DVH, file="DVHprev.rds")
            return(DvHcurr)
        })
    })

    output$DVHinfo <- renderPrint({
        dvh <- DVH()$DVH
        if(!is.null(dvh)) {
            if(input$DVHverbose == '1') {
                print(dvh, verbose=FALSE)
            } else {
                print(dvh, verbose=TRUE)
            }
        } else {
            NULL
        }
    })

    output$metrSelStruct <- renderUI({
        dvh <- DVH()$DVH
        if(!is.null(dvh)) {
            checkboxGroupInput("metrSelStruct",
                               label=h5("Select structures"),
                               choices=getStrIDs(dvh, what="structure", choices=TRUE),
                               selected=1)
        } else {
            NULL
        }
    })

    output$metrSelPat <- renderUI({
        dvh <- DVH()$DVH
        if(!is.null(dvh)) {
            checkboxGroupInput("metrSelPat",
                               label=h5("Select patients"),
                               choices=getStrIDs(dvh, what="patient", choices=TRUE),
                               selected=1)
        } else {
            NULL
        }
    })

    output$plotSelStruct <- renderUI({
        dvh <- DVH()$DVH
        if(!is.null(dvh)) {
            checkboxGroupInput("plotSelStruct",
                               label=h5("Select structures"),
                               choices=getStrIDs(dvh, what="structure", choices=TRUE),
                               selected=1)
        } else {
            NULL
        }
    })

    output$plotSelPat <- renderUI({
        dvh <- DVH()$DVH
        if(!is.null(dvh)) {
            checkboxGroupInput("plotSelPat",
                               label=h5("Select patients"),
                               choices=getStrIDs(dvh, what="patient", choices=TRUE),
                               selected=1)
        } else {
            NULL
        }
    })

    output$metrics <- renderDataTable({
        dvh        <- DVH()$DVH
        sortOpts   <- c('1'="observed", '2'="structure", '3'="metric", '4'="patID")
        splitOpts  <- c('1'="structure", '2'="metric", '3'="patID")
        selMetrics <- if(length(input$metrInput) > 0) {
            metrRaw <- unlist(strsplit(input$metrInput, "[[:blank:],]"))
            metrRaw[nzchar(metrRaw)]
        } else {
            c("DMEAN", "D1CC", "V10%")
        }
        selStruct <- if(length(input$metrSelStruct) > 0) {
            if(!is.null(dvh)) {
                getStrIDs(dvh, what="structure")[as.numeric(input$metrSelStruct)]
            } else {
                NULL
            }
        } else {
            NULL
        }
        selPat <- if(length(input$metrSelPat) > 0) {
            if(!is.null(dvh)) {
                getStrIDs(dvh, what="patient")[as.numeric(input$metrSelPat)]
            } else {
                NULL
            }
        } else {
            NULL
        }
        interp  <- "linear" # c("linear", "spline", "ksmooth")[as.numeric(input$metrInterp)]
        EUDa    <- if(input$metrEUDa  != "") { as.numeric(input$metrEUDa)  } else { NULL }
        EUDfd   <- if(input$metrEUDfd != "") { as.numeric(input$metrEUDfd) } else { NULL }
        EUDab   <- if(input$metrEUDab != "") { as.numeric(input$metrEUDab) } else { NULL }

        NTCPtype    <- c("probit", "logit", "poisson")[as.numeric(input$metrNTCPtype)]
        NTCPtd50    <- if(input$metrNTCPtd50    != "") { as.numeric(input$metrNTCPtd50)    } else { NULL }
        NTCPn       <- if(input$metrNTCPn       != "") { as.numeric(input$metrNTCPn)       } else { NULL }
        NTCPm       <- if(input$metrNTCPm       != "") { as.numeric(input$metrNTCPm)       } else { NULL }
        NTCPgamma50 <- if(input$metrNTCPgamma50 != "") { as.numeric(input$metrNTCPgamma50) } else { NULL }

        sortSel <- input$metrSortBy
        sortBy <- if(length(sortSel) > 0) {
            sortOpts[sortSel]
        } else {
            NULL
        }
        if(!is.null(dvh)) {
            argL <- list(x=dvh,
                         metric=selMetrics,
                         patID=selPat,
                         structure=selStruct,
                         sortBy=sortBy,
                         interp=interp,
                         EUDa=EUDa, EUDfd=EUDfd, EUDab=EUDab,
                         NTCPtype=NTCPtype, NTCPtd50=NTCPtd50, NTCPn=NTCPn, NTCPm=NTCPm, NTCPgamma50=NTCPgamma50,
                          TCPtype=NTCPtype, TCPtcd50=NTCPtd50,  TCPn=NTCPn,  TCPm=NTCPm,  TCPgamma50=NTCPgamma50)
            argL <- Filter(Negate(is.null), argL)
            metr <- do.call(getMetric, argL)
            metr$observed <- round(metr$observed, 2)
            metr
        } else {
            NULL
        }
    })#, options=list(pageLength=25))

    output$saveMetrics <- downloadHandler(
        filename=function() { "metrics.txt" },
        content=function(file) {
            dvh       <- DVH()$DVH
            sortOpts  <- c('1'="observed", '2'="structure", '3'="metric", '4'="patID")
            splitOpts <- c('1'="metric", '2'="structure", '3'="patID")
            selMetrics <- if(length(input$metrInput) > 0) {
                metrRaw <- unlist(strsplit(input$metrInput, "[[:blank:],]"))
                metrRaw[nzchar(metrRaw)]
            } else {
                "DMEAN"
            }
            selPat <- if(length(input$metrSelPat) > 0) {
                if(!is.null(dvh)) {
                    getStrIDs(dvh, what="patient")[as.numeric(input$metrSelPat)]
                } else {
                    NULL
                }
            } else {
                NULL
            }
            selStruct <- if(length(input$metrSelStruct) > 0) {
                if(!is.null(dvh)) {
                    getStrIDs(dvh, what="structure")[as.numeric(input$metrSelStruct)]
                } else {
                    NULL
                }
            } else {
                NULL
            }

            interp  <- "linear"
            EUDa    <- if(input$metrEUDa  != "") { as.numeric(input$metrEUDa)  } else { NULL }
            EUDfd   <- if(input$metrEUDfd != "") { as.numeric(input$metrEUDfd) } else { NULL }
            EUDab   <- if(input$metrEUDab != "") { as.numeric(input$metrEUDab) } else { NULL }

            NTCPtype    <- c("probit", "logit", "poisson")[as.numeric(input$metrNTCPtype)]
            NTCPtd50    <- if(input$metrNTCPtd50    != "") { as.numeric(input$metrNTCPtd50)    } else { NULL }
            NTCPn       <- if(input$metrNTCPn       != "") { as.numeric(input$metrNTCPn)       } else { NULL }
            NTCPm       <- if(input$metrNTCPm       != "") { as.numeric(input$metrNTCPm)       } else { NULL }
            NTCPgamma50 <- if(input$metrNTCPgamma50 != "") { as.numeric(input$metrNTCPgamma50) } else { NULL }

            sortSel <- input$metrSortBy
            sortBy <- if(length(sortSel) > 0) {
                sortOpts[sortSel]
            } else {
                NULL
            }
            argL <- list(x=dvh,
                         metric=selMetrics,
                         patID=selPat,
                         structure=selStruct,
                         sortBy=sortBy,
                         interp=interp,
                         EUDa=EUDa, EUDfd=EUDfd, EUDab=EUDab,
                         NTCPtype=NTCPtype, NTCPtd50=NTCPtd50, NTCPn=NTCPn, NTCPgamma50=NTCPgamma50,
                          TCPtype=NTCPtype, TCPtcd50=NTCPtd50,  TCPn=NTCPn,  TCPgamma50=NTCPgamma50)
            argL <- Filter(Negate(is.null), argL)
            metr <- do.call(getMetric, argL)
            dec <- c('1'=".",  '2'=",")[input$saveMetrDec]
            sep <- c('1'="\t", '2'=" ", '3'=",", '4'=";")[input$saveMetrSep]
            write.table(metr, file=file, dec=dec, sep=sep, row.names=FALSE)
        },
        contentType='text/plain' # MIME type
    )

#     output$DVHplotOrg <- renderPlot({
#         dvh <- DVH()$DVH
#         if(!is.null(dvh)) {
#             argL <- list(x=dvh,
#                          patID=    getStrIDs(dvh, what="patient")[  as.numeric(input$plotSelPat)],
#                          structure=getStrIDs(dvh, what="structure")[as.numeric(input$plotSelStruct)],
#                          thresh=input$plotThreshVol,
#                          rel=input$plotPlotVol == '1')
#             do.call(showDVH, argL)
#         } else {
#             NULL
#         }
#     })

    output$DVHplot <- renderUI({
        dvh       <- DVH()
        byPat     <- input$plotByPat == '1'
        selPat    <- getStrIDs(dvh$DVH, what="patient")[  as.numeric(input$plotSelPat)]
        selStruct <- getStrIDs(dvh$DVH, what="structure")[as.numeric(input$plotSelStruct)]
        plotOutputL <- if(byPat) {
            ## 1 diagram per patient
            lapply(seq_along(dvh$DVH), function(i) {
                iMax <- sum(selPat %in% names(dvh$DVH))
                if(i <= iMax) {    # avoid creating empty plots
                    plotOutput(paste0("DVHplot", i))
                } else {
                    NULL
                }
            })
        } else {
            ## 1 diagram per structure
            lapply(seq_along(dvh$DVHbyStruct), function(i) {
                iMax <- sum(selStruct %in% names(dvh$DVHbyStruct))
                if(i <= iMax) {    # avoid creating empty plots
                    plotOutput(paste0("DVHplot", i))
                } else {
                    NULL
                }
            })
        }

        ## weed out NULL components, convert the list to a tagList and return
        plotOutputL <- Filter(Negate(is.null), plotOutputL)
        do.call(tagList, plotOutputL)
    })

    ## adapted from Winston Chang: https://gist.github.com/wch/5436415
    ## call renderPlot maxNID times and discard those not required (return NULL)
    ## problem: number of plots should be dynamic -> reactive context required
    for(i in seq_len(maxNID)) {
    ## need local so that each item gets its own number. Without it, the value
    ## of i in the renderPlot() will be the same across all instances, because
    ## of when the expression is evaluated.
        local({
            localI <- i
            ## DVH plot
            output[[paste0("DVHplot", localI)]] <- renderPlot({
                dvh   <- DVH()
                rel   <- input$plotPlotVol == '1'
                cumul <- input$plotType    == '1'
                byPat <- input$plotByPat   == '1'
                selPat    <- getStrIDs(dvh$DVH, what="patient")[  as.numeric(input$plotSelPat)]
                selStruct <- getStrIDs(dvh$DVH, what="structure")[as.numeric(input$plotSelStruct)]
                if(( byPat && (sum(selPat    %in% names(dvh$DVH))         < localI)) ||
                   (!byPat && (sum(selStruct %in% names(dvh$DVHbyStruct)) < localI))) {
                        NULL
                } else {
                    ## 1 diagram per patient/structure
                    x <- if(byPat) {
                        sharedNames <- intersect(selPat, names(dvh$DVH))
                        dvh$DVH[sharedNames][[localI]]
                    } else {
                        sharedNames <- intersect(selStruct, names(dvh$DVHbyStruct))
                        dvh$DVHbyStruct[sharedNames][[localI]]
                    }

                    showDVH(x=x,
                            cumul=cumul,
                            byPat=byPat,
                            patID=selPat,
                            structure=selStruct,
                            thresh=input$plotThreshVol,
                            rel=rel,
                            addMSD=input$plotMSD)
                }
            })

            ## constraint plot
            output[[paste0("constraintPlot", localI)]] <- renderPlot({
                dvh    <- DVH()
                constr <- DVHconstr()
                rel    <- input$constrPlotVol == "1"
                byPat  <- input$constrByPat   == "1"
                ## 1 diagram per patient/structure
                ## restrict DVH and constr to the same IDs/structures
                xConstrSub <- if(byPat) {
                    DVHmetrics:::harmoConstrDVH.DVHLstLst(dvh$DVH,         constr=constr, byPat=byPat)
                } else {
                    DVHmetrics:::harmoConstrDVH.DVHLstLst(dvh$DVHbyStruct, constr=constr, byPat=byPat)
                }

                if(is.null(constr) || (length(xConstrSub$x) < localI)) {
                    NULL
                } else {
                    showConstraint(x=xConstrSub$x[[localI]],
                                   constr=xConstrSub$constr[[localI]],
                                   byPat=byPat,
                                   thresh=input$constrThreshVol,
                                   rel=rel)
                }
            })
        })
    }

    output$saveDVHPDF <- downloadHandler(
        filename=function() { "DVH.pdf" },
        content=function(file) {
            dvh   <- DVH()
            byPat <- input$plotByPat   == '1'
            rel   <- input$plotPlotVol == "1"
            cumul <- input$plotType    == "1"
            selPat    <- getStrIDs(dvh$DVH, what="patient")[  as.numeric(input$plotSelPat)]
            selStruct <- getStrIDs(dvh$DVH, what="structure")[as.numeric(input$plotSelStruct)]
            x <- if(byPat) {
                dvh$DVH
            } else {
                dvh$DVHbyStruct
            }

            argL <- list(x=x,
                         cumul=cumul,
                         byPat=byPat,
                         patID=selPat,
                         structure=selStruct,
                         thresh=input$plotThreshVol,
                         rel=rel,
                         addMSD=input$plotMSD)
            pdf(file, width=7, height=5)
            do.call(showDVH, argL)
            dev.off()
        },
        contentType='application/pdf' # MIME type
    )

    makeDVHJPG <- function(fName, dvh, byPat, rel, cumul, selPat, selStruct, thresh) {
        argL <- list(x=dvh,
                     cumul=cumul,
                     byPat=byPat,
                     patID=selPat,
                     structure=selStruct,
                     thresh=thresh,
                     rel=rel)
        jpeg(fName, quality=100, width=700, height=500)
        do.call(showDVH, argL)
        dev.off()
    }

    output$saveDVHJPG <- downloadHandler(
        filename=function() { "DVH_JPGs.zip" },
        content=function(file) {
            dvh    <- DVH()
            byPat  <- input$plotByPat   == '1'
            rel    <- input$plotPlotVol == "1"
            cumul  <- input$plotType    == "1"
            thresh <- input$plotThreshVol
            addMSD <- input$plotMSD
            x <- if(byPat) {
                dvh$DVH
            } else {
                dvh$DVHbyStruct
            }

            selPat    <- getStrIDs(dvh$DVH, what="patient")[  as.numeric(input$plotSelPat)]
            selStruct <- getStrIDs(dvh$DVH, what="structure")[as.numeric(input$plotSelStruct)]
            nFiles <- if(byPat) {
                sum(selPat    %in% names(x))
            } else {
                sum(selStruct %in% names(x))
            }

            ## write all JPEGs into temporary directory and zip them
            fNames <- paste0("DVH_", sprintf("%02d", seq_len(nFiles)), ".jpg")
            tmpdir <- tempdir()
            setwd(tmpdir)

            for(i in seq_along(fNames)) {
                fN         <- fNames[i]
                selPatI    <- if(byPat) { selPat[i] } else { selPat }
                selStructI <- if(byPat) { selStruct } else { selStruct[i] }
                makeDVHJPG(fN, dvh=x, byPat=byPat, rel=rel, cumul=cumul,
                           selPat=selPatI, selStruct=selStructI,
                           addMSD=addMSD, thresh=thresh)
            }

            zip(zipfile=file, files=fNames)
            if(file.exists(paste0(file, ".zip"))) {
                file.rename(paste0(file, ".zip"), file)
            }
        },
        contentType = "application/zip"
    )

    output$saveDVHMSD <- downloadHandler(
        filename=function() { "DVH_M_SD.txt" },
        content=function(file) {
            dvh    <- DVH()
            byPat  <- input$plotByPat   == '1'
            cumul  <- input$plotType    == '1'
            x <- if(byPat) {
                dvh$DVH
            } else {
                dvh$DVHbyStruct
            }

            selPat    <- getStrIDs(dvh$DVH, what="patient")[  as.numeric(input$plotSelPat)]
            selStruct <- getStrIDs(dvh$DVH, what="structure")[as.numeric(input$plotSelStruct)]

            argL <- list(x=x,
                         fun=list(mean=mean, median=median, min=min, max=max, sd=sd),
                         byPat=byPat,
                         thin=1,
                         cumul=cumul,
                         patID=selPat,
                         structure=selStruct,
                         interp="linear",
                         fixed=TRUE)
            argL   <- Filter(Negate(is.null), argL)
            DVHMSD <- do.call(getMeanDVH, argL)
            dec    <- c('1'=".",  '2'=",")[input$saveDVHDec]
            sep    <- c('1'="\t", '2'=" ", '3'=",", '4'=";")[input$saveDVHSep]
            write.table(DVHMSD, file=file, dec=dec, sep=sep, row.names=FALSE)
        },
        contentType='text/plain' # MIME type
    )

    ## reactive conductor
    DVHconstr <- reactive({
        input$applyConstraints
        isolate({
            constr <- if(input$constrIn == '1') {
                dataConstr
            } else if(input$constrIn == '2') {
                if(is.null(input$constrUpload)) {
                    NULL
                } else {
                    dec <- c('1'=".",  '2'=",")[input$constrDec]
                    sep <- c('1'="\t", '2'=" ", '3'=",", '4'=";")[input$constrSep]
                    ## peek at data to see if we have a header
                    line1 <- readLines(input$constrUpload$datapath, n=1)
                    if(grepl("constraint", line1, ignore.case=TRUE)) {
                        read.table(file=input$constrUpload$datapath,
                                   header=TRUE, sep=sep, dec=dec, stringsAsFactors=FALSE)
                    } else {
                        read.table(file=input$constrUpload$datapath,
                                   header=FALSE, sep=sep, dec=dec, stringsAsFactors=FALSE)
                    }
                }
            } else if(input$constrIn == '3') {
                dec <- c('1'=".",  '2'=",")[input$constrPasteDec]
                sep <- c('1'="\t", '2'=" ", '3'=",", '4'=";")[input$constrPasteSep]
                ## peek at data to see if we have a header
                line1 <- readLines(textConnection(input$constrPaste), n=1)
                if(grepl("constraint", line1, ignore.case=TRUE)) {
                    read.table(file=textConnection(input$constrPaste),
                               header=TRUE, sep=sep, dec=dec, stringsAsFactors=FALSE)
                } else {
                    read.table(file=textConnection(input$constrPaste),
                               header=FALSE, sep=sep, dec=dec, stringsAsFactors=FALSE)
                }
            } else {
                NULL
            }

            if(!is.null(constr)) {
                cNames <- tolower(names(constr))
                cNames[cNames == "patid"] <- "patID"
                constr <- setNames(constr, cNames)
            }

            return(constr)
        })
    })

    output$constraints <- renderDataTable({
        constr <- DVHconstr()
        dvh    <- DVH()$DVH
        outSel <- constrOutInv[input$constrOut]
        interp <- "linear" #c("linear", "spline", "ksmooth")[as.numeric(input$constrInterp)]
        EUDa   <- if(input$constrEUDa  != "") { as.numeric(input$constrEUDa)  } else { NULL }
        EUDfd  <- if(input$constrEUDfd != "") { as.numeric(input$constrEUDfd) } else { NULL }
        EUDab  <- if(input$constrEUDab != "") { as.numeric(input$constrEUDab) } else { NULL }

        NTCPtype    <- c("probit", "logit", "poisson")[as.numeric(input$constrNTCPtype)]
        NTCPtd50    <- if(input$constrNTCPtd50    != "") { as.numeric(input$constrNTCPtd50)    } else { NULL }
        NTCPn       <- if(input$constrNTCPn       != "") { as.numeric(input$constrNTCPn)       } else { NULL }
        NTCPm       <- if(input$constrNTCPm       != "") { as.numeric(input$constrNTCPm)       } else { NULL }
        NTCPgamma50 <- if(input$constrNTCPgamma50 != "") { as.numeric(input$constrNTCPgamma50) } else { NULL }

        sortOpts <- c('1'="compliance", '2'="dstMin", '3'="deltaV", '4'="deltaD",
                      '5'="observed", '6'="patID", '7'="structure", '8'="constraint")
        sortSel  <- input$constrSortBy
        sortBy   <- if(length(sortSel) > 0) {
            sortOpts[sortSel]
        } else {
            "none"
        }

        if(!is.null(constr) && !is.null(dvh)) {
            argL <- list(x=dvh,
                         constr=constr, byPat=TRUE, interp=interp,
                         semSign=input$constrSemSign, sortBy=sortBy,
                         EUDa=EUDa, EUDfd=EUDfd, EUDab=EUDab,
                         NTCPtype=NTCPtype, NTCPtd50=NTCPtd50, NTCPn=NTCPn, NTCPm=NTCPm, NTCPgamma50=NTCPgamma50,
                          TCPtype=NTCPtype, TCPtcd50=NTCPtd50,  TCPn=NTCPn,  TCPm=NTCPm,  TCPgamma50=NTCPgamma50)
            argL <- Filter(Negate(is.null), argL)
            x <- do.call(checkConstraint, argL)
            x$observed  <- round(x$observed,  2)
            x$deltaVpc  <- round(x$deltaVpc,  2)
            x$deltaDpc  <- round(x$deltaDpc,  2)
            x$dstMin    <- round(x$dstMin,    2)
            x$dstMinRel <- round(x$dstMinRel, 2)
            x[ , outSel]
        } else {
            NULL
        }

    })#, options=list(pageLength=25))

#     output$constraintPlotOrg <- renderPlot({
#         constr <- DVHconstr()
#         dvh    <- DVH()$DVH
#         if(!is.null(constr) && !is.null(dvh)) {
#             argL <- list(x=dvh,
#                          thresh=input$constrThreshVol,
#                          rel=input$constrPlotVol == '1',
#                          constr=constr)
#             do.call(showConstraint, argL)
#         } else {
#             NULL
#         }
#     })

    output$constraintPlot <- renderUI({
        dvh    <- DVH()
        constr <- DVHconstr()
        byPat  <- input$constrByPat == "1"
        xConstrSub <- if(byPat) {
            DVHmetrics:::harmoConstrDVH.DVHLstLst(dvh$DVH,         constr=constr, byPat=byPat)
        } else {
            DVHmetrics:::harmoConstrDVH.DVHLstLst(dvh$DVHbyStruct, constr=constr, byPat=byPat)
        }

        plotOutputL <- lapply(seq_along(xConstrSub$x), function(i) {
            plotOutput(paste0("constraintPlot", i)) })

        ## convert the list to a tagList and return
        do.call(tagList, plotOutputL)
    })

    output$saveConstrTxt <- downloadHandler(
        filename=function() { "constraints.txt" },
        content=function(file) {
            interp <-"linear"
            EUDa   <- if(input$constrEUDa  != "") { as.numeric(input$constrEUDa)  } else { NULL }
            EUDfd  <- if(input$constrEUDfd != "") { as.numeric(input$constrEUDfd) } else { NULL }
            EUDab  <- if(input$constrEUDab != "") { as.numeric(input$constrEUDab) } else { NULL }

            NTCPtype    <- c("probit", "logit", "poisson")[as.numeric(input$constrNTCPtype)]
            NTCPtd50    <- if(input$constrNTCPtd50    != "") { as.numeric(input$constrNTCPtd50)    } else { NULL }
            NTCPn       <- if(input$constrNTCPn       != "") { as.numeric(input$constrNTCPn)       } else { NULL }
            NTCPm       <- if(input$constrNTCPm       != "") { as.numeric(input$constrNTCPm)       } else { NULL }
            NTCPgamma50 <- if(input$constrNTCPgamma50 != "") { as.numeric(input$constrNTCPgamma50) } else { NULL }

            argL <- list(x=DVH()$DVH,
                         constr=DVHconstr(), interp=interp,
                         EUDa=EUDa, EUDfd=EUDfd, EUDab=EUDab,
                         NTCPtype=NTCPtype, NTCPtd50=NTCPtd50, NTCPn=NTCPn, NTCPm=NTCPm, NTCPgamma50=NTCPgamma50,
                          TCPtype=NTCPtype, TCPtcd50=NTCPtd50,  TCPn=NTCPn,  TCPm=NTCPm,  TCPgamma50=NTCPgamma50)
            argL   <- Filter(Negate(is.null), argL)
            constr <- do.call(checkConstraint, argL)
            dec    <- c('1'=".",  '2'=",")[input$saveConstrDec]
            sep    <- c('1'="\t", '2'=" ", '3'=",", '4'=";")[input$saveConstrSep]
            write.table(constr, file=file, dec=dec, sep=sep, row.names=FALSE)
        },
        contentType='text/plain' # MIME type
    )

    output$saveConstrPDF <- downloadHandler(
        filename=function() { "constraints.pdf" },
        content=function(file) {
            rel   <- input$constrPlotVol == '1'
            byPat <- input$constrByPat   == '1'
            pdf(file, width=7, height=5)
            showConstraint(x=DVH()$DVH,
                           constr=DVHconstr(),
                           byPat=byPat,
                           thresh=input$constrThreshVol,
                           rel=rel)
            dev.off()
        },
        contentType='application/pdf' # MIME type
    )

    makeConstrJPG <- function(fName, dvh, constr, byPat, thresh, rel) {
        argL <- list(x=dvh,
                     constr=constr,
                     byPat=byPat,
                     thresh=thresh,
                     rel=rel)
        jpeg(fName, quality=100, width=700, height=500)
        do.call(showConstraint, argL)
        dev.off()
    }

    output$saveConstrJPG <- downloadHandler(
        filename=function() { "constraints.zip" },
        content=function(file) {
            dvh    <- DVH()
            constr <- DVHconstr()
            rel    <- input$constrPlotVol == '1'
            byPat  <- input$constrByPat   == '1'
            thresh <- input$constrThreshVol
            x <- if(byPat) {
                dvh$DVH
            } else {
                dvh$DVHbyStruct
            }

            xConstrSub <- DVHmetrics:::harmoConstrDVH.DVHLstLst(x, constr=constr, byPat=byPat)
            nFiles     <- length(xConstrSub$x)

            ## write all JPEGs into temporary directory and zip them
            fNames <- paste0("DVH_", sprintf("%02d", seq_len(nFiles)), ".jpg")
            tmpdir <- tempdir()
            setwd(tmpdir)

            for(i in seq_along(fNames)) {
                fN <- fNames[i]
                makeConstrJPG(fN, dvh=xConstrSub$x, constr=xConstrSub$constr,
                              byPat=byPat, rel=rel, thresh=thresh)
            }

            zip(zipfile=file, files=fNames)
            if(file.exists(paste0(file, ".zip"))) {
                file.rename(paste0(file, ".zip"), file)
            }
        },
        contentType = "application/zip"
    )

    output$BED <- renderPrint({
        if(input$BEDtype %in% c('1', '2')) {
            D <- if(input$BED_BED_D  != "") {
                as.numeric(strsplit(input$BED_BED_D, "[[:blank:]]")[[1]])
            } else {
                NULL
            }
            fd <- if(input$BED_BED_FD != "") {
                as.numeric(strsplit(input$BED_BED_FD, "[[:blank:]]")[[1]])
            } else {
                NULL
            }
            fn <- if(input$BED_BED_FN != "") {
                as.numeric(strsplit(input$BED_BED_FN, "[[:blank:]]")[[1]])
            } else {
                NULL
            }
            ab <- if(input$BED_BED_AB != "") {
                as.numeric(strsplit(input$BED_BED_AB, "[[:blank:]]")[[1]])
            } else {
                NULL
            }
        } else if(input$BEDtype == '3') {
            D1 <- if(input$BED_IED_D1  != "") {
                as.numeric(strsplit(input$BED_IED_D1, "[[:blank:]]")[[1]])
            } else {
                NULL
            }
            D2 <- if(input$BED_IED_D2  != "") {
                as.numeric(strsplit(input$BED_IED_D2, "[[:blank:]]")[[1]])
            } else {
                NULL
            }
            fd1 <- if(input$BED_IED_FD1 != "") {
                as.numeric(strsplit(input$BED_IED_FD1, "[[:blank:]]")[[1]])
            } else {
                NULL
            }
            fd2 <- if(input$BED_IED_FD2 != "") {
                as.numeric(strsplit(input$BED_IED_FD2, "[[:blank:]]")[[1]])
            } else {
                NULL
            }
            fn1 <- if(input$BED_IED_FN1 != "") {
                as.numeric(strsplit(input$BED_IED_FN1, "[[:blank:]]")[[1]])
            } else {
                NULL
            }
            fn2 <- if(input$BED_IED_FN2 != "") {
                as.numeric(strsplit(input$BED_IED_FN2, "[[:blank:]]")[[1]])
            } else {
                NULL
            }
            ab <- if(input$BED_IED_AB != "") {
                as.numeric(strsplit(input$BED_IED_AB, "[[:blank:]]")[[1]])
            } else {
                NULL
            }
        }

        if(input$BEDtype == '1') {
            getBED(D=D, fd=fd, fn=fn, ab=ab)
        } else if(input$BEDtype == '2') {
            getEQD2(D=D, fd=fd, fn=fn, ab=ab)
        } else if(input$BEDtype == '3') {
            getIsoEffD(D1=D1, D2=D2, fd1=fd1, fd2=fd2, ab=ab)
        }
    })
})
