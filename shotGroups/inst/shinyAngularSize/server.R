source("helper.R")

shinyServer(function(input, output) {
    #####---------------------------------------------------------------------------
    ## distance to target, unit distance, unit xy-coords - UI element
    #####---------------------------------------------------------------------------
    output$unitDstXY <- renderUI({
        dstTarget <- 100
        ## angular size -> dst to target, and unit dst, but not unit xy
        inputPanel()
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
    ## absolute to angular size
    #####---------------------------------------------------------------------------

    output$abs2Ang <- renderPrint({
        dstTrgt1     <- input$dstTrgt1
        unitDstTrgt1 <- unitsDstInv[input$unitDst1]
        szeAbs1      <- as.numeric(strsplit(input$angszeAbs1, "[[:blank:]]")[[1]])
        unitAbs1     <- unitsAbsInv[input$angszeUnitAbs1]
        unitAngOut1  <- unitsAngInv[input$angszeUnitAngOut1]

        getMOAmult <- function(uao) {
            x <- getMOA(szeAbs1, dst=dstTrgt1,
                   conversion=paste0(unitDstTrgt1, "2", unitAbs1, collapse=""),
                   type=uao)
            signif(x, 4)
        }

        values <- Map(getMOAmult, unitAngOut1)
        x <- setNames(as.data.frame(values), unitAngOut1)
        if(!any(is.na(szeAbs1)) && !any(duplicated(szeAbs1))) {
            rownames(x) <- szeAbs1
        } else {
            NULL
        }

        x
    })

    #####---------------------------------------------------------------------------
    ## angular to absolute size
    #####---------------------------------------------------------------------------

    output$ang2Abs <- renderPrint({
        dstTrgt2     <- input$dstTrgt2
        unitDstTrgt2 <- unitsDstInv[input$unitDst2]
        szeAng2     <- as.numeric(strsplit(input$angszeAng2, "[[:blank:]]")[[1]])
        unitAng2    <- unitsAngInv[input$angszeUnitAng2]
        unitAbsOut2 <- unitsAbsInv[input$angszeUnitAbsOut2]

        fromMOAmult <- function(uao) {
            x <- fromMOA(szeAng2, dst=dstTrgt2,
                    conversion=paste0(unitDstTrgt2, "2", uao, collapse=""),
                    type=unitAng2)
            signif(x, 4)
        }

        values <- Map(fromMOAmult, unitAbsOut2)
        x <- setNames(as.data.frame(values), unitAbsOut2)
        if(!any(is.na(szeAng2)) && !any(duplicated(szeAng2))) {
            rownames(x) <- szeAng2
        } else {
                NULL
        }

        x
    })

    output$angAbs2Dist <- renderPrint({
        szeAbs3     <- as.numeric(strsplit(input$angszeAbs3, "[[:blank:]]")[[1]])
        szeAng3     <- as.numeric(strsplit(input$angszeAng3, "[[:blank:]]")[[1]])
        unitAbs3    <- unitsAbsInv[input$angszeUnitAbs3]
        unitAng3    <- unitsAngInv[input$angszeUnitAng3]
        unitDstOut3 <- unitsAbsInv[input$angszeUnitDstOut3]

        getDistMult <- function(udo) {
            x <- getDistance(szeAbs3,
                        angular=szeAng3,
                        conversion=paste0(udo, "2", unitAbs3),
                        type=unitAng3)
            signif(x, 4)
        }

        values <- Map(getDistMult, unitDstOut3)
        x <- setNames(as.data.frame(values), unitDstOut3)
        rNames <- paste0(szeAbs3, "_", szeAng3)
        if(!any(is.na(szeAbs3) | is.na(szeAng3)) && !any(duplicated(rNames))) {
            rownames(x) <- rNames
        } else {
            NULL
        }

        x
    })
})
