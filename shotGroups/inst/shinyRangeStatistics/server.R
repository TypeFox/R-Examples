source("helper.R")

shinyServer(function(input, output) {
    #####---------------------------------------------------------------------------
    ## Range statistics -> Rayleigh sigma
    #####---------------------------------------------------------------------------
    output$range2sigma <- renderPrint({
        rangeStatStats <- as.numeric(strsplit(input$rangeStatStats, "[[:blank:]]")[[1]])
        stat           <- rangeStatInv[input$rangeStatType1]
        n              <- input$rangeStatN
        nGroups        <- input$rangeStatNGroups
        CIlevel        <- input$rangeStatCILevel
        dstTrgt1       <- input$dstTrgt1
        unitDstTrgt1   <- unitsDstInv[input$unitDstTrgt1]
        unitXY1        <- unitsXYInv[input$unitXY1]
        conversion     <- paste0(unitDstTrgt1, "2", unitXY1, collapse="")
        range2sigma(rangeStatStats, stat=stat, n=n, nGroups=nGroups,
                    CIlevel=CIlevel, collapse=TRUE, dstTarget=dstTrgt1,
                    conversion=conversion)
    })

    ## save output to file
    output$saveShape <- downloadHandler(
        filename=function() { "range2sigma.txt" },
        content=function(file) {
            writeLines(fileName(), con=file)
            "ABC"
        },
        contentType='text/plain' # MIME type
    )

    #####---------------------------------------------------------------------------
    ## Efficiency - required number of groups
    #####---------------------------------------------------------------------------
    output$effNGroups <- renderPrint({
        stat    <- rangeStatSigInv[input$effStatType1]
        n       <- as.numeric(strsplit(input$effN1, "[[:blank:]]")[[1]])
        CIlevel <- input$effCILevel1
        CIwidth <- input$effCIWidth1
        efficiency(n=n, CIlevel=CIlevel, CIwidth=CIwidth, stat=stat)
    })
# var tips = ['shots per group',
#             'required number of groups with n shots each (including fractional groups)',
#             'required number of groups with n shots each (full groups only)',
#             'required total number of shots, assuming we're shooting groups of size n each (including fractional groups)',
#             'required total number of shots, assuming we're shooting groups of size n each (full groups only)',
#             'desired CI level (coverage probability)',
#             'desired CI width (as fraction of the mean)'],
#     header = table.columns().header();
#     for (var i = 0; i < tips.length; i++) {
#         $(header[i]).attr('title', tips[i]);
#     }
# "))

    #####---------------------------------------------------------------------------
    ## Efficiency - achievable CI width
    #####---------------------------------------------------------------------------
    output$effCIWidth <- renderPrint({
        stat    <- rangeStatSigInv[input$effStatType2]
        n       <- as.numeric(strsplit(input$effN2, "[[:blank:]]")[[1]])
        nGroups <- input$effNGroups2
        CIlevel <- input$effCILevel2
        out <- efficiency(n=n, nGroups=nGroups, CIlevel=CIlevel, stat=stat)
        out <- cbind(out, E=out[ , "CIwidth"]/2)
        rownames(out) <- NULL
        out
    })
})
