pla.plots <-
function (dataframe,
          sampleLabels = levels(unlist(dataframe["Sample"])),
          indexOfReference = 1,
          design = "blocks",
          main = "Parallel Line Model",
          setPdf = FALSE,
          pdfName = "PlaPlots.pdf",
          nc = 4,
          mfrow = c(2, nc),
          oma = c(1, 0, 2, 0),
          mar = c(3.5, 3.5, 0, 0) + 0.6,
          joinReplicates = TRUE,
          showRho = FALSE,
          plots = "default",
          xlab = "Log(Dosis)",
          ylab = "Response: Response",
          pch = 14 + as.numeric(unlist(dataframe["Replicate"])),
          cex = 2,
          lwd = 4,
          colTst = "black",
          colRef = "grey",
          colRho = "grey10",
          colNrm = "grey70",
          tests = NULL) 
{
    selectS <- !is.na(match(as.character(unlist(dataframe["Sample"])),
                            sampleLabels))
    dataframe <- dataframe[selectS, ]
    design <- .string2design(design)
    if (any(25 < pch & pch < 32))
        pch <- pch %% 25
    latin <- FALSE
    blocks <- FALSE
    crd <- FALSE
    if ((design == "lsd"))
        latin <- TRUE
    if ((design == "rbd"))
        blocks <- TRUE
    if (design == "crd") 
        crd <- TRUE
    OK <- .checkPlaFrame(dataframe, design = design,
                         latin = latin, blocks = blocks, crd = crd)

    if (any(plots == "default")) {
        restricted             <- TRUE
        unrestricted           <- TRUE
        twoway                 <- (nc == 4)
        stepadjusted           <- TRUE
        qplot                  <- TRUE
        histogram              <- (nc == 4)
        residualsBoxSample     <- FALSE
        residualsBoxSampleStep <- FALSE
        residualsSample        <- FALSE
        residualsRow           <- FALSE
        residualsColumn        <- FALSE
        residualsBlock         <- FALSE
        residualsPlate         <- TRUE
        residualsStep          <- TRUE
    } else {
        restricted             <- any(plots == "all") |
            any(plots == "restricted")
        unrestricted           <- any(plots == "all") |
            any(plots == "unrestricted")
        twoway                 <- any(plots == "all") |
            any(plots == "twoway")
        stepadjusted           <- any(plots == "all") | 
            any(plots == "stepadjusted")

        qplot                  <- any(plots == "all") | 
            any(plots == "qplot")
        histogram              <- any(plots == "all") | 
            any(plots == "histogram")
        residualsBoxSample     <- any(plots == "ALL") | 
            any(plots == "residualsBoxSample")
        residualsBoxSampleStep <- any(plots == "all") | 
            any(plots == "residualsBoxSampleStep")

        residualsSample        <- any(plots == "all") | 
            any(plots == "residualsSample")
        residualsRow           <- any(plots == "all") | 
            any(plots == "residualsRow")
        residualsColumn        <- any(plots == "all") | 
            any(plots == "residualsColumn")
        residualsBlock         <- any(plots == "all") | 
            any(plots == "residualsBlock")

        residualsPlate         <- any(plots == "all") | 
            any(plots == "residualsPlate")
        residualsStep          <- any(plots == "all") | 
            any(plots == "residualsStep")
    }

    if (OK) {
        if (setPdf) 
            pdf(pdfName)
        if (!is.null(mfrow)) par(mfrow = mfrow)
        if (!is.null(oma))   par(oma = oma)
        if (!is.null(mar))   par(mar = mar)
        dataJitter <- jitterSteps(dataframe)
        nSamples <- length(levels(dataframe[,"Sample"]))-1
        if (nSamples+1 != length(sampleLabels))
            warning(
                paste0("Length of argument 'treamentLabels' does not ",
                       "correspond to number of treatements in data"))
        if (nSamples > length(colTst)) colTst <- rep(colTst, nSamples)

        LMrsimple <- lm(Response ~ factor(Sample) + Z, data = dataframe)
        if (crd) 
            LMr <- LMrsimple
        if (latin) 
            LMr <- lm(Response ~ factor(Row) + factor(Column) + factor(Sample) + 
                      Z, data = dataframe)
        if (blocks) 
            LMr <- lm(Response ~ factor(Replicate) + factor(Sample) + 
                      Z, data = dataframe)

        ## Restricted:
        if (restricted) {
            color <- ifelse(is.null(tests), "black",
                            ifelse(tests$Regression == "Test failed!",
                                   "red", "black"))
            par(fg = color)
            plotSamples(dataJitter, LMrsimple,
                        sampleLabels = sampleLabels, 
                        indexOfReference = indexOfReference,
                        pdfName = "None",
                        joinReplicates = FALSE, showRho = showRho, 
                        xlab = xlab, ylab = ylab, pch = 16,
                        cex = cex, lwd = lwd, 
                        colTst = colTst, colRef = colRef, colRho = colRho,
                        sub = "Restricted Model")
            par(fg = "black")
        }

        ## Unrestricted:
        if (unrestricted) {
            LMusimple <- lm(Response ~ factor(Sample) + Z + factor(Sample):Z, 
                            data = dataframe)
            if (crd) 
                LMu <- LMusimple
            if (latin) 
                LMu <- lm(Response ~ factor(Row) + factor(Column) +
                          factor(Sample) + 
                          Z + factor(Sample):Z, data = dataframe)
            if (blocks) 
                LMu <- lm(Response ~ factor(Replicate) + factor(Sample) + 
                          Z + factor(Sample):Z, data = dataframe)
            color <- ifelse(is.null(tests), "black",
                            ifelse(tests$Parallelity == "Test failed!",
                                   "red", "black"))
            par(fg = color)
            plotSamples(dataJitter, LMusimple,
                        sampleLabels = sampleLabels, 
                        indexOfReference = indexOfReference,
                        pdfName = "None",
                        joinReplicates = FALSE, showRho = FALSE, 
                        xlab = xlab, ylab = ylab, pch = 19,
                        cex = cex, lwd = lwd, 
                        colTst = colTst, colRef = colRef, colRho = colRho,
                        sub = "Unrestricted Model")
            par(fg = "black")
        }

        ## Interaction, twoway, ...:
        if (twoway) {
            LMmsimple <- lm(Response ~ factor(Sample) +
                            Z + factor(Sample):factor(Z), 
                            data = dataframe)
            if (crd) 
                LMm <- LMmsimple
            if (latin) 
                LMm <- lm(Response ~ factor(Row) + factor(Column) +
                          factor(Sample) + 
                          Z + factor(Sample):factor(Z), data = dataframe)
            if (blocks) 
                LMm <- lm(Response ~ factor(Replicate) + factor(Sample) + 
                          Z + factor(Sample):factor(Z), data = dataframe)
            color <- ifelse(is.null(tests), "black",
                            ifelse(tests$Linearity == "Test failed!",
                                   "red", "black"))
            par(fg = color)
            plotSamples(dataJitter, LMmsimple,
                        sampleLabels = sampleLabels, 
                        indexOfReference = indexOfReference,
                        pdfName = "None", joinReplicates = joinReplicates,
                        showRho = FALSE, 
                        xlab = xlab, ylab = ylab, pch = pch,
                        cex = cex, lwd = lwd, 
                        colTst = colTst, colRef = colRef, colRho = colRho, 
                        sub = "Mean graph")
            par(fg = "black")
        }

        ## Mean-adjusted / Stepadjusted Interaction, twoway, ...:
        if (stepadjusted) {
            mn <- mean(unlist(dataJitter["Response"]), na.rm = TRUE)
            AdjustmentStep <- mn - fitted(lm(Response ~ factor(Dilution), 
                                             data = dataJitter))
            stepAdjusted <- dataJitter
            i <- !is.na(stepAdjusted["Response"])
            stepAdjusted[i, "Response"] <-
                dataJitter[i, "Response"] + AdjustmentStep
            LMqsimple <- lm(Response ~ factor(Sample) + Z +
                            factor(Sample):factor(Z), 
                            data = stepAdjusted)
            if (crd) 
                LMq <- LMqsimple
            if (latin) 
                LMq <- lm(Response ~ factor(Row) + factor(Column) +
                          factor(Sample) + 
                          Z + factor(Sample):factor(Z), data = dataframe)
            if (blocks) 
                LMq <- lm(Response ~ factor(Replicate) + factor(Sample) + 
                          Z + factor(Sample):factor(Z), data = dataframe)
            color <- ifelse(is.null(tests), "black",
                            ifelse(tests$Linearity == "Test failed!",
                                   "red", "black"))
            par(fg = color)
            plotSamples(stepAdjusted, LMqsimple,
                        sampleLabels = sampleLabels, 
                        indexOfReference = indexOfReference,
                        pdfName = "None",
                        joinReplicates = joinReplicates, showRho = FALSE,
                        xlab = xlab,
                        ylab =
                        "Respons: Response adjusted by mean for dilutionstep", 
                        pch = pch,
                        cex = cex, lwd = lwd, colTst = colTst,
                        colRef = colRef, colRho = colRho,
                        sub = "Mean-Mean graph")
            par(fg = "black")
        }

        PlateJitter <- as.numeric(unlist(dataJitter["Replicate"])) + 
            (as.numeric(unlist(dataJitter["Dilution"])) - 2) * 0.05
        i <- !is.na(dataJitter["Response"])
        pchI <- pch
        if (length(pch) == length(unlist(dataJitter["Response"])))
            pchI <- pch[i]
        dataF <- cbind(dataJitter[i,],
                       Fitted = fitted(LMr),
                       Residual = residuals(LMr), 
                       PlateJitter = PlateJitter[i])

        ## Qplot:
        if (qplot) {
            qqnorm(residuals(LMr), col = "black", cex = 0.75, cex.lab = 0.75, 
                   cex.main = 0.75, ylab = "Residual, Restricted model")
            qqline(residuals(LMr), col = colNrm)
        }

        ## Histogram:
        if (histogram) {
            m <- mean(residuals(LMr))
            std <- sqrt(var(residuals(LMr)))
            hist(residuals(LMr), density = 20, breaks = 8, prob = TRUE, 
                 xlab = "Residual, Restricted model", cex.lab = 0.75, 
                 main = "")
            x <- 1
            curve(dnorm(x, mean = m, sd = std), col = colNrm, lwd = 2, 
                  add = TRUE)
        }

        ## Residuals by Sample (Box):
        if (residualsBoxSample & any(dimnames(dataF)[[2]] == "Sample")) {
            plot(Residual ~ Sample, pch = pchI,
                 col = c(colRef, colTst)[as.numeric(unlist(dataF["Sample"]))], 
                 xlab = "Sample", ylab = "Residuals, Restricted model", 
                 cex.lab = 0.75, sub = "Restricted Model", data = dataF)
        }


        ## Residuals by SampleStep (Box):
        if (residualsBoxSampleStep &
            any(dimnames(dataF)[[2]] == "SampleStep")) {
            plot(Residual ~ SampleStep, pch = pchI,
                 col = c(colRef, colTst)[as.numeric(unlist(dataF["Sample"]))], 
                 xlab = "Sample * Step", ylab = "Residuals, Restricted model", 
                 cex.lab = 0.75, sub = "Restricted Model", data = dataF)
        }

        ## Residuals by Sample:
        if (residualsSample & any(dimnames(dataF)[[2]] == "Sample")) {
            dataF <- cbind(dataF,
                           SampleNumeric = as.numeric(
                               unlist(dataJitter["Sample"])) + 
                           (as.numeric(
                               unlist(dataJitter["Dilution"])) - 2) * 0.05)
            plot(Residual ~ SampleNumeric, pch = pchI,
                 col = c(colRef, colTst)[as.numeric(unlist(dataF["Sample"]))], 
                 xlab = "Sample", ylab = "Residuals, Restricted model", 
                 cex.lab = 0.75, sub = "Restricted Model", data = dataF)
        }

        ## Residuals by Row:
        if (residualsRow & any(dimnames(dataF)[[2]] == "Row")) {
            dataF <- cbind(dataF,
                           JitterRow = as.numeric(unlist(dataJitter["Row"])) + 
                           (as.numeric(
                               unlist(dataJitter["Dilution"])) - 2) * 0.05)
            plot(Residual ~ JitterRow, pch = pchI,
                 col = c(colRef, colTst)[as.numeric(unlist(dataF["Sample"]))], 
                 xlab = "Row", ylab = "Residuals, Restricted model", 
                 cex.lab = 0.75, sub = "Restricted Model", data = dataF)
        }

        ## Residuals by Column:
        if (residualsColumn & any(dimnames(dataF)[[2]] == "Column")) {
            dataF <- cbind(dataF,
                           JitterColumn = as.numeric(
                               unlist(dataJitter["Column"])) + 
                           (as.numeric(
                               unlist(dataJitter["Dilution"])) - 2) * 0.05)
            plot(Residual ~ JitterColumn, pch = pchI,
                 col = c(colRef, colTst)[as.numeric(unlist(dataF["Sample"]))], 
                 xlab = "Column", ylab = "Residuals, Restricted model", 
                 cex.lab = 0.75, sub = "Restricted Model", data = dataF)
        }

        ## Residuals by block:
        if (residualsBlock & any(dimnames(dataF)[[2]] == "Block")) {
            dataF <- cbind(dataF,
                           JitterBlock = as.numeric(
                               unlist(dataJitter["Block"])) + 
                           (as.numeric(
                               unlist(dataJitter["Dilution"])) - 2) * 0.05)
            plot(Residual ~ JitterBlock, pch = pchI,
                 col = c(colRef, colTst)[as.numeric(unlist(dataF["Sample"]))], 
                 xlab = "Block", ylab = "Residuals, Restricted model", 
                 cex.lab = 0.75, sub = "Restricted Model", data = dataF)
        }

        ## Residuals by Plate / Replicate:
        if (residualsPlate) {
            plot(Residual ~ PlateJitter, pch = pchI,
                 col = c(colRef, colTst)[as.numeric(unlist(dataF["Sample"]))], 
                 xlab = "Replicate", ylab = "Residuals, Restricted model", 
                 cex.lab = 0.75, sub = "Restricted Model", data = dataF)
        }

        ## Residuals by Step:
        if (residualsStep) {
            plot(Residual ~ Zjitter, pch = pchI,
                 col = c(colRef, colTst)[as.numeric(unlist(dataF["Sample"]))], 
                 xlab = xlab, ylab = "Residuals, Restricted model",
                 cex.lab = 0.75, 
                 sub = "Restricted Model", data = dataF)
        }

        title(main = main, outer = TRUE)
        title(sub = "Jens Henrik Badsberg", line = 0, adj = 0.99, 
              cex.sub = 0.5, outer = TRUE)
        title(sub = date(), line = 0, adj = 0.02, cex.sub = 0.5, 
              outer = TRUE)

        invisible(list(data = dataF, lm = LMrsimple))
    }
}
