fit <- function (object, ...) 
    {
        if (any(names(list(...)) == "show"))
            stop(
                "Use \"print\"! - it will also return the vaues for \"fit\"")
        else
            pla.fit(object@data,
                    sampleLabels = object@sampleLabels,
                    indexOfReference = object@indexOfReference,
                    StdName = object@StdName,
                    design = object@design,
                    dfAdj = object@dfAdjustment,
                    dr = object@dilutionRatio,
                    factor = object@factor, # => Can not give as argument!
                    alpha = object@alpha,
                    main = object@assayTitle,
                    show = FALSE,
                    ...,
                    returnPotencyEstimates = TRUE)
    }

pla.fit <-
function (data,
          sampleLabels = levels(unlist(data["Sample"])),
          indexOfReference = 1,
          StdName = sampleLabels[indexOfReference],
          design = "blocks",
          dfAdj = 0,
          dr = 2,
          factor = 1,
          alpha = 0.05,
          main = "Parallel Line Model",
          tag = "PLA",
          expectedAnova = NULL,
          expectedPotency = NULL,
          formatTests = "long",
          show = FALSE,
          sink = FALSE,
          Sweave = FALSE,
          printPotencyEstimates = TRUE,
          returnPotencyEstimates = TRUE) 
{

    DoPlots         <- FALSE
    showHead        <- TRUE
    showAnova       <- TRUE
    showTests       <- TRUE
    showRegressions <- TRUE
    showRatios      <- TRUE
    showPhEur       <- TRUE
    showRsquare     <- TRUE
    showSlope       <- TRUE
    showLogPotency  <- TRUE
    showPotency     <- TRUE
    if (is.logical(show)) {
        Show <- show
    } else {
        Show <- TRUE
        if (is.character(show)) {
            showHead        <- any(show == "all") | any(show == "head")
            showAnova       <- any(show == "all") | any(show == "anova")
            showTests       <- any(show == "all") | any(show == "tests")
            showRegressions <- any(show == "all") | any(show == "regressions")
            showRatios      <- any(show == "all") | any(show == "ratios")
            showPhEur       <- any(show == "all") | any(show == "pheur")
            showRsquare     <- any(show == "all") | any(show == "rsquare")
            showSlope       <- any(show == "all") | any(show == "slope")
            showLogPotency  <- any(show == "all") | any(show == "logpotency")
            showPotency     <- any(show == "all") | any(show == "potency")
            if (showPhEur) {
                showRsquare     <- TRUE
                showSlope       <- TRUE
                showLogPotency  <- TRUE
                showPotency     <- TRUE
            }
        }
    }

    newPage <- function(outer = FALSE) {
        if (Sweave) {
            cat("\\end{Soutput}\n")
            cat("\\end{Schunk}\n")
            cat("\\newpage\n")
            cat("\\begin{Schunk}\n")
            cat("\\begin{Soutput}\n")
        }
    }
    myLine <- function(char = "=", length = 80, start = "", end = "\n") {
        if (Show) {
            cat(start)
            cat(paste(rep(char, 80), collapse = ""))
            cat(end)
        }
    }
    myCat <- function(object) {
        if (Show) {
            cat(object)
        }
    }
    notEmpty <- function(X) {
        ifelse(length(X) > 0, TRUE, FALSE)
    }

    checkTest <- function(Rname, lessThan = FALSE, label = "",
                          name = "", chr = "",
                          ok = "Test passed!",
                          failed = "Test failed!") {
        Alpha <- .extractAlpha(alpha, name)
        AlphaMin <- min(Alpha)
        AlphaMax <- max(Alpha)
        test <- NULL
        test <- anovaTable[Rname, ]
        if (lessThan) 
            passed <- ifelse(test[1, "Pr(>F)"] < AlphaMin, ok,
                             ifelse(test[1, "Pr(>F)"] < AlphaMax,
                                    ifelse(length(Alpha) > 2,
                                           paste0("Unknown: ",
                                                  length(which(
                                                      test[1,
                                                           "Pr(>F)"] > Alpha))),
                                           "Unknown"),
                                    failed))
        else passed <- ifelse(test[1, "Pr(>F)"] > AlphaMax, ok,
                              ifelse(test[1, "Pr(>F)"] > AlphaMin,
                                     ifelse(length(Alpha) > 2,
                                            paste0("Unknown: ",
                                                   length(which(
                                                       test[1,
                                                            "Pr(>F)"] < Alpha))),
                                            "Unknown"),
                                     failed))
        dn <- dimnames(test)[[2]]
        names(passed) <- "Validity"
        test <- cbind(test, qf(1 - ifelse(lessThan, AlphaMin, AlphaMax),
                               test[1, "Df"], df),
                      Passed = passed)
        dimnames(test)[[2]] <- c(dn, "F(critical)", label)
        dimnames(test)[[1]] <- paste("Test of", name, sep = " ")
        if (Show & ((formatTests == "long") |
                    (formatTests == "both"))) {
            if (Show) 
                print(test)
            myLine(chr)
        }
        dimnames(test)[[1]] <- name
        dimnames(test)[[2]] <- c(dn, "F(critical)", "Validity")
        return(test)
    }

    ## if (StdName != "S") 
    ##    warning("Name of 'StdName' is not 'S'")
    selectS <- !is.na(match(as.character(unlist(data["Sample"])), sampleLabels))
    data <- data[selectS, ]
    latin <- FALSE
    blocks <- FALSE
    crd <- FALSE
    design <- .string2design(design)
    if ((design == "lsd"))
        latin <- TRUE
    if ((design == "rbd"))
        blocks <- TRUE
    if (design == "crd") 
        crd <- TRUE
    OK <- .checkPlaFrame(data, design = design,
                         latin = latin, blocks = blocks, crd = crd)
    if (sink) 
        sink(paste(tag, "-Result.txt", sep = ""))
    if (showHead) {
        myLine()
        myCat(paste(c("Dilution ration:", dr, "\n"), sep = " "))
        myCat(paste(c("Factor:", factor, "\n"), sep = " "))
        myCat(paste(c("Significance level alpha:", alpha * 100, "%\n"),
                    sep = " "))
        myLine()
    }
    if (OK) {
        if (showRegressions) {
            if (crd) 
                RestrictedModel <- lm(Response ~ -1 + factor(Sample) + 
                                      Z, data = data)
            if (latin) 
                RestrictedModel <- lm(Response ~ -1 + factor(Sample) + 
                                      factor(Row) + factor(Column) + Z,
                                      data = data)
            if (blocks) 
                RestrictedModel <- lm(Response ~ -1 + factor(Sample) + 
                                      factor(Replicate) + Z, data = data)
            newPage()
            myCat("\n")
            myLine()
            myCat(
                paste0("Regression, Restricted model (Common Slope), ",
                       "with adjusting for 'blocks'.:\n"))
            myLine("-")
            if (Show) 
                print(summary(RestrictedModel))
            myLine("-")
            myCat("Slope:\n")
            if (Show) 
                print(coef(RestrictedModel)["Z"])
            myLine()
            myCat("\n")
            newPage()
            myLine()
            myCat(
                paste0("Regression, Unrestricted Model (Different slopes), ",
                       "with adjusting for 'blocks':\n"))
            myLine("-")
            if (crd) 
                Unrestricted <- lm(Response ~ -1 + factor(Sample) + 
                                   factor(Sample):Z, data = data)
            if (latin) 
                Unrestricted <- lm(Response ~ -1 + factor(Sample) + 
                                   factor(Row) + factor(Column) +
                                   factor(Sample):Z, 
                                   data = data)
            if (blocks) 
                Unrestricted <- lm(Response ~ -1 + factor(Sample) + 
                                   factor(Replicate) + factor(Sample):Z,
                                   data = data)
            if (Show) 
                print(summary(Unrestricted))
            myLine("-")
            myCat("Slopes:\n")
            if (Show) {
                namesCoefs <- names(coef(Unrestricted))
                coefPos <-
                    ((substr(namesCoefs, 1,
                             nchar("factor(Sample)")) == "factor(Sample)") &
                     unlist(lapply(strsplit(namesCoefs, ":"), length)) == 2)
                print(coef(Unrestricted)[coefPos])

            }
            myLine()
            myCat("\n")
        }

        if (showAnova) {
            newPage()
            myCat("\n")
            myLine("-")
        }

        if (crd) 
            FeillerModel <- lm(Response ~ factor(Sample) + Z, data = data)
        if (latin) 
            FeillerModel <- lm(Response ~ factor(Sample) + 
                               factor(Row) + factor(Column) + Z, data = data)
        if (blocks) 
            FeillerModel <- lm(Response ~ factor(Sample) + 
                               factor(Replicate) + Z, data = data)

        if (latin) {
            InteractionModel <- lm(Response ~ factor(Row) + factor(Column) + 
                                   factor(Sample) * Dilution +
                                   factor(Sample):factor(Dilution), 
                                   data = data)
            if (showAnova)
                myCat("Latin Square Design Analysis (LSD):\n")
        }
        if (crd) {
            InteractionModel <- lm(Response ~
                                   factor(Sample) * Dilution +
                                   factor(Sample):factor(Dilution), 
                                   data = data)
            if (showAnova)
                myCat("Completly Randomized Design Analysis (CRD):\n")
        }
        if (blocks) {
            InteractionModel <- lm(Response ~ factor(Replicate) +
                                   factor(Sample) * Dilution +
                                   factor(Sample):factor(Dilution),
                                   data = data)
            if (showAnova)
                myCat("Randomized Block Design Analysis (RBD):\n")
        }

        if (showAnova) {
            myLine("-")
            if (Show) 
                print(anova(InteractionModel))
            myLine("-")
        }

        anovaTable <- anova(InteractionModel)
        class(anovaTable) <- "data.frame"
        pheur <- pheur325(data, lmInteraction = InteractionModel,
                          lmRestricted = FeillerModel, dfAdj = dfAdj,
                          dr = dr, sampleLabels = sampleLabels,
                          StdName = StdName, factor = factor,
                          alpha = .extractAlpha(alpha, "Confidence"))
        {

            if (Show & ((formatTests == "long") |
                        (formatTests == "both"))) {
                newPage()
                myCat("\n")
                myLine("-")
            }

            testResiduals <- anovaTable["Residuals", ]
            df <- testResiduals[1, "Df"]
            test <- anovaTable["factor(Sample)", ]
            dn <- dimnames(test)[[2]]
            dimnames(test)[[1]] <- "Test of Preparation:"
            test <- cbind(test, DFresiduals = df)
            if (Show & ((formatTests == "long") |
                        (formatTests == "both"))) {
                if (Show) 
                    print(test)
                myLine("~")
            }
            passed <- "-"
            names(passed) <- "Validity"
            Anova <- cbind(anovaTable["factor(Sample)", ], NA, passed)
            dimnames(Anova)[[1]] <- "Preparation:"
            dimnames(Anova)[[2]] <- c(dn, "F(critical)", "Validity")

            if (latin) {
                tRow <- checkTest("factor(Row)", FALSE,
                                  "This test passes if F(Blocks) < F(critical)", 
                                  "Rows:", "~",
                                  ok = "(Test passed)", failed = "(Test failed)")
                Anova <- rbind(Anova, tRow)
                tCol <- checkTest("factor(Column)", FALSE,
                                  "This test passes if F(Blocks) < F(critical)", 
                                  "Columns:", "~",
                                  ok = "(Test passed)", failed = "(Test failed)")
                Anova <- rbind(Anova, tCol)
            }

            if (blocks) {
                tBlc <- checkTest("factor(Replicate)", FALSE,
                                  "This test passes if F(Blocks) < F(critical)", 
                                  "Blocks:", "~",
                                  ok = "(Test passed)", failed = "(Test failed)")
                Anova <- rbind(Anova, tBlc)
            }
            tReg <- checkTest("Dilution", TRUE,
                              "This test passes if F(regression) > F(critical)", 
                              "Regression:", "~")
            Anova <- rbind(Anova, tReg)
            OkReg <- tReg[, 7]
            tLin <- checkTest("factor(Sample):factor(Dilution)", FALSE,
                        "This test passes if F(non-linearity) < F(critical)", 
                              "Linearity:", "~")
            Anova <- rbind(Anova, tLin)
            OkLin <- tLin[, 7]
            tPar <- checkTest("factor(Sample):Dilution", FALSE,
                        "This test passes if F(non-parallelity) < F(critical)", 
                              "Parallelism:",
                              ifelse((formatTests == "long"), "=", " "))
            Anova <- rbind(Anova, tPar)
            OkPar <- tPar[, 7]

            if (Show & ((formatTests == "long") |
                        (formatTests == "BOTH"))) {
                ## newPage()
                ## myCat("\n")
            }

            if (Show & ((formatTests == "short") |
                        (formatTests == "both"))) {
                if (formatTests == "both")
                    myLine("+")
                print(Anova[, -6])
                myLine()
                ## newPage()
                ## myCat("\n")
            }
        }

        {
            SS <- unlist(pheur$SS)
            reg <- unlist(pheur$reg)
            K <- pheur$K
            KP <- pheur$KP
        }

        Tests <- list(Regression = OkReg, Linearity = OkLin,
                      Parallelity = OkPar)
        Fits <- list(Tests = Tests, AT = anovaTable, SS = SS, reg = reg, 
                     K = K, KP = KP)
        if (notEmpty(expectedPotency)) {
            de <- dim(expectedPotency)
            e1 <- de[1]
            dk <- dim(KP)
            if (is.null(dk))
                KPX <- matrix(rep(KP[2:4], rep(e1, 3)), ncol = 3) else {
                    ## when length(factor) == 1:
                    k1 <- dk[1]
                    ## cat(paste(" --- fit pla: ", e1, k1, e1/k1, " --- \n"))
                    KPX <- matrix(rep(c(KP[,2:4]), rep(e1/k1, k1*3)), ncol = 3)
                }
            relativePotency <- KPX / expectedPotency
        }
        else relativePotency <- matrix(ncol = 0, nrow = 0)
        if (notEmpty(expectedAnova)) {
            relativeAnova <-
                rbind(as.numeric(anovaTable["factor(Sample)", ]                /
                                 expectedAnova["Preparations",    ]),
                      as.numeric(anovaTable["Dilution", ]                      /
                                 expectedAnova["Regression",      ]),
                      as.numeric(anovaTable["factor(Sample):Dilution", ]       /
                                 expectedAnova["Non-parallelism", ]),
                      as.numeric(anovaTable["factor(Sample):factor(Dilution)",]/
                                 expectedAnova["Non-linearity",   ]), 
                      as.numeric(anovaTable["Residual", ]                      /
                                 expectedAnova["Residual_error", ] ))
        }
        else relativeAnova <- matrix(ncol = 0, nrow = 0)

        if (showRatios) {
            newPage()
            if (notEmpty(expectedAnova)) {
                myLine()
                myCat("Ratio of ANOVA tables:\n")
                if (Show) 
                    print(relativeAnova)
                myLine()
                myCat("\n")
            }

            if (notEmpty(expectedPotency)) {
                myLine()
                myCat("Ratio of Potency tables:\n")
                if (Show) 
                    print(relativePotency)
                myLine()
                myCat("\n")
            }
        }

        if (showPhEur) {
            myLine()
            myCat(
                paste0("Sums of Squares, Slope, Variance, C, ",
                       "Potency and Confidence-intervals:\n"))
            myLine("-")
            if (printPotencyEstimates & Show) {
                digits <- 8
                if (showRsquare)
                    print(SS,  digits = digits)
                if (showSlope)
                    print(reg, digits = digits)
                if (showLogPotency)
                    print(K,   digits = digits)
                if (showPotency)
                    print(KP,  digits = digits)
                myLine("-")
            }
            myLine()
        }

        if (DoPlots) {
            pla.plots(data, pdfName = paste(tag, "-Plots.pdf", sep = ""),
                      sampleLabels = sampleLabels, 
                      design = design, main = main, tests = Tests)
        }
        Args <- list(data = data, sampleLabels = sampleLabels,
                     design = design, dfAdj = dfAdj, dr = dr,
                     main = main, alpha = alpha, factor = factor)
        if (returnPotencyEstimates) 
            invisible(new("plaFit", inpArgs = Args, design = design, 
                          lm = InteractionModel, pheur = pheur, anova = Anova,
                          tests = Tests, relAnova = relativeAnova, 
                          relPotency = relativePotency))
    }
}
