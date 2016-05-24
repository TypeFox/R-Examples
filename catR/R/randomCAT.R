randomCAT<-function (trueTheta, itemBank, model = NULL, responses = NULL, 
    genSeed = NULL, maxItems = 50, cbControl = NULL, nAvailable = NULL, 
    start = list(fixItems = NULL, seed = NULL, nrItems = 1, theta = 0, 
        D = 1, randomesque = 1, startSelect = "MFI"), test = list(method = "BM", 
        priorDist = "norm", priorPar = c(0, 1), range = c(-4, 
            4), D = 1, parInt = c(-4, 4, 33), itemSelect = "MFI", 
        infoType = "observed", randomesque = 1, AP = 1, constantPatt = NULL), 
    stop = list(rule = "length", thr = 20, alpha = 0.05), final = list(method = "BM", 
        priorDist = "norm", priorPar = c(0, 1), range = c(-4, 
            4), D = 1, parInt = c(-4, 4, 33), alpha = 0.05), 
    allTheta = FALSE, save.output = FALSE, output = c("path", 
        "name", "csv")) 
{
    if (missing(trueTheta)) {
        if (is.null(responses)) 
            stop("'trueTheta' was not provided!", call. = FALSE)
        trueTheta <- NA
    }
    else trueTheta <- trueTheta
    if (!testList(start, type = "start")$test) 
        stop(testList(start, type = "start")$message, call. = FALSE)
    if (!testList(test, type = "test")$test) 
        stop(testList(test, type = "test")$message, call. = FALSE)
    if (!testList(stop, type = "stop")$test) 
        stop(testList(stop, type = "stop")$message, call. = FALSE)
    if (!testList(final, type = "final")$test) 
        stop(testList(final, type = "final")$message, call. = FALSE)
    if (!is.null(model) & !is.null(test$constantPatt)) 
        stop("Treating constant patterns with specific adjustment not yet implemented with polytomous models", 
            call. = FALSE)
    if (!is.null(responses)) 
        assigned.responses <- TRUE
    else assigned.responses <- FALSE
    if (!is.null(cbControl)) {
        prov <- breakBank(itemBank)
        itemBank <- prov$itemPar
        cbGroup <- prov$cbGroup
        if (!test.cbList(cbControl, cbGroup)$test) 
            stop(test.cbList(cbControl, cbGroup)$message, call. = FALSE)
    }
    else {
        cbGroup <- NULL
        itemBank <- as.matrix(itemBank)
    }
    internalCAT <- function() {
        startList <- list(fixItems = start$fixItems, seed = start$seed, 
            nrItems = NULL, theta = start$theta, D = 1, randomesque = 1, 
            startSelect = "MFI")
        startList$nrItems <- ifelse(is.null(start$nrItems), 1, 
            start$nrItems)
        if(is.null(start$theta)) startList$theta <- 0
        startList$D <- ifelse(is.null(start$D), 1, start$D)
        startList$randomesque <- ifelse(is.null(start$randomesque ), 
            1, start$randomesque )
        startList$startSelect <- ifelse(is.null(start$startSelect), 
            "MFI", start$startSelect)
        start <- startList
        stopList <- list(rule = NULL, thr = 20, alpha = 0.05)
        stopList$rule <- ifelse(is.null(stop$rule), "length", 
            stop$rule)
        stopList$thr <- ifelse(is.null(stop$thr), 20, stop$thr)
        stopList$alpha <- ifelse(is.null(stop$alpha), 0.05, stop$alpha)
        stop <- stopList
        if (stop$rule == "length" & stop$thr > maxItems) 
            stop("'maxItems' is smaller than test length criterion 'stop$thr'", 
                call. = FALSE)
        testList <- list(method = NULL, priorDist = NULL, priorPar = c(0, 
            1), range = c(-4, 4), D = 1, parInt = c(-4, 4, 33), 
            itemSelect = "MFI", infoType = "observed", randomesque = 1, 
            AP = 1, rule = stop$rule, thr = stop$thr, constantPatt = NULL)
        testList$method <- ifelse(is.null(test$method), "BM", 
            test$method)
        testList$priorDist <- ifelse(is.null(test$priorDist), 
            "norm", test$priorDist)
        if (!is.null(test$priorPar)) {
            testList$priorPar[1] <- test$priorPar[1]
            testList$priorPar[2] <- test$priorPar[2]
        }
        if (!is.null(test$range)) {
            testList$range[1] <- test$range[1]
            testList$range[2] <- test$range[2]
        }
        testList$D <- ifelse(is.null(test$D), 1, test$D)
        if (!is.null(test$parInt)) {
            testList$parInt[1] <- test$parInt[1]
            testList$parInt[2] <- test$parInt[2]
            testList$parInt[3] <- test$parInt[3]
        }
        testList$itemSelect <- ifelse(is.null(test$itemSelect), 
            "MFI", test$itemSelect)
        testList$infoType <- ifelse(is.null(test$infoType), "observed", 
            test$infoType)
        testList$randomesque <- ifelse(is.null(test$randomesque), 
            1, test$randomesque)
        testList$AP <- ifelse(is.null(test$AP), 1, test$AP)
        if (!is.null(test$constantPatt)) 
            testList$constantPatt <- test$constantPatt
        test <- testList
        finalList <- list(method = NULL, priorDist = NULL, priorPar = c(0, 
            1), range = c(-4, 4), D = 1, parInt = c(-4, 4, 33), 
            alpha = 0.05)
        finalList$method <- ifelse(is.null(final$method), "BM", 
            final$method)
        finalList$priorDist <- ifelse(is.null(final$priorDist), 
            "norm", final$priorDist)
        if (is.null(final$priorPar) == FALSE) {
            finalList$priorPar[1] <- final$priorPar[1]
            finalList$priorPar[2] <- final$priorPar[2]
        }
        if (!is.null(final$range)) {
            finalList$range[1] <- final$range[1]
            finalList$range[2] <- final$range[2]
        }
        finalList$D <- ifelse(is.null(final$D), 1, final$D)
        if (!is.null(final$parInt)) {
            finalList$parInt[1] <- final$parInt[1]
            finalList$parInt[2] <- final$parInt[2]
            finalList$parInt[3] <- final$parInt[3]
        }
        finalList$alpha <- ifelse(is.null(final$alpha), 0.05, 
            final$alpha)
        final <- finalList
        if (stop$rule == "length") 
            maxLength <- min(c(maxItems, stop$thr, nrow(itemBank)))
        else maxLength <- min(maxItems, nrow(itemBank))
        if (stop$rule == "classification" & (test$itemSelect == 
            "progressive" | test$itemSelect == "proportional")) 
            stop("'classification' rule cannot be considered with neither 'progressive' nor 'proportional' item selection rules!", 
                call. = FALSE)
        pr0 <- startItems(itemBank = itemBank, model = model, 
            fixItems = start$fixItems, seed = start$seed, nrItems = start$nrItems, 
            theta = start$theta, D = start$D, randomesque = start$randomesque, 
            startSelect = start$startSelect, nAvailable = nAvailable)
        ITEMS <- pr0$items
        PAR <- rbind(pr0$par)
        if (is.null(ITEMS)) 
            PATTERN <- NULL
        else {
            if (!is.null(responses)) 
                PATTERN <- responses[ITEMS]
            else PATTERN <- genPattern(trueTheta, PAR, model = model, 
                D = test$D, seed = genSeed)
        }
        if (!is.null(ITEMS)) 
            TH <- thetaEst(PAR, PATTERN, model = model, D = test$D, 
                method = test$method, priorDist = test$priorDist, 
                priorPar = test$priorPar, range = test$range, 
                parInt = test$parInt, current.th = start$theta, 
                constantPatt = test$constantPatt, bRange = range(itemBank[, 
                  2]))
        else TH <- start$theta
        if (!is.null(ITEMS)) 
            SETH <- semTheta(TH, PAR, x = PATTERN, model = model, 
                D = test$D, method = test$method, priorDist = test$priorDist, 
                priorPar = test$priorPar, parInt = test$parInt, 
                constantPatt = test$constantPatt)
        else SETH <- NA
        thProv <- TH
        enter.cat <- TRUE
        if (!is.na(SETH)) {
            if (length(PATTERN) >= maxLength | (stop$rule == 
                "precision" & SETH <= stop$thr) | (stop$rule == 
                "classification" & (TH - qnorm(1 - stop$alpha/2) * 
                SETH >= stop$thr | TH + qnorm(1 - stop$alpha/2) * 
                SETH <= stop$thr))) 
                enter.cat <- FALSE
        }
        if (!enter.cat) {
            finalEst <- thetaEst(PAR, PATTERN, model = model, 
                D = final$D, method = final$method, priorDist = final$priorDist, 
                priorPar = final$priorPar, range = final$range, 
                parInt = final$parInt)
            seFinal <- semTheta(finalEst, PAR, x = PATTERN, model = model, 
                D = final$D, method = final$method, priorDist = final$priorDist, 
                priorPar = final$priorPar, parInt = final$parInt)
            confIntFinal <- c(finalEst - qnorm(1 - final$alpha/2) * 
                seFinal, finalEst + qnorm(1 - final$alpha/2) * 
                seFinal)
            endWarning <- FALSE
            RES <- list(trueTheta = trueTheta, model = model, 
                maxItems = maxItems, testItems = ITEMS, itemPar = PAR, 
                pattern = PATTERN, thetaProv = TH, seProv = SETH, 
                thFinal = finalEst, seFinal = seFinal, ciFinal = confIntFinal, 
                genSeed = genSeed, startFixItems = start$fixItems, 
                startSeed = start$seed, startNrItems = start$nrItems, 
                startTheta = start$theta, startD = start$D, startRandomesque = start$randomesque, 
                startSelect = start$startSelect, 
                provMethod = test$method, provDist = test$priorDist, 
                provPar = test$priorPar, provRange = test$range, 
                provD = test$D, itemSelect = test$itemSelect, 
                infoType = test$infoType, randomesque = test$randomesque, 
                AP = test$AP, constantPattern = test$constantPatt, 
                cbControl = cbControl, cbGroup = cbGroup, stopRule = stop$rule, 
                stopThr = stop$thr, stopAlpha = stop$alpha, endWarning = endWarning, 
                finalMethod = final$method, finalDist = final$priorDist, 
                finalPar = final$priorPar, finalRange = final$range, 
                finalD = final$D, finalAlpha = final$alpha, save.output = save.output, 
                output = output, assigned.responses = assigned.responses)
            class(RES) <- "cat"
        }
        else {
            repeat {
                pr <- nextItem(itemBank, model = model, theta = thProv, 
                  out = ITEMS, x = PATTERN, criterion = test$itemSelect, 
                  method = test$method, parInt = test$parInt, 
                  priorDist = test$priorDist, priorPar = test$priorPar, 
                  infoType = test$infoType, D = test$D, range = test$range, 
                  randomesque = test$randomesque, AP = test$AP, 
                  cbControl = cbControl, cbGroup = cbGroup, maxItems = maxItems, 
                  rule = test$rule, thr = test$thr, SETH = SETH[length(SETH)], 
                  nAvailable = nAvailable)
                ITEMS <- c(ITEMS, pr$item)
                PAR <- rbind(PAR, pr$par)
                if (!is.null(responses)) 
                  PATTERN <- c(PATTERN, responses[pr$item])
                else PATTERN <- c(PATTERN, genPattern(trueTheta, 
                  pr$par, model = model, D = test$D, seed = genSeed))
                thProv <- thetaEst(PAR, PATTERN, model = model, 
                  D = test$D, method = test$method, priorDist = test$priorDist, 
                  priorPar = test$priorPar, range = test$range, 
                  parInt = test$parInt, current.th = TH[length(TH)], 
                  constantPatt = test$constantPatt, bRange = range(itemBank[, 
                    2]))
                TH <- c(TH, thProv)
                seProv <- semTheta(thProv, PAR, x = PATTERN, 
                  model = model, D = test$D, method = test$method, 
                  priorDist = test$priorDist, priorPar = test$priorPar, 
                  parInt = test$parInt, constantPatt = test$constantPatt)
                SETH <- c(SETH, seProv)
                if ((length(ITEMS) >= maxLength) | (stop$rule == 
                  "precision" & seProv <= stop$thr) | (stop$rule == 
                  "classification" & (thProv - qnorm(1 - stop$alpha/2) * 
                  seProv >= stop$thr | thProv + qnorm(1 - stop$alpha/2) * 
                  seProv <= stop$thr))) 
                  break
            }
            finalEst <- thetaEst(PAR, PATTERN, model = model, 
                D = final$D, method = final$method, priorDist = final$priorDist, 
                priorPar = final$priorPar, range = final$range, 
                parInt = final$parInt)
            seFinal <- semTheta(finalEst, PAR, x = PATTERN, model = model, 
                D = final$D, method = final$method, priorDist = final$priorDist, 
                priorPar = final$priorPar, parInt = final$parInt)
            confIntFinal <- c(finalEst - qnorm(1 - final$alpha/2) * 
                seFinal, finalEst + qnorm(1 - final$alpha/2) * 
                seFinal)
            if ((stop$rule == "length" & length(ITEMS) < stop$thr) | 
                (stop$rule == "precision" & seProv > stop$thr) | 
                (stop$rule == "classification" & thProv - qnorm(1 - 
                  stop$alpha/2) * seProv < stop$thr & thProv + 
                  qnorm(1 - stop$alpha/2) * seProv > stop$thr)) 
                endWarning <- TRUE
            else endWarning <- FALSE
            RES <- list(trueTheta = trueTheta, model = model, 
                maxItems = maxItems, testItems = ITEMS, itemPar = PAR, 
                pattern = PATTERN, thetaProv = TH, seProv = SETH, 
                thFinal = finalEst, seFinal = seFinal, ciFinal = confIntFinal, 
                genSeed = genSeed, startFixItems = start$fixItems, 
                startSeed = start$seed, startNrItems = start$nrItems, 
                startTheta = start$theta, startD = start$D, startRandomesque = start$randomesque, 
                startSelect = start$startSelect, 
                provMethod = test$method, provDist = test$priorDist, 
                provPar = test$priorPar, provRange = test$range, 
                provD = test$D, itemSelect = test$itemSelect, 
                infoType = test$infoType, randomesque = test$randomesque, 
                AP = test$AP, constantPattern = test$constantPatt, 
                cbControl = cbControl, cbGroup = cbGroup, stopRule = stop$rule, 
                stopThr = stop$thr, stopAlpha = stop$alpha, endWarning = endWarning, 
                finalMethod = final$method, finalDist = final$priorDist, 
                finalPar = final$priorPar, finalRange = final$range, 
                finalD = final$D, finalAlpha = final$alpha, save.output = save.output, 
                output = output, assigned.responses = assigned.responses)
            class(RES) <- "cat"
        }
        if (allTheta & start$nrItems > 0) {
            nra <- length(RES$pattern) - length(RES$thetaProv)
            if (nra > 0) {
                prov.th <- prov.se <- NULL
                for (k in 1:nra) {
                  prov.par <- rbind(RES$itemPar[1:k, ])
                  prov.th[k] <- thetaEst(prov.par, RES$pattern[1:k], 
                    model = model, D = test$D, method = test$method, 
                    priorDist = test$priorDist, priorPar = test$priorPar, 
                    range = test$range, parInt = test$parInt, 
                    constantPatt = test$constantPatt, bRange = range(itemBank[, 
                      2]))
                  prov.se[k] <- semTheta(prov.th[k], prov.par, 
                    RES$pattern[1:k], model = model, D = test$D, 
                    method = test$method, priorDist = test$priorDist, 
                    priorPar = test$priorPar, parInt = test$parInt, 
                    constantPatt = test$constantPatt)
                }
                RES$thetaProv <- c(prov.th, RES$thetaProv)
                RES$seProv <- c(prov.se, RES$seProv)
            }
        }
        return(RES)
    }
    resToReturn <- internalCAT()
    if (save.output) {
        if (output[1] == "path") 
            wd <- paste(getwd(), "/", sep = "")
        else wd <- output[1]
        if (output[3] == "csv") 
            fileName <- paste(wd, output[2], ".csv", sep = "")
        else fileName <- paste(wd, output[2], ".txt", sep = "")
        capture.output(resToReturn, file = fileName)
    }
    return(resToReturn)
}






print.cat<-function (x, ...) 
{
  if (!x$assigned.responses) {
cat("Random generation of a CAT response pattern", "\n")
if (is.null(x$genSeed)) cat("  without fixing the random seed", "\n", "\n")
else cat("  with random seed equal to", x$genSeed, "\n", "\n")
}

  else cat("Post-hoc simulation of a full bank provided response pattern", "\n","\n")
  if (is.null(x$model)){
    if (min(x$itemPar[,4])<1) mod<-"Four-Parameter Logistic model"
    else{
      if (max(x$itemPar[,3])>0) mod<-"Three-Parameter Logistic model"
      else{
        if (length(unique(x$itemPar[,1]))>1) mod<-"Two-Parameter Logistic model"
        else 
          mod<-"One-Parameter Logistic (Rasch) model"
      }
    }
  }
  else{
    if (x$model=="GRM") mod<-"Graded Response Model"
    if (x$model=="MGRM") mod<-"Modified Graded Response Model"
    if (x$model=="PCM") mod<-"Partial Credit Model"
    if (x$model=="GPCM") mod<-"Generalized Partial Credit Model"
    if (x$model=="RSM") mod<-"Rating Scale Model"
    if (x$model=="NRM") mod<-"Nominal Response Model"
  }
  cat(" Item bank calibrated under", mod,"\n","\n")
  if (!is.na(x$trueTheta)) cat(" True ability level:", round(x$trueTheta, 2), "\n", "\n")
else cat(" True ability level was not provided", "\n", "\n")
  cat(" Starting parameters:", "\n")
if (x$startNrItems==0){
cat("   No early item was selected","\n")
cat("   Starting ability level:", round(x$startTheta,3),"\n")
}
else{
  if (is.null(x$startFixItems)) {
if (!is.null(x$startSeed)) nr1<-x$startNrItems
    else nr1 <- length(x$startTheta)
}
  else nr1 <- length(x$startFixItems)
  if (x$startSelect == "progressive" | x$startSelect == "proportional")   
    cat("   Number of early items:", 1, "\n")
  else
    cat("   Number of early items:", nr1, "\n")
  if (!is.null(x$startFixItems)) 
    met1 <- "Chosen by administrator"
  else {
    if (!is.null(x$startSeed)) 
      met1 <- "Random selection in item bank"
    else {
      if (x$startSelect == "bOpt") {
        if (x$startNrItems == 1) 
          met1 <- "matching item difficulty to starting ability"
        else met1 <- "matching item difficulties to starting abilities"
      }
      if (x$startSelect == "thOpt") {
        if (x$startNrItems == 1) 
          met1 <- "matching ability level with maximum information to starting ability"
        else met1 <- "matching ability level with maximum information to starting abilities"
      }
      if (x$startSelect == "progressive" | x$startSelect == "proportional") 
                  met1 <- "random selection of the first item"
      if (x$startSelect == "MFI") {
        if (x$startNrItems == 1) 
          met1 <- "maximum informative item for starting ability"
        else met1 <- "maximum informative items for starting abilities"
      }
    }
  }
  if (nr1 == 1) 
    cat("   Early item selection:", met1, "\n")
  else cat("   Early items selection:", met1, "\n")
  if (!is.null(x$startFixItems)) {
    if (length(x$startFixItems) == 1) 
      met1bis <- paste("   Item administered: ", x$startFixItems, sep = "")
    else {
      met1bis <- paste("   Items administered: ", x$startFixItems[1], sep = "")
      if (length(x$startFixItems) == 2) 
        met1bis <- paste(met1bis, " and ", x$startFixItems[2], sep = "")
      else {
        for (i in 2:(length(x$startFixItems) - 1)) met1bis <- paste(met1bis, ", ", x$startFixItems[i], sep = "")
        met1bis <- paste(met1bis, " and ", x$startFixItems[length(x$startFixItems)], sep = "")
      }
    }
    cat(met1bis, "\n")
  }
  if (is.null(x$startFixItems) & !is.null(x$startSeed)) {
    if (x$startNrItems == 1) 
      met1bis <- paste("   Item administered: ", x$testItems[1], sep = "")
    else {
      met1bis <- paste("   Items administered: ", x$testItems[1], sep = "")
      if (x$startNrItems == 2) 
        met1bis <- paste(met1bis, " and ", x$testItems[2], sep = "")
      else {
        for (i in 2:(x$startNrItems - 1)) met1bis <- paste(met1bis, ", ", x$testItems[i], sep = "")
        met1bis <- paste(met1bis, " and ", x$testItems[x$startNrItems], sep = "")
      }
    }
    cat(met1bis, "\n")
  }
  if (is.null(x$startFixItems) & is.null(x$startSeed)) {
    retain <- x$startTheta
    x$startTheta <- round(x$startTheta,3)
    if (x$startSelect != "progressive" & x$startSelect != "proportional") {  
      if (length(x$startTheta) == 1) 
        met1bis <- paste("   Starting ability: ", x$startTheta, sep = "")
      else {
        met1bis <- paste("   Starting abilities: ", sort(x$startTheta)[1], sep = "")
        if (length(x$startTheta) == 2) 
          met1bis <- paste(met1bis, " and ", sort(x$startTheta)[2], sep = "")
        else {
          for (i in 2:(length(x$startTheta) - 1)) met1bis <- paste(met1bis, ", ", sort(x$startTheta)[i], sep = "")
          met1bis <- paste(met1bis, " and ", sort(x$startTheta)[length(x$startTheta)], sep = "")
        }
      }
    cat(met1bis, "\n")
    }
    x$startTheta <- retain
  }
}
  cat("\n", "Adaptive test parameters:", "\n")
  itemSel <- switch(x$itemSelect, MFI = "maximum Fisher information", 
                    Urry = "Urry's procedure", MLWI = "Maximum likelihood weighted information (MLWI)", 
                    MPWI = "Maximum posterior weighted information (MPWI)", 
                    MEI = "Maximum expected information (MEI)", MEPV = "Minimum Expected Posterior Variance (MEPV)", 
                    random = "Random selection", progressive = "Progressive method",
                    proportional = "Proportional method", thOpt = "Optimal theta selection",
                    KL = "Kullback-Leibler (KL) information", KLP = "Posterior Kullback-Leibler (KLP) information",
                    GDI="Global discrimination index (GDI)", GDIP = "Posterior global discrimination index (GDIP)")
  cat("   Next item selection method:", itemSel, "\n")
if (x$itemSelect=="proportional" | x$itemSelect=="progressive")
  cat("     Acceleration parameter for", x$itemSelect,"method:", x$AP, "\n")
if (x$itemSelect=="KLP" | x$itemSelect=="MPWI"){
   met3 <- switch(x$provDist, norm = paste("N(", round(x$provPar[1],2), ",", round(x$provPar[2]^2, 2), ") prior", sep = ""), 
                   unif = paste("U(", round(x$provPar[1], 2), ",", round(x$provPar[2], 2), ") prior", sep = ""), Jeffreys = "Jeffreys' prior")
  cat("     Prior ability distribution for", x$itemSelect, "method:", met3, "\n")
}
  if (x$itemSelect == "MEI" & is.null(x$model)) {
    infTyp <- switch(x$infoType, observed = "observed information function", 
                     Fisher = "Fisher information function")
    cat("     Type of information:", infTyp, "\n")
  }
  met2 <- switch(x$provMethod, BM = "Bayes modal (MAP) estimator", 
                 WL = "Weighted likelihood estimator", ML = "Maximum likelihood estimator", 
                 EAP = "Expected a posteriori (EAP) estimator")
  if (x$provMethod == "BM" | x$provMethod == "EAP") {
    met3 <- switch(x$provDist, norm = paste("N(", round(x$provPar[1],2), ",", round(x$provPar[2]^2, 2), ") prior", sep = ""), 
                   unif = paste("U(", round(x$provPar[1], 2), ",", round(x$provPar[2], 2), ") prior", sep = ""), Jeffreys = "Jeffreys' prior")
  }
  if (x$provMethod == "ML") 
    ra1 <- paste("[", round(x$provRange[1], 2), ",", round(x$provRange[2], 
                                                           2), "]", sep = "")
  cat("   Provisional ability estimator:", met2, "\n")
  if (x$provMethod == "BM" | x$provMethod == "EAP") 
    cat("     Provisional prior ability distribution:", met3, 
        "\n")
  if (x$provMethod == "ML") 
    cat("   Provisional range of ability values:", ra1, "\n")
if (!is.null(x$model) | is.null(x$constantPattern)) adj<-"none"
else adj<-switch(x$constantPattern,fixed4="fixed .4 stepsize",
fixed7="fixed .7 stepsize",
var="variable stepsize")
 cat("   Ability estimation adjustment for constant pattern:", adj, "\n")
  cat("\n")
  cat(" Stopping rule:", "\n")
  met4 <- switch(x$stopRule, length = "length of test", precision = "precision of ability estimate", 
                 classification = paste("classification based on ", 100 * 
                                          (1 - x$stopAlpha), "% confidence interval", sep = ""))
  cat("   Stopping criterion:", met4, "\n")
  switch(x$stopRule, precision = cat("   Maximum SE value:", 
                                     round(x$stopThr, 2), "\n"), classification = cat("   Classification threshold:", 
                                                                                      round(x$stopThr, 2), "\n"))
  cat("   Maximum test length:", ifelse(x$stopRule == "length", 
                                        min(c(x$maxItems, x$stopThr)), x$maxItems), "items", 
      "\n")
  cat("\n", "Randomesque method:", "\n")
 if (is.null(x$startFixItems) & is.null(x$startSeed))
cat("   Number of 'randomesque' starting items: ",x$startRandomesque,"\n",sep="")
else cat("   Number of 'randomesque' starting items: irrelevant","\n",sep="")
  if (x$startSelect == "progressive" | x$startSelect == "proportional")     
    cat("   Number of 'randomesque' test items: ",1,"\n",sep="")
  else cat("   Number of 'randomesque' test items: ",x$randomesque,"\n",sep="")
  cat("\n", "Content balancing control:", "\n")
  if (is.null(x$cbControl)) cat("   No control for content balancing","\n")
  else{
    cat("   Expected proportions of items per subgroup:","\n","\n")
    mat<-rbind(round(x$cbControl$props/sum(x$cbControl$props),3))
    rownames(mat)<-""
    colnames(mat)<-x$cbControl$names
    print(format(mat, justify = "right"), quote = FALSE)
    cat("\n")
  }
  cat("\n", "Adaptive test details:", "\n")
  mat <- rbind(as.character(1:length(x$testItems)), as.character(x$testItems), 
               round(x$pattern, 0))
  nra <- length(x$pattern) - length(x$thetaProv)
if (nra<0) mat <-rbind(mat, c(round(x$thetaProv[2:length(x$thetaProv)], 3)), 
               c(round(x$seProv[2:length(x$seProv)], 3)))
else 
  mat <- rbind(mat, c(rep(NA, nra), round(x$thetaProv, 3)), 
               c(rep(NA, nra), round(x$seProv, 3)))
  rownames(mat) <- c("Nr", "Item", "Resp.", "Est.", "SE")
  colnames(mat) <- rep("", ncol(mat))
  if (x$startSelect != "progressive" & x$startSelect != "proportional") {
    if (nra==0 & x$startNrItems>1) {
      numb <- x$startNrItems-1
      for (i in 4:5) {
        for (j in 1:numb) {
          mat[i,j]<-paste("(",mat[i,j],")",sep="")
        }
      }
    }
  }
  print(format(mat, justify = "right"), quote = FALSE)
  cat("\n")
  if (x$endWarning) 
    cat("WARNING: stopping rule was not satisfied after", 
        length(x$pattern), "items!", "\n", "\n")
  cat(" Final results:", "\n")
  met <- switch(x$finalMethod, BM = "Bayes modal (MAP) estimator", 
                WL = "Weighted likelihood estimator", ML = "Maximum likelihood estimator", 
                EAP = "Expected a posteriori (EAP) estimator")
  if (x$finalMethod == "BM" | x$finalMethod == "EAP") {
    met2 <- switch(x$finalDist, norm = paste("N(", round(x$finalPar[1], 
                                                         2), ",", round(x$finalPar[2]^2, 2), ") prior", sep = ""), 
                   unif = paste("U(", round(x$finalPar[1], 2), ",", 
                                round(x$finalPar[2], 2), ") prior", sep = ""), 
                   Jeffreys = "Jeffreys' prior")
  }
  if (x$finalMethod == "ML") 
    ra1 <- paste("[", round(x$finalRange[1], 2), ",", round(x$finalRange[2], 
                                                            2), "]", sep = "")
  cat("   Length of adaptive test:", length(x$testItems), "items", 
      "\n")
  cat("   Final ability estimator:", met, "\n")
  if (x$finalMethod == "BM" | x$finalMethod == "EAP") 
    cat("   Final prior distribution:", met2, "\n")
  if (x$finalMethod == "ML") 
    cat("   Final range of ability values:", ra1, "\n")
  cat("   Final ability estimate (SE):", round(x$thFinal, 3), 
      paste("(", round(x$seFinal, 3), ")", sep = ""), "\n")
  cat(paste("   ", (1 - x$finalAlpha) * 100, "% confidence interval: [", 
            round(x$ciFinal[1], 3), ",", round(x$ciFinal[2], 3), 
            "]", sep = ""), "\n")
  if (x$stopRule == "classification") {
    if (x$ciFinal[1] > x$stopThr) 
      mess <- paste("ability is larger than ", round(x$stopThr, 
                                                     2), sep = "")
    else {
      if (x$ciFinal[2] < x$stopThr) 
        mess <- paste("ability is smaller than ", round(x$stopThr, 
                                                        2), sep = "")
      else mess <- paste("ability is not different from ", 
                         round(x$stopThr, 2), sep = "")
    }
    cat("   Final subject classification:", mess, "\n")
  }
  if (!is.null(x$cbControl)){
    cat("   Proportions of items per subgroup (expected and observed):","\n","\n")
    mat<-rbind(round(x$cbControl$props/sum(x$cbControl$props),3))
    nr<-NULL
    for (i in 1:length(x$cbControl$names)) nr[i]<-length(x$testItems[x$cbGroup[x$testItems]==x$cbControl$names[i]])
    nr<-nr/sum(nr)
    mat<-rbind(mat,round(nr,3))
    rownames(mat)<-c("Exp.","Obs.")
    colnames(mat)<-x$cbControl$names
    print(format(mat, justify = "right"), quote = FALSE)
    cat("\n")
    cat("   Items administered per subgroup:","\n","\n")
    for (i in 1:length(x$cbControl$names)) {
      if (length(x$testItems[x$cbGroup[x$testItems]==x$cbControl$names[i]])==0) mess<-"none"
      else{
        its<-sort(x$testItems[x$cbGroup[x$testItems]==x$cbControl$names[i]])
        mess<-its[1]
        if (length(its)>1){
          for (j in 2:length(its)) mess<-paste(mess,", ",its[j],sep="")
        }
      }
      cat("   ",x$cbControl$names[i],": ",mess,"\n",sep="")
    }
  }
  if (!x$save.output) 
    cat("\n","Output was not captured!", "\n")
  else {
    if (x$output[1] == "path") 
      wd <- paste(getwd(), "/", sep = "")
    else wd <- x$output[1]
if (x$output[3]=="csv") fileName <- paste(wd, x$output[2], ".csv", sep = "")
else fileName <- paste(wd, x$output[2], ".txt", sep = "")
    cat("\n","Output was captured and saved into file", "\n", 
        " '", fileName, "'", "\n", "on ",as.character(Sys.Date()),"\n","\n", sep = "")
  }
}



plot.cat<-function (x, ci = FALSE, alpha = 0.05, trueTh = TRUE, classThr = NULL, 
    save.plot = FALSE, save.options = c("path", "name", "pdf"), 
    ...) 
{
    if (!is.logical(trueTh)) 
        stop("'trueTh' must be either TRUE or FALSE", call. = FALSE)
    if (!is.logical(ci)) 
        stop("'ci' must be either TRUE or FALSE", call. = FALSE)
    if (!is.numeric(classThr) & !is.null(classThr)) 
        stop("'classThr' must be either a  numeric threshold or NULL", 
            call. = FALSE)
    internalCAT <- function() {
        res <- x
        X <- 1:length(res$testItems)
        nra <- length(res$pattern) - length(res$thetaProv)
        if (nra < 0) {
            indic <- 2:length(res$thetaProv)
            Y <- res$thetaProv[indic]
            r1 <- res$thetaProv[indic] - qnorm(1 - alpha/2) * 
                res$seProv[indic]
            r2 <- res$thetaProv[indic] + qnorm(1 - alpha/2) * 
                res$seProv[indic]
        }
        else {
            Y <- c(rep(NA, nra), res$thetaProv)
            r1 <- res$thetaProv - qnorm(1 - alpha/2) * res$seProv
            r2 <- res$thetaProv + qnorm(1 - alpha/2) * res$seProv
        }
        r1[abs(r1)==Inf]<-0
        r2[abs(r2)==Inf]<-0
        if (!is.na(res$trueTheta)) vectRange <- c(res$thetaProv, res$trueTheta)
else vectRange <- res$thetaProv
        if (ci) 
            vectRange <- c(vectRange, r1, r2)
        if (!is.null(classThr)) 
            vectRange <- c(vectRange, classThr)
        ra <- range(vectRange)
        ra[1] <- ra[1] - 0.2
        ra[2] <- ra[2] + 0.2
        if (nra >= 0) {
            r1 <- c(rep(NA, nra), r1)
            r2 <- c(rep(NA, nra), r2)
        }
        plot(X, Y, type = "o", xlab = "Item", ylab = "Ability estimate", 
            ylim = ra, cex = 0.7)
        if (ci) {
            for (i in 1:length(X)) {
                lines(rep(i, 2), c(r1[i], r2[i]), lty = 3)
                lines(c(i - 0.2, i + 0.2), rep(r1[i], 2))
                lines(c(i - 0.2, i + 0.2), rep(r2[i], 2))
            }
        }
        if (trueTh & !is.na(res$trueTheta)) 
            abline(h = res$trueTheta)
        if (!is.null(classThr)) 
            abline(h = classThr, lty = 2)
    }
    internalCAT()
    if (save.plot) {
        plotype <- NULL
        if (save.options[3] == "pdf") 
            plotype <- 1
        if (save.options[3] == "jpeg") 
            plotype <- 2
        if (is.null(plotype)) 
            cat("Invalid plot type (should be either 'pdf' or 'jpeg').", 
                "\n", "The plot was not captured!", "\n")
        else {
            if (save.options[1] == "path") 
                wd <- paste(getwd(), "/", sep = "")
            else wd <- save.options[1]
            nameFile <- paste(wd, save.options[2], switch(plotype, 
                `1` = ".pdf", `2` = ".jpg"), sep = "")
            if (plotype == 1) {
                {
                  pdf(file = nameFile)
                  internalCAT()
                }
                dev.off()
            }
            if (plotype == 2) {
                {
                  jpeg(filename = nameFile)
                  internalCAT()
                }
                dev.off()
            }
            cat("The plot was captured and saved into", "\n", 
                " '", nameFile, "'", "\n", "\n", sep = "")
        }
    }
    else cat("The plot was not captured!", "\n", sep = "")
}



