simulateRespondents<-function (thetas, itemBank, responsesMatrix = NULL, model = NULL, 
    genSeed = NULL, maxItems = 50, cbControl = NULL, rmax = 1, 
    Mrmax = "restricted", start = list(fixItems = NULL, seed = NULL, 
        nrItems = 1, theta = 0, randomesque = 1, startSelect = "MFI"), 
    test = list(method = "BM", priorDist = "norm", priorPar = c(0, 
        1), range = c(-4, 4), D = 1, parInt = c(-4, 4, 33), itemSelect = "MFI", 
        infoType = "observed", randomesque = 1, AP = 1, constantPatt = NULL), 
    stop = list(rule = "length", thr = 20, alpha = 0.05), final = list(method = "BM", 
        priorDist = "norm", priorPar = c(0, 1), range = c(-4, 
            4), D = 1, parInt = c(-4, 4, 33), alpha = 0.05), 
    save.output = FALSE, output = c("", "catR", "csv")) 
{
    if (length(thetas) == 1) {
        if (!is.null(responsesMatrix)) {
            resp <- as.matrix(responsesMatrix)
            res <- randomCAT(trueTheta = thetas, itemBank = itemBank, 
                responses = resp[1, ], model = model, maxItems = maxItems, 
                cbControl = cbControl, start = start, test = test, 
                stop = stop, final = final, save.output = save.output, 
                output = output, allTheta = TRUE)
        }
        else {
            res <- randomCAT(trueTheta = thetas, itemBank = itemBank, 
                model = model, maxItems = maxItems, cbControl = cbControl, 
                start = start, test = test, stop = stop, final = final, 
                save.output = save.output, output = output, allTheta = TRUE)
        }
        return(res)
    }
    else {
        if (!is.null(genSeed)) {
            if (length(genSeed) != length(thetas)) 
                stop("'thetas' and 'genSeed' must have the same length!", 
                  call. = FALSE)
        }
        internalSR <- function() {
            start.time <- Sys.time()
            respondents <- length(thetas)
            nAvailable <- NULL
            if (respondents < 1) 
                stop(paste("Length of 'thetas' has and invalid value: ", 
                  length(respondents), sep = ""))
            bank_size <- nrow(itemBank)
            estimatedThetas <- NULL
            vItemExposure <- NULL
            kpar <- rep(1, bank_size)
            last_shown <- -1
            exposure <- rep(0, bank_size)
            exposureRates <- rep(0, bank_size)
            numberItems <- NULL
            totalSeFinal <- c()
            thrOK <- NULL
            if (!is.null(responsesMatrix)) 
                resp <- as.matrix(responsesMatrix)
            if (stop$rule == "length") 
                itemsRow = stop$thr
            if (stop$rule == "precision" | stop$rule == "classification") 
                itemsRow = maxItems
            row.head1 <- rep("items.administrated", itemsRow)
            row.head1[1:length(row.head1)] <- paste(row.head1[1:length(row.head1)], 
                1:length(row.head1), sep = ".")
            row.head2 <- rep("responses", itemsRow)
            row.head2[1:length(row.head2)] <- paste(row.head2[1:length(row.head2)], 
                1:length(row.head2), sep = ".")
            row.head3 <- rep("estimated.theta", itemsRow)
            row.head3[1:length(row.head3)] <- paste(row.head3[1:length(row.head3)], 
                1:length(row.head3), sep = ".")
            row.head <- c()
            row.head <- c(row.head, "respondent", "true.theta", 
                row.head1, row.head2, "start.theta", row.head3)
            responses.df <- data.frame()
            for (i in 1:respondents) {
                if (!floor((i - 1) * 10/respondents) == last_shown) {
                  last_shown <- floor((i - 1) * 10/respondents)
                  cat("Simulation process: ", last_shown * 10, 
                    "%\n")
                }
                if (rmax < 1) {
                  nAvailable <- rep(1, bank_size)
                  if (Mrmax == "restricted") 
                    nAvailable[exposureRates >= rmax] <- 0
                  if (Mrmax == "IE") {
                    kpar[(exposureRates/kpar) <= rmax] <- 1
                    kpar[(exposureRates/kpar) > rmax] <- rmax * 
                      kpar[(exposureRates/kpar) > rmax]/exposureRates[(exposureRates/kpar) > 
                      rmax]
                    nAvailable[runif(bank_size) > kpar] <- 0
                  }
                }
                if (!is.null(responsesMatrix)) {
                  rCAT <- randomCAT(trueTheta = thetas[i], itemBank = itemBank, 
                    responses = resp[i, ], model = model, genSeed = genSeed[i], 
                    maxItems = maxItems, cbControl = cbControl, 
                    nAvailable = nAvailable, start = start, test = test, 
                    stop = stop, final = final, allTheta = TRUE)
                }
                else {
                  rCAT <- randomCAT(trueTheta = thetas[i], itemBank = itemBank, 
                    model = model, genSeed = genSeed[i], maxItems = maxItems, 
                    cbControl = cbControl, nAvailable = nAvailable, 
                    start = start, test = test, stop = stop, 
                    final = final, allTheta = TRUE)
                }
                estimatedThetas <- c(estimatedThetas, rCAT$thFinal)
                vItemExposure <- c(vItemExposure, rCAT$testItems)
                exposure[rCAT$testItems[1:length(rCAT$testItems)]] <- exposure[rCAT$testItems[1:length(rCAT$testItems)]] + 
                  1
                exposureRates = exposure/i
                numberItems <- c(numberItems, length(rCAT$testItems))
                totalSeFinal <- c(totalSeFinal, rCAT$seFinal)
                if (rCAT$stopRule == "precision") 
                  thrOK <- c(thrOK, rCAT$seFinal <= stop$thr)
                if (rCAT$stopRule == "classification") 
                  thrOK <- c(thrOK, (rCAT$trueTheta - stop$thr) * 
                    (rCAT$thFinal - stop$thr) >= 0)
                if (rCAT$stopRule == "length") 
                  thrOK <- c(thrOK, 1)
                items.administrated <- rep(-99, itemsRow)
                responses <- rep(-99, itemsRow)
                provisional.theta <- rep(-99, itemsRow)
                items.administrated[1:length(rCAT$testItems)] <- rCAT$testItems
                responses[1:length(rCAT$pattern)] <- rCAT$pattern
                if (rCAT$startNrItems == 0) 
                  rCAT$thetaProv <- rCAT$thetaProv[2:length(rCAT$thetaProv)]
                provisional.theta[1:length(rCAT$thetaProv)] <- rCAT$thetaProv
                row <- c()
                row <- c(row, i, rCAT$trueTheta, items.administrated, 
                  responses, rCAT$startTheta, provisional.theta)
                responses.df <- rbind(responses.df, row)
            }
            colnames(responses.df) <- row.head
            final.values.df <- data.frame(thetas, estimatedThetas, 
                totalSeFinal, numberItems)
            colnames(final.values.df) <- c("true.theta", "estimated.theta", 
                "final.SE", "total.items.administrated")
            resCor <- cor(thetas, estimatedThetas)
            RMSE <- sqrt(sum((estimatedThetas - thetas)^2)/respondents)
            bias <- sum(estimatedThetas - thetas)/respondents
            testLength = sum(exposureRates)
            position_min <- which(exposureRates == min(exposureRates))
            position_max <- which(exposureRates == max(exposureRates))
            overlap <- sum(exposureRates^2)/testLength
            condTheta <- rep(0, 10)
            condRMSE <- rep(0, 10)
            condBias <- rep(0, 10)
            condnItems <- rep(0, 10)
            condSE <- rep(0, 10)
            condthrOK <- rep(0, 10)
            ndecile <- rep(0, 10)
            for (z in 1:10) {
                if (z < 10) 
                  subset <- which(findInterval(thetas, quantile(thetas, 
                    seq(0, 1, 0.1))) == z)
                else subset <- c(which(findInterval(thetas, quantile(thetas, 
                  seq(0, 1, 0.1))) == z), which(thetas == max(thetas)))
                condTheta[z] <- mean(thetas[subset])
                condRMSE[z] <- sqrt(sum((estimatedThetas[subset] - 
                  thetas[subset])^2)/length(subset))
                condBias[z] <- sum(estimatedThetas[subset] - 
                  thetas[subset])/length(subset)
                condnItems[z] <- mean(numberItems[subset])
                condSE[z] <- mean(totalSeFinal[subset])
                condthrOK[z] <- mean(thrOK[subset])
                ndecile[z] <- length(subset)
            }
            finish.time <- Sys.time()
            cat("Simulation process: ", 100, "%\n")
            res <- list(thetas = thetas, itemBank = itemBank, 
                responsesMatrix = responsesMatrix, model = model, 
                genSeed = genSeed, maxItems = maxItems, cbControl = cbControl, 
                rmax = rmax, Mrmax = Mrmax, start = start, test = test, 
                stop = stop, final = final, save.output = save.output, 
                output = output, estimatedThetas = estimatedThetas, 
                correlation = resCor, bias = bias, RMSE = RMSE, 
                thrOK = thrOK, exposureRates = exposureRates, 
                testLength = testLength, overlap = overlap, numberItems = numberItems, 
                condTheta = condTheta, condBias = condBias, condRMSE = condRMSE, 
                condnItems = condnItems, condSE = condSE, condthrOK = condthrOK, 
                ndecile = ndecile, final.values.df = final.values.df, 
                responses.df = responses.df, start.time = start.time, 
                finish.time = finish.time)
            class(res) <- "catResult"
            return(res)
        }
        resToReturn <- internalSR()
        if (save.output) {
            if (output[2] != "") 
                output[2] <- paste0(output[2], ".")
            if (output[1] == "") 
                wd <- paste(getwd(), "/", sep = "")
            else wd <- output[1]
            fileName1 <- paste(wd, output[2], "main.", output[3], 
                sep = "")
            fileName2 <- paste(wd, output[2], "responses.", output[3], 
                sep = "")
            fileName3 <- paste(wd, output[2], "tables.", output[3], 
                sep = "")
            fileName4 <- paste(wd, output[2], "deciles.", output[3], 
                sep = "")
            capture.output(resToReturn, file = fileName1)
            if (output[3] == "csv") 
                sep <- ";"
            else sep <- "\t"
            write.table(resToReturn$responses.df, fileName2, 
                quote = FALSE, sep = sep, row.names = FALSE)
            write.table(resToReturn$final.values.df, file = fileName3, 
                sep = sep, quote = FALSE, row.names = FALSE)
        }
        return(resToReturn)
    }
}




print.catResult <- function(x, ...) {
  if (is.null(x$responsesMatrix))
    simulation <- FALSE else
    simulation <- TRUE
  if (!simulation) {
    cat("\n","** Simulation of multiple examinees **", "\n", "\n")
if (is.null(x$genSeed)) cat("Random seed was not fixed", "\n", "\n")
else cat("Random seed was fixed (see argument 'genSeed')", "\n", "\n")
}
  else
    cat("\n","** Post-hoc simulation of multiple examinees **", "\n", "\n")
  if (difftime(x$finish.time, x$start.time, units="hours") > 1) {
    dif_time <- difftime(x$finish.time, x$start.time, units="hours")    
    units <- "hours"
  }
  else if (difftime(x$finish.time, x$start.time, units="mins") > 1) {
    dif_time <- difftime(x$finish.time, x$start.time, units="mins")    
    units <- "minutes"
  }
  else {
    dif_time <- difftime(x$finish.time, x$start.time, units="secs")    
    units <- "seconds"
  }
  cat("Simulation time:", round(dif_time, digits=4), units ,  "\n", "\n")
  cat("Number of simulees:", length(x$thetas), "\n")
  cat("Item bank size:",length(x$itemBank[,1]), "items", "\n")
  if (is.null(x$model))
    cat("IRT model:","4PL", "\n", "\n")
  else
    cat("IRT model:",x$model, "\n", "\n")
  cat("Item selection criterion:", x$test$itemSelect,"\n")
  #cat("Starting selection rule:", x$start$startSelect  , "\n")
  cat("Stopping rule:", x$stop$rule, "\n")
  if (x$stop$rule=="length")
    cat("\t", "Test length:", x$stop$thr, "\n")
  if (x$stop$rule=="precision")
    cat("\t", "Standard error:", x$stop$thr, "\n")
  if (x$stop$rule=="classification")
    cat("\t", "Ability level:", x$stop$thr, "\n")
  cat("rmax:",  x$rmax, "\n")
  if (x$rmax < 1)
    cat("\t", "Restriction method:",  x$Mrmax, "\n")
  cat("\n")
  cat("Mean test length:", x$testLength, "items", "\n")
    if(!simulation)
    cat("Correlation(true thetas,estimated thetas):", round(x$correlation, 4), "\n")
  else
    cat("Correlation(assigned thetas,CAT estimated thetas):", round(x$correlation, 4) ,"\n")
  cat("RMSE:", round(x$RMSE,4), "\n")
  cat("Bias:", round(x$bias,4), "\n")
  if (x$stop$rule != "length")
    cat("Proportion of simulees that satisfy the stop criterion:", mean(x$thrOK), "\n", "\n")
  cat("Maximum exposure rate:", round(max(x$exposureRates),4),"\n")
  cat("Number of item(s) with maximum exposure rate:", length(which(x$exposureRates==max(x$exposureRates))), "\n")
  cat("Minimum exposure rate:", round(min(x$exposureRates),4),"\n")
  cat("Number of item(s) with minimum exposure rate:", length(which(x$exposureRates==min(x$exposureRates))), "\n")
  cat("Item overlap rate:", round(x$overlap,4), "\n","\n")
  cat("Conditional results", "\n")
  
  condDeciles <- data.frame()
  decTheta <- c("Mean Theta", round(x$condTheta,3))
  decRMSE <- c("RMSE", round(x$condRMSE,3))
  deccondBias <- c("Mean bias", round(x$condBias,3))
  deccondnItems <- c("Mean test length", round(x$condnItems,3))
  deccondSE <- c("Mean standard error", round(x$condSE,3))
  if (x$stop$rule == "precision") 
    deccondthrOK <- c("Proportion stop satisfying SE", round(x$condthrOK,3))
  else
    deccondthrOK <- c("Proportion correct classification", round(x$condthrOK,3))
  decndecile <- c("Number of simulees", x$ndecile)
  if (x$stop$rule == "length") 
    condDeciles <- rbind(condDeciles, decTheta, decRMSE, deccondBias, decndecile)
  if (x$stop$rule == "precision") 
    condDeciles <- rbind(condDeciles, decTheta, decRMSE, deccondBias, deccondnItems, deccondSE, deccondthrOK, decndecile)
  if (x$stop$rule == "classification") 
    condDeciles <- rbind(condDeciles, decTheta, decRMSE, deccondBias, deccondnItems, deccondSE, deccondthrOK, decndecile)
  colnames(condDeciles) <- c("Measure","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10")  
  print(condDeciles, row.names=FALSE)
  cat("\n")
  if (x$save.output) {
    if (x$output[2] != "")
      x$output[2] <- paste0(x$output[2],".")
    if (x$output[1]=="")
      wd <- paste(getwd(), "/",sep = "")
    else
      wd <- x$output[1]
    fileName1 <- paste(wd, x$output[2], "main.", x$output[3], sep = "")
    fileName2 <- paste(wd, x$output[2], "responses.", x$output[3], sep = "")
    fileName3 <- paste(wd, x$output[2], "tables.", x$output[3], sep = "")
    cat("These results were saved in files:","\n",fileName1,"\n",fileName2,"\n",fileName3,"\n")
  }
  else cat("These results can be saved by setting 'save.output' to TRUE","\n"," in the 'simulateRespondents' function","\n")  
}




plot.catResult <- function(x, type="all", deciles="theta", save.plot=FALSE, save.options=c("","plot","pdf"), res=300, ...) {
  type <- switch(type, all="all",trueEst="trueEst",expRate="expRate",cumExpRate="cumExpRate",
            cumNumberItems="cumNumberItems",expRatePara="expRatePara",condBias="condBias",condRMSE="condRMSE",
            numberItems="numberItems",sError="sError",condThr="condThr")
  if (is.null(x$responsesMatrix))
    simulation <- FALSE else
    simulation <- TRUE
  if (length(x$model)==0)
    x$model <- "dicho"
  if (deciles!="theta" & deciles!="deciles")
    stop("'deciles' must be either 'theta' or 'deciles'",call.=FALSE)
  if (is.null(type))
    stop("invalid 'type' argument",call.=FALSE)
  if ((type=="cumNumberItems" | type=="condThr" | type=="numberItems") & x$stop$rule=="length")
    stop("mismatch between 'type' value and 'length' stopping rule",call.=FALSE)
  if (type=="expRatePara" & (x$model=="PCM" | x$model=="NRM"))
    stop("mismatch between 'expRatePara' type value and 'PCM' or 'NRM' model",call.=FALSE)
  internalCAT<-function() {
    if (deciles == "theta")
      xline <- x$condTheta
    else
      xline <- 1:10
    if(simulation == TRUE) {
      xAccuracy <- "Assigned Theta"
      yAccuracy <- "CAT Estimated Theta"
      if (deciles == "deciles")
        thetasDeciles <- "Assigned Theta Deciles"
      else
        thetasDeciles <- "Full Bank Estimated Theta"
    } else {
      xAccuracy <- "True Theta"
      yAccuracy <- "Estimated Theta"
      if (deciles == "deciles")
        thetasDeciles <- "True Theta Deciles"
      else
        thetasDeciles <- "True Theta"
    }
    plot.trueEst <- function(x, ...) {
      plot(x$thetas, x$estimatedThetas, main="Accuracy", xlab=xAccuracy, ylab=yAccuracy)
      abline(lm(x$estimatedThetas ~ x$thetas), col="red")
    }
    plot.expRate <- function(x, ...) {
      plot(sort(x$exposureRates, decreasing = TRUE), type="l", main="Exposure Rates", xlab="Item Rank", ylab="Item Exposure Rate")
    }
    plot.cumExpRate <- function(x, ...) {
      plot(cumsum(sort(x$exposureRates,decreasing = TRUE))/x$testLength, type="l", main="Cumulative Exposure Rates", xlab="Item Rank", ylab="Cumulative item Exposure Rate")
    }
    plot.cumNumberItems <- function(x, ...) {
      respondents<-length(x$numberItems)
      plot(seq(1,respondents,1)*100/respondents, sort(x$numberItems,decreasing = FALSE), type="l", main="Test Length", xlab="Percentage of Examinees", ylab="Test Length")
    }
    plot.expRatePara <- function(x, ...) {
      plot(x$itemBank[,1], x$exposureRates, main="Exposure and a Parameter", xlab="Discrimination Parameter", pch=20, ylab="Item Exposure Rate")
    }
    plot.condBias <- function(x, ...) {
      plot(xline, x$condBias, type="o", main="Conditional Bias", xlab=thetasDeciles, ylab="Bias")
    }
    plot.condRMSE <- function(x, ...) {
      plot(xline, x$condRMSE, type="o", main="Conditional RMSE", xlab=thetasDeciles, ylab="RMSE")
    }
    plot.numberItems <-function(x, ...){
      plot(xline, x$condnItems, type="o", main="Conditional Test Length", xlab=thetasDeciles, ylab="Test Length")
    }
    plot.sError <-function(x, ...){
      plot(xline, x$condSE, type="o", main="Conditional Standard Error", xlab=thetasDeciles, ylab="Standard Error")
    }
    plot.condThr <-function(x, ...){
      if (x$stop$rule=="precision")
        plot(xline, x$condthrOK, type="o", main="Stop Satisfying SE", xlab=thetasDeciles, ylab="Proportion")
      else
        plot(xline, x$condthrOK, type="o", main="Correct Classification", xlab=thetasDeciles, ylab="Proportion")
    }
    if (type=="all"){
      if (x$stop$rule == "precision" | x$stop$rule == "classification") {  
        par(mfrow=c(3,3))
        plot.trueEst(x) 
        plot.condBias(x)
        plot.condRMSE(x)
        plot.expRate(x) 
        plot.condThr(x)
        plot.cumNumberItems(x)
        if (x$model!="PCM" & x$model!="NRM")
          plot.expRatePara(x)
        plot.numberItems(x)
        plot.sError(x)
        par(mfrow=c(1,1))
      }
      else {
        par(mfrow=c(2,3))
        plot.trueEst(x) 
        plot.condBias(x)
        plot.condRMSE(x)
        plot.expRate(x) 
        plot.cumExpRate(x)
        if (x$model!="PCM" | x$model!="NRM")
          plot.expRatePara(x)
        par(mfrow=c(1,1))
        }
    }
    if (type=="trueEst")
      plot.trueEst(x) 
    if (type=="expRate")
      plot.expRate(x) 
    if (type=="cumExpRate")
      plot.cumExpRate(x) 
    if (type=="cumNumberItems")
      plot.cumNumberItems(x) 
    if (type=="expRatePara")
      plot.expRatePara(x) 
    if (type=="condBias")
      plot.condBias(x) 
    if (type=="condRMSE")
      plot.condRMSE(x) 
    if (type=="numberItems")
      plot.numberItems(x) 
    if (type=="sError")
      plot.sError(x) 
    if (type=="condThr") 
      plot.condThr(x) 
    
  }  
  internalCAT()
  if (save.plot) {
    plotype <- NULL
  if (save.options[3] == "pdf") 
    plotype <- 1
  if (save.options[3] == "jpeg") 
    plotype <- 2
  if (is.null(plotype)) 
    cat("Invalid plot type (should be either 'pdf' or 'jpeg').", "\n", "The plot was not captured!", "\n")
  else {
    if (save.options[1] == "") 
      wd <- paste(getwd(), "/", sep = "")
    else wd <- save.options[1]
    nameFile <- paste(wd, save.options[2], switch(plotype, `1` = ".pdf", `2` = ".jpg"), sep = "")
    if (plotype == 1) {
      if (type=="all" & x$stop$rule == "precision") {
        {
        pdf(file = nameFile,width=15,height=15)
        internalCAT()
        }
        dev.off()
      }
      if (type=="all" & x$stop$rule == "classification") {
        {
        pdf(file = nameFile,width=15,height=15)
        internalCAT()
        }
        dev.off()
      }
      if (type=="all" & x$stop$rule == "length"){
        {
        pdf(file = nameFile,width=10,height=5)
        internalCAT()
        }
        dev.off()
      }
      if (type!="all"){
        {
        pdf(file = nameFile)
        internalCAT()
        }
        dev.off()
      }
    }
    if (plotype == 2) {
      if (type=="all" & x$stop$rule == "precision"){
        {
        jpeg(filename = nameFile,width=24,height=24,units="cm",res=res)
        internalCAT()
        }
        dev.off()
      }
      if (type=="all" & x$stop$rule == "classification"){
        {
        jpeg(filename = nameFile,width=24,height=24,units="cm",res=res)
        internalCAT()
        }
        dev.off()
      }
      if (type=="all" & x$stop$rule == "length"){
        {
        jpeg(filename = nameFile,width=24,height=16,units="cm",res=res)
        internalCAT()
        }
        dev.off()
      }
      if (type!="all"){
        {
        jpeg(filename = nameFile)
        internalCAT()
        }
        dev.off()
      }
    }
    cat("The plot was captured and saved into", "\n", " '", nameFile, "'", "\n", "\n", sep = "")
    }
  }
  else cat("The plot was not captured!", "\n", sep = "")
}

