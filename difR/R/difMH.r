# DIF MANTEL-HAENZEL

difMH<-function (Data, group, focal.name, anchor=NULL, MHstat = "MHChisq", correct = TRUE, exact=FALSE,
    alpha = 0.05, purify = FALSE, nrIter = 10, save.output = FALSE, 
    output = c("out", "default")) 
{
    internalMH <- function() {
        if (length(group) == 1) {
            if (is.numeric(group)) {
                gr <- Data[, group]
                DATA <- Data[, (1:ncol(Data)) != group]
                colnames(DATA) <- colnames(Data)[(1:ncol(Data)) != 
                  group]
            }
            else {
                gr <- Data[, colnames(Data) == group]
                DATA <- Data[, colnames(Data) != group]
                colnames(DATA) <- colnames(Data)[colnames(Data) != 
                  group]
            }
        }
        else {
            gr <- group
            DATA <- Data
        }
        Group <- rep(0, nrow(DATA))
        Group[gr == focal.name] <- 1
        Q <- switch(MHstat, MHChisq = qchisq(1 - alpha, 1), logOR = qnorm(1 - 
            alpha/2))
        if (is.null(Q)) stop("'MHstat' argument not valid", call. = FALSE)
if (!is.null(anchor)){
dif.anchor<-anchor
if (is.numeric(anchor)) ANCHOR<-anchor
else{
ANCHOR<-NULL
for (i in 1:length(anchor)) ANCHOR[i]<-(1:ncol(DATA))[colnames(DATA)==anchor[i]]
}
}
else {
ANCHOR<-1:ncol(DATA)
dif.anchor<-NULL
}
if (exact){
        if (!purify | !is.null(anchor)) {
            PROV <- mantelHaenszel(DATA, Group, correct = correct,exact=exact,anchor=ANCHOR)
                STATS <- PROV$resMH
            if (min(PROV$Pval) >=alpha) DIFitems <- "No DIF item detected"
            else DIFitems <- (1:ncol(DATA))[PROV$Pval < alpha]
            RES <- list(MH = STATS, Pval=PROV$Pval, alpha = alpha, DIFitems = DIFitems, 
                correct = correct, exact=exact, purification = purify, names = colnames(DATA), 
                anchor.names=dif.anchor,save.output = save.output, output = output)
if (!is.null(anchor)) {
RES$MH[ANCHOR]<-NA
RES$Pval[ANCHOR]<-NA
for (i in 1:length(RES$DIFitems)){
if (sum(RES$DIFitems[i]==ANCHOR)==1) RES$DIFitems[i]<-NA
}
RES$DIFitems<-RES$DIFitems[!is.na(RES$DIFitems)]
}
        }
        else {
            nrPur <- 0
            difPur <- NULL
            noLoop <- FALSE
            prov1 <- mantelHaenszel(DATA, Group, correct = correct,exact=exact)
            stats1 <- prov1$resMH
            if (min(prov1$Pval)>=alpha) {
                DIFitems <- "No DIF item detected"
                noLoop <- TRUE
            }
            else {
                dif <- (1:ncol(DATA))[prov1$Pval<alpha]
                difPur <- rep(0, length(stats1))
                difPur[dif] <- 1
                repeat {
                  if (nrPur >= nrIter) 
                    break
                  else {
                    nrPur <- nrPur + 1
                    nodif <- NULL
                    if (is.null(dif)) 
                      nodif <- 1:ncol(DATA)
                    else {
                      for (i in 1:ncol(DATA)) {
                        if (sum(i == dif) == 0) 
                          nodif <- c(nodif, i)
                      }
                    }
                    prov2 <- mantelHaenszel(DATA, Group, correct = correct, 
                      anchor = nodif,exact=exact)
                    stats2 <- prov2$resMH
                    if (min(prov2$Pval)>=alpha) dif2 <- NULL
                    else dif2 <- (1:ncol(DATA))[prov2$Pval<alpha]
                    difPur <- rbind(difPur, rep(0, ncol(DATA)))
                    difPur[nrPur + 1, dif2] <- 1
                    if (length(dif) != length(dif2)) 
                      dif <- dif2
                    else {
                      dif <- sort(dif)
                      dif2 <- sort(dif2)
                      if (sum(dif == dif2) == length(dif)) {
                        noLoop <- TRUE
                        break
                      }
                      else dif <- dif2
                    }
                  }
                }
                stats1 <- stats2
                prov1 <- prov2
                DIFitems <- (1:ncol(DATA))[prov1$Pval<alpha]
            }
            if (!is.null(difPur)) {
                ro <- co <- NULL
                for (ir in 1:nrow(difPur)) ro[ir] <- paste("Step", 
                  ir - 1, sep = "")
                for (ic in 1:ncol(difPur)) co[ic] <- paste("Item", 
                  ic, sep = "")
                rownames(difPur) <- ro
                colnames(difPur) <- co
            }
            RES <- list(MH = stats1, Pval=prov1$Pval, alpha = alpha, DIFitems = DIFitems, 
                correct = correct, exact=exact, purification = purify, nrPur = nrPur, 
                difPur = difPur, convergence = noLoop, names = colnames(DATA), 
                anchor.names=NULL, save.output = save.output, output = output)
        }
}
else{
        if (!purify | !is.null(anchor)) {
            PROV <- mantelHaenszel(DATA, Group, correct = correct,exact=exact,anchor=ANCHOR)
            if (MHstat == "MHChisq") 
                STATS <- PROV$resMH
            else STATS <- log(PROV$resAlpha)/sqrt(PROV$varLambda)
            if (max(abs(STATS),na.rm=TRUE) <= Q) 
                DIFitems <- "No DIF item detected"
            else DIFitems <- (1:ncol(DATA))[is.na(STATS)==FALSE & abs(STATS) > Q]
            RES <- list(MH = STATS, alphaMH = PROV$resAlpha, 
                varLambda = PROV$varLambda, MHstat = MHstat, 
                alpha = alpha, thr = Q, DIFitems = DIFitems, 
                correct = correct, exact=exact, purification = purify, names = colnames(DATA), 
                anchor.names=dif.anchor, save.output = save.output, output = output)
if (!is.null(anchor)) {
RES$MH[ANCHOR]<-NA
RES$alphaMH[ANCHOR]<-NA
RES$varLambda[ANCHOR]<-NA
for (i in 1:length(RES$DIFitems)){
if (sum(RES$DIFitems[i]==ANCHOR)==1) RES$DIFitems[i]<-NA
}
RES$DIFitems<-RES$DIFitems[!is.na(RES$DIFitems)]
}
        }
        else {
            nrPur <- 0
            difPur <- NULL
            noLoop <- FALSE
            prov1 <- mantelHaenszel(DATA, Group, correct = correct,exact=exact)
            if (MHstat == "MHChisq") 
                stats1 <- prov1$resMH
            else stats1 <- log(prov1$resAlpha)/sqrt(prov1$varLambda)
            if (max(abs(stats1),na.rm=TRUE) <= Q) {
                DIFitems <- "No DIF item detected"
                noLoop <- TRUE
            }
            else {
                dif <- (1:ncol(DATA))[is.na(stats1)==FALSE & abs(stats1) > Q]
                difPur <- rep(0, length(stats1))
                difPur[dif] <- 1
                repeat {
                  if (nrPur >= nrIter) 
                    break
                  else {
                    nrPur <- nrPur + 1
                    nodif <- NULL
                    if (is.null(dif) == TRUE) 
                      nodif <- 1:ncol(DATA)
                    else {
                      for (i in 1:ncol(DATA)) {
                        if (sum(i == dif) == 0) 
                          nodif <- c(nodif, i)
                      }
                    }
                    prov2 <- mantelHaenszel(DATA, Group, correct = correct, 
                      anchor = nodif,exact=exact)
                    if (MHstat == "MHChisq") 
                      stats2 <- prov2$resMH
                    else stats2 <- log(prov2$resAlpha)/sqrt(prov2$varLambda)
                    if (max(abs(stats2),na.rm=TRUE) <= Q) 
                      dif2 <- NULL
                    else dif2 <- (1:ncol(DATA))[is.na(stats2)==FALSE & abs(stats2) > 
                      Q]
                    difPur <- rbind(difPur, rep(0, ncol(DATA)))
                    difPur[nrPur + 1, dif2] <- 1
                    if (length(dif) != length(dif2)) 
                      dif <- dif2
                    else {
                      dif <- sort(dif)
                      dif2 <- sort(dif2)
                      if (sum(dif == dif2) == length(dif)) {
                        noLoop <- TRUE
                        break
                      }
                      else dif <- dif2
                    }
                  }
                }
                stats1 <- stats2
                prov1 <- prov2
                DIFitems <- (1:ncol(DATA))[is.na(stats1)==FALSE & abs(stats1) > Q]
            }
            if (is.null(difPur) == FALSE) {
                ro <- co <- NULL
                for (ir in 1:nrow(difPur)) ro[ir] <- paste("Step", 
                  ir - 1, sep = "")
                for (ic in 1:ncol(difPur)) co[ic] <- paste("Item", 
                  ic, sep = "")
                rownames(difPur) <- ro
                colnames(difPur) <- co
            }
            RES <- list(MH = stats1, alphaMH = prov1$resAlpha, 
                varLambda = prov1$varLambda, MHstat = MHstat, 
                alpha = alpha, thr = qchisq(1 - alpha, 1), DIFitems = DIFitems, 
                correct = correct, exact=exact,purification = purify, nrPur = nrPur, 
                difPur = difPur, convergence = noLoop, names = colnames(DATA), 
                anchor.names=NULL, save.output = save.output, output = output)
        }
}
        class(RES) <- "MH"
        return(RES)
    }
    resToReturn <- internalMH()
    if (save.output == TRUE) {
        if (output[2] == "default") 
            wd <- paste(getwd(), "/", sep = "")
        else wd <- output[2]
        fileName <- paste(wd, output[1], ".txt", sep = "")
        capture.output(resToReturn, file = fileName)
    }
    return(resToReturn)
}



# METHODS
plot.MH<-function (x, pch = 8, number = TRUE, col = "red", save.plot=FALSE,save.options=c("plot","default","pdf"),...) 
{
if (x$exact) stop("Error: plot is not available with exact Mantel-Haenszel test",call.=FALSE)
internalMH<-function(){
    res <- x
    if (res$MHstat == "MHChisq") 
        yl <- c(0, max(c(res$MH, res$thr) + 1,na.rm=TRUE))
    else yl <- c(min(c(res$MH, -res$thr) - 0.5,na.rm=TRUE), max(c(res$MH, 
        res$thr) + 0.5,na.rm=TRUE))
    ytitle = switch(res$MHstat, MHChisq = "MH Chi-square statistic", 
        logOR = "log OR statistic")
    if (!number) {
        plot(res$MH, xlab = "Item", ylab = ytitle, ylim = yl, 
            pch = pch, main = "Mantel-Haenszel")
        if (!is.character(res$DIFitems)) 
            points(res$DIFitems, res$MH[res$DIFitems], pch = pch, 
                col = col)
    }
    else {
        plot(res$MH, xlab = "Item", ylab = ytitle, ylim = yl, 
            col = "white", main = "Mantel-Haenszel")
        text(1:length(res$MH), res$MH, 1:length(res$MH))
        if (!is.character(res$DIFitems)) 
            text(res$DIFitems, res$MH[res$DIFitems], res$DIFitems, 
                col = col)
    }
    abline(h = res$thr)
    if (res$MHstat == "logOR") 
        abline(h = -res$thr)
}
internalMH()
if (save.plot){
plotype<-NULL
if (save.options[3]=="pdf") plotype<-1
if (save.options[3]=="jpeg") plotype<-2
if (is.null(plotype)) cat("Invalid plot type (should be either 'pdf' or 'jpeg').","\n","The plot was not captured!","\n")
else {
if (save.options[2]=="default") wd<-paste(getwd(),"/",sep="")
else wd<-save.options[2]
fileName<-paste(wd,save.options[1],switch(plotype,'1'=".pdf",'2'=".jpg"),sep="")
if (plotype==1){
{
pdf(file=fileName)
internalMH()
}
dev.off()
}
if (plotype==2){
{
jpeg(filename=fileName)
internalMH()
}
dev.off()
}
cat("The plot was captured and saved into","\n"," '",fileName,"'","\n","\n",sep="")
}
}
else cat("The plot was not captured!","\n",sep="")
}



print.MH<-function (x, ...) 
{
    res <- x
    cat("\n")
    cat("Detection of Differential Item Functioning using Mantel-Haenszel method", 
        "\n")
    if (res$correct & !res$exact) 
        corr <- "with "
    else corr <- "without "
    if (res$purification & is.null(res$anchor.names)) 
        pur <- "with "
    else pur <- "without "
    cat(corr, "continuity correction and ", pur, "item purification", 
        "\n", "\n", sep = "")
    if (res$exact) cat("Results based on exact inference","\n","\n")
    else cat("Results based on asymptotic inference","\n","\n")
    if (res$purification & is.null(res$anchor.names)) {
        if (res$nrPur <= 1) 
            word <- " iteration"
        else word <- " iterations"
        if (!res$convergence) {
            cat("WARNING: no item purification convergence after ",res$nrPur, word, "\n", sep = "")
            loop <- NULL
            for (i in 1:res$nrPur) loop[i] <- sum(res$difPur[1,] == res$difPur[i + 1, ])
            if (max(loop) != length(res$MH)) 
                cat("(Note: no loop detected in less than ", 
                  res$nrPur, word, ")", "\n", sep = "")
            else cat("(Note: loop of length ", min((1:res$nrPur)[loop == 
                length(res$MH)]), " in the item purification process)", 
                "\n", sep = "")
            cat("WARNING: following results based on the last iteration of the purification", 
                "\n", "\n")
        }
        else cat("Convergence reached after ", res$nrPur, word, 
            "\n", "\n", sep = "")
    }
if (is.null(res$anchor.names)) {
itk<-1:length(res$MH)
cat("No set of anchor items was provided", "\n", "\n")
}
else {
itk<-(1:length(res$MH))[!is.na(res$MH)]
cat("Anchor items (provided by the user):", "\n")
if (is.numeric(res$anchor.names)) mm<-res$names[res$anchor.names]
else mm<-res$anchor.names
mm <- cbind(mm)
rownames(mm) <- rep("", nrow(mm))
colnames(mm) <- ""
print(mm, quote = FALSE)
cat("\n", "\n")
}
if (res$exact) met<-"Exact statistic:"
else    met <- switch(res$MHstat, MHChisq = "Mantel-Haenszel Chi-square statistic:", 
        logOR = "Log odds-ratio statistic:")
    cat(met, "\n", "\n")
if (res$exact)  pval<-round(res$Pval,4)
else{
    if (res$MHstat == "MHChisq") pval <- round(1 - pchisq(res$MH, 1), 4)
    else pval <- round(2 * (1 - pnorm(abs(res$MH))), 4)
}
    symb <- symnum(pval, c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
        "**", "*", ".", ""))
    if (!res$exact) m1 <- cbind(round(res$MH[itk], 4), pval[itk])
else m1 <- cbind(round(res$MH[itk]), pval[itk])
    m1 <- noquote(cbind(format(m1, justify = "right"), symb[itk]))
    if (!is.null(res$names)) rownames(m1) <- res$names[itk]
    else {
        rn <- NULL
        for (i in 1:nrow(m1)) rn[i] <- paste("Item", i, sep = "")
        rownames(m1) <- rn[itk]
    }
    colnames(m1) <- c("Stat.", "P-value", "")
    print(m1)
    cat("\n")
    cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ", 
        "\n")
if (res$exact) cat("\n","Significance level: ",res$alpha, "\n", "\n", sep = "")
else cat("\n", "Detection threshold: ", round(res$thr, 4), " (significance level: ", 
        res$alpha, ")", "\n", "\n", sep = "")
    if (is.character(res$DIFitems)) cat("Items detected as DIF items:", res$DIFitems, "\n", "\n")
    else {
        cat("Items detected as DIF items:", "\n")
   if (!is.null(res$names)) m2 <- res$names
    else {
        rn <- NULL
        for (i in 1:length(res$MH)) rn[i] <- paste("Item", i, sep = "")
        m2 <- rn
    }
        m2 <- cbind(m2[res$DIFitems])
        rownames(m2) <- rep("", nrow(m2))
        colnames(m2) <- ""
        print(m2, quote = FALSE)
        cat("\n", "\n")
    }
if (!res$exact){
    cat("Effect size (ETS Delta scale):", "\n", "\n")
    cat("Effect size code:", "\n")
    cat(" 'A': negligible effect", "\n")
    cat(" 'B': moderate effect", "\n")
    cat(" 'C': large effect", "\n", "\n")
    r2 <- round(-2.35 * log(res$alphaMH), 4)
    symb1 <- symnum(abs(r2), c(0, 1, 1.5, Inf), symbols = c("A", 
        "B", "C"))
    matR2 <- cbind(round(res$alphaMH[itk], 4), r2[itk])
    matR2 <- noquote(cbind(format(matR2, justify = "right"), 
        symb1[itk]))
    if (!is.null(res$names)) rownames(matR2) <- res$names[itk]
    else {
        rn <- NULL
        for (i in 1:nrow(matR2)) rn[i] <- paste("Item", i, sep = "")
        rownames(matR2) <- rn[itk]
    }
    colnames(matR2) <- c("alphaMH", "deltaMH", "")
    print(matR2)
    cat("\n")
    cat("Effect size codes: 0 'A' 1.0 'B' 1.5 'C'", "\n")
    cat(" (for absolute values of 'deltaMH')", "\n", "\n")
}
    if (!x$save.output) cat("Output was not captured!","\n")
    else {
if (x$output[2]=="default") wd<-paste(getwd(),"/",sep="")
else wd<-x$output[2]
fileName<-paste(wd,x$output[1],".txt",sep="")
cat("Output was captured and saved into file","\n"," '",fileName,"'","\n","\n",sep="")
}
}


