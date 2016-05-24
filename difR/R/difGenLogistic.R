difGenLogistic<-function (Data, group, focal.names, anchor=NULL, match="score",type = "both", criterion = "LRT", 
    alpha = 0.05, purify = FALSE, nrIter = 10, save.output = FALSE, 
    output = c("out", "default")) 
{
if (purify & match[1]!="score") stop("purification not allowed when matching variable is not 'score'",call.=FALSE)
    internalGenLog <- function() {
        if (length(focal.names) == 1) 
            RES <- difLogistic(Data, group, focal.name = focal.names, member.type="group", match=match,
                type = type, alpha = alpha, purify = purify, 
                nrIter = nrIter, save.output = save.output, output = output)
        else {
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
            DF <- length(focal.names)
            for (i in 1:DF) Group[gr == focal.names[i]] <- i
            Q <- switch(type, both = qchisq(1 - alpha, 2 * DF), 
                udif = qchisq(1 - alpha, DF), nudif = qchisq(1 - 
                  alpha, DF))
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
            if (!purify| match[1]!="score" | !is.null(anchor)) {
                PROV <- genLogistik(DATA, Group, match=match,type = type, 
                  criterion = criterion,anchor=ANCHOR)
                STATS <- PROV$stat
                deltaR2 <- PROV$deltaR2
                covMat <- PROV$covMat
                if (max(STATS) <= Q) {
                  DIFitems <- "No DIF item detected"
                  logitPar <- PROV$parM1
                }
                else {
                  DIFitems <- (1:ncol(DATA))[STATS > Q]
                  logitPar <- PROV$parM1
                  for (idif in 1:length(DIFitems)) logitPar[DIFitems[idif], 
                    ] <- PROV$parM0[DIFitems[idif], ]
                }
                RES <- list(genLogistik = STATS, logitPar = logitPar, 
                  parM0 = PROV$parM0, covMat = covMat, deltaR2 = deltaR2, 
                  alpha = alpha, thr = Q, DIFitems = DIFitems, match=PROV$match,
                  type = type, purification = purify, names = colnames(DATA), 
                  anchor.names=dif.anchor,focal.names = focal.names, criterion = criterion, 
                  save.output = save.output, output = output)
if (!is.null(anchor) & match[1]=="score") {
RES$genLogistik[ANCHOR]<-NA
RES$logitPar[ANCHOR,]<-NA
RES$parM0[ANCHOR,]<-NA
RES$covMat[,,ANCHOR]<-NA
RES$deltaR2[ANCHOR]<-NA
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
                prov1 <- genLogistik(DATA, Group, match=match,type = type, 
                  criterion = criterion)
                stats1 <- prov1$stat
                deltaR2 <- prov1$deltaR2
                covMat <- prov1$covMat
                if (max(stats1) <= Q) {
                  DIFitems <- "No DIF item detected"
                  logitPar <- prov1$parM1
                  noLoop <- TRUE
                }
                else {
                  dif <- (1:ncol(DATA))[stats1 > Q]
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
                      prov2 <- genLogistik(DATA, Group, anchor = nodif, match=match,
                        type = type, criterion = criterion)
                      stats2 <- prov2$stat
                      deltaR2 <- prov2$deltaR2
                      covMat <- prov2$covMat
                      if (max(stats2) <= Q) 
                        dif2 <- NULL
                      else dif2 <- (1:ncol(DATA))[stats2 > Q]
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
                  prov1 <- prov2
                  stats1 <- stats2
                  deltaR2 <- deltaR2
                  covMat <- covMat
                  DIFitems <- (1:ncol(DATA))[stats1 > Q]
                  logitPar <- prov1$parM1
                  for (idif in 1:length(DIFitems)) logitPar[DIFitems[idif], 
                    ] <- prov1$parM0[DIFitems[idif], ]
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
                RES <- list(genLogistik = stats1, logitPar = logitPar, 
                  parM0 = prov1$parM0, covMat = covMat, deltaR2 = deltaR2, 
                  alpha = alpha, thr = Q, DIFitems = DIFitems, match=prov1$match,
                  type = type, purification = purify, nrPur = nrPur, 
                  difPur = difPur, convergence = noLoop, names = colnames(DATA), 
                  anchor.names=NULL,focal.names = focal.names, criterion = criterion, 
                  save.output = save.output, output = output)
            }
            class(RES) <- "genLogistic"
        }
        return(RES)
    }
    resToReturn <- internalGenLog()
    if (save.output) {
        if (output[2] == "default") 
            wd <- paste(getwd(), "/", sep = "")
        else wd <- output[2]
        fileName <- paste(wd, output[1], ".txt", sep = "")
        capture.output(resToReturn, file = fileName)
    }
    return(resToReturn)
}



# METHODS
plot.genLogistic <- function (x, plot = "lrStat", item = 1, itemFit="best",pch = 8, number = TRUE, col = "red", 
colIC = rep("black",length(x$focal.names)+1), 
ltyIC = 1:(length(x$focal.names)+1), title=NULL, save.plot=FALSE,save.options=c("plot","default","pdf"),ref.name=NULL,...) 
{
internalGenLog<-function(){
    res <- x
    plotType <- switch(plot, lrStat=1, itemCurve=2)
    if (is.null(plotType)) return("Error: misspecified 'type' argument")
    else {
	if (plotType==1){
    	 if (!number) {
        plot(res$genLogistik, xlab = "Item", ylab = "Generalized logistic regression statistic", 
            ylim = c(0, max(c(res$genLogistik, res$thr) + 1,na.rm=TRUE)), pch = pch, 
            main = ifelse(is.null(title),"Generalized logistic regression",title))
        if (!is.character(res$DIFitems)) 
            points(res$DIFitems, res$genLogistik[res$DIFitems], 
                pch = pch, col = col)
       }
       else {
        plot(res$genLogistik, xlab = "Item", ylab = "Generalized logistic regression statistic", 
            ylim = c(0, max(c(res$genLogistik, res$thr) + 1,na.rm=TRUE)), col = "white", 
            main = ifelse(is.null(title),"Generalized logistic regression",title))
        text(1:length(res$genLogistik), res$genLogistik, 1:length(res$genLogistik))
        if (!is.character(res$DIFitems)) 
            text(res$DIFitems, res$genLogistik[res$DIFitems], res$DIFitems, 
                col = col)
       }
       abline(h = res$thr)
	}
	else {
      it <- ifelse(is.character(item) | is.factor(item),
	           (1:length(res$names))[res$names==item],item)
if (is.na(res$logitPar[it,1])) stop("Selected item is an anchor item!",call.=FALSE)
	if (itemFit=="best") logitPar <- res$logitPar[it,]
      else logitPar<-res$parM0[it,]
	s <- seq(0,length(res$genLogistik),0.1)
 	expit <- function(t) exp(t)/(1+exp(t))
      mainName <- ifelse(is.null(title),ifelse(is.character(res$names[it]),res$names[it],paste("Item ", it, sep="")),title)
      plot(s, expit(logitPar[1]+logitPar[2]*s), col = colIC[1], type = "l",
		lty = ltyIC[1], ylim = c(0,1), xlab = "Score", ylab = "Probability",
            main=mainName)
      if (itemFit=="null" | (itemFit=="best" & !is.character(res$DIFitems) & sum(res$DIFitems==it)==1)){
	for (i in 1:length(res$focal.names))
           lines(s, expit(logitPar[1]+logitPar[2]*s+logitPar[2+i]+logitPar[2+length(res$focal.names)+i]*s),
		  col = colIC[1+i], lty = ltyIC[1+i])
      legnames <- ifelse(is.null(ref.name),"Reference",ref.name)    
      if (is.character(res$focal.names) | is.factor(res$focal.names)) legnames <- c(legnames,res$focal.names)
      else{
      for (t in 1:length(res$focal.names)) legnames <- c(legnames,paste("Focal ",
                res$focal.names[t],sep=""))
      }
	legend(0, 1, legnames, col = colIC, lty = ltyIC, bty = "n")
	}
      }
    }
}
internalGenLog()
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
internalGenLog()
}
dev.off()
}
if (plotype==2){
{
jpeg(filename=fileName)
internalGenLog()
}
dev.off()
}
cat("The plot was captured and saved into","\n"," '",fileName,"'","\n","\n",sep="")
}
}
else cat("The plot was not captured!","\n",sep="")
}



print.genLogistic<-function (x, ...) 
{
    res <- x
    cat("\n")
    mess1 <- switch(res$type, both = " both types of ", nudif = " nonuniform ", 
        udif = " uniform ")
    cat("Detection of", mess1, "Differential Item Functioning", 
        "\n", "using Generalized logistic regression method,", 
        "\n", sep = "")
    if (res$purification & is.null(res$anchor.names) & res$match=="score") 
        pur <- "with "
    else pur <- "without "
    cat(pur, "item purification ", sep = "")
    cat("and with ", length(res$focal.names), " focal groups", 
        "\n", "\n", sep = "")
    if (is.character(res$focal.names)| is.factor(res$focal.names)) {
        cat("Focal groups:", "\n")
        nagr <- cbind(res$focal.names)
        rownames(nagr) <- rep("", nrow(nagr))
        colnames(nagr) <- ""
        print(nagr, quote = FALSE)
        cat("\n")
    }
    cat("DIF flagging criterion:", ifelse(res$criterion == "Wald", 
        "Wald test", "Likelihood ratio test"), "\n", "\n")
    if (res$purification & is.null(res$anchor.names) & res$match=="score") {
        if (res$nrPur <= 1) 
            word <- " iteration"
        else word <- " iterations"
        if (!res$convergence) {
            cat("WARNING: no item purification convergence after ", 
                res$nrPur, word, "\n", sep = "")
            loop <- NULL
            for (i in 1:res$nrPur) loop[i] <- sum(res$difPur[1, 
                ] == res$difPur[i + 1, ])
            if (max(loop) != length(res$genLogistik)) 
                cat("(Note: no loop detected in less than ", 
                  res$nrPur, word, ")", "\n", sep = "")
            else cat("(Note: loop of length ", min((1:res$nrPur)[loop == 
                length(res$genLogistik)]), " in the item purification process)", 
                "\n", sep = "")
            cat("WARNING: following results based on the last iteration of the purification", 
                "\n", "\n")
        }
        else cat("Convergence reached after ", res$nrPur, word, 
            "\n", "\n", sep = "")
    }
    if (res$match=="score") cat("Matching variable: test score","\n","\n") 
    else cat("Matching variable: specified matching variable","\n","\n")
if (is.null(res$anchor.names) | res$match!="score") {
itk<-1:length(res$genLogistik)
cat("No set of anchor items was provided", "\n", "\n")
}
else {
itk<-(1:length(res$genLogistik))[!is.na(res$genLogistik)]
cat("Anchor items (provided by the user):", "\n")
if (is.numeric(res$anchor.names)) mm<-res$names[res$anchor.names]
else mm<-res$anchor.names
mm <- cbind(mm)
rownames(mm) <- rep("", nrow(mm))
colnames(mm) <- ""
print(mm, quote = FALSE)
cat("\n", "\n")
}
    cat("Generalized Logistic regression statistic:", "\n", "\n")
    nGroups <- length(res$focal.names)
    df <- switch(res$type, both = 2 * nGroups, udif = nGroups, 
        nudif = nGroups)
    pval <- round(1 - pchisq(res$genLogistik, df), 4)
    symb <- symnum(pval, c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
        "**", "*", ".", ""))
    m1 <- cbind(round(res$genLogistik[itk], 4), pval[itk])
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
    cat("\n", "Detection threshold: ", round(res$thr, 4), " (significance level: ", 
        res$alpha, ")", "\n", "\n", sep = "")
    if (is.character(res$DIFitems)) 
        cat("Items detected as DIF items:", res$DIFitems, "\n", 
            "\n")
    else {
        mess2 <- switch(res$type, both = " ", nudif = " nonuniform ", 
            udif = " uniform ")
        cat("Items detected as", mess2, "DIF items:", "\n", sep = "")
   if (!is.null(res$names)) m2 <- res$names
    else {
        rn <- NULL
        for (i in 1:length(res$genLogistik)) rn[i] <- paste("Item", i, sep = "")
        m2 <- rn
    }
        m2 <- cbind(m2[res$DIFitems])
        rownames(m2) <- rep("", nrow(m2))
        colnames(m2) <- ""
        print(m2, quote = FALSE)
        cat("\n", "\n")
    }
    cat("Effect size (Nagelkerke's R^2):", "\n", "\n")
    cat("Effect size code:", "\n")
    cat(" 'A': negligible effect", "\n")
    cat(" 'B': moderate effect", "\n")
    cat(" 'C': large effect", "\n", "\n")
    r2 <- round(res$deltaR2, 4)
    symb1 <- symnum(r2, c(0, 0.13, 0.26, 1), symbols = c("A", 
        "B", "C"))
    symb2 <- symnum(r2, c(0, 0.035, 0.07, 1), symbols = c("A", 
        "B", "C"))
    matR2 <- noquote(cbind(format(r2[itk], justify = "right"), symb1[itk], 
        symb2[itk]))
    if (!is.null(res$names)) rownames(matR2) <- res$names[itk]
    else {
        rn <- NULL
        for (i in 1:nrow(matR2)) rn[i] <- paste("Item", i, sep = "")
        rownames(matR2) <- rn[itk]
    }
    colnames(matR2) <- c("R^2", "ZT", "JG")
    print(matR2)
    cat("\n")
    cat("Effect size codes:", "\n")
    cat(" Zumbo & Thomas (ZT): 0 'A' 0.13 'B' 0.26 'C' 1", "\n")
    cat(" Jodoin & Gierl (JG): 0 'A' 0.035 'B' 0.07 'C' 1", 
        "\n")
    if (!x$save.output) 
        cat("\n", "Output was not captured!", "\n")
    else {
        if (x$output[2] == "default") 
            wd <- paste(getwd(), "/", sep = "")
        else wd <- x$output[2]
        fileName <- paste(wd, x$output[1], ".txt", sep = "")
        cat("\n", "Output was captured and saved into file", 
            "\n", " '", fileName, "'", "\n", "\n", sep = "")
    }
}

