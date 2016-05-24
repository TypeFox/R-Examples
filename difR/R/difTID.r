# DIF TRANSFORMED ITEM DIFFICULTIES (ANGOFF's DELTA METHOD)

difTID<-function (Data, group, focal.name, anchor=NULL,props = NULL, thrTID = 1.5, purify = FALSE, nrIter = 10, save.output=FALSE, output=c("out","default")) 
{
internalTID<-function(){
if (is.null(props)){
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

PROPS<-matrix(NA,ncol(DATA),2)
for (i in 1:ncol(DATA)){
PROPS[i,1]<-mean(DATA[,i][Group==0],na.rm=TRUE)
PROPS[i,2]<-mean(DATA[,i][Group==1],na.rm=TRUE)
}
itNames<-colnames(DATA)
}

else 
{
PROPS<-props
itNames<-rownames(props)
}
Q <- thrTID
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
    if (!purify | !is.null(anchor)) {
        PROV <- trItemDiff(PROPS,anchor=ANCHOR)
        STATS <- PROV$dist
        if (max(abs(STATS)) <= Q) 
            DIFitems <- "No DIF item detected"
        else DIFitems <- (1:ncol(DATA))[abs(STATS) > Q]
        RES <- list(Dj = STATS, prop=PROPS,delta=PROV$delta,axisPar=PROV$pars, thr = Q, DIFitems = DIFitems, 
            purification = purify, names = itNames, anchor.names=dif.anchor,save.output=save.output,output=output)
if (!is.null(anchor)) {
RES$Dj[ANCHOR]<-NA
RES$prop[ANCHOR,]<-NA
RES$delta[ANCHOR,]<-NA
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
        prov1 <- trItemDiff(PROPS)
        stats1 <- prov1$dist
        if (max(abs(stats1)) <= Q) {
            DIFitems <- "No DIF item detected"
            noLoop <- TRUE
        }
        else {
            dif <- (1:ncol(DATA))[abs(stats1) > Q]
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
                  prov2 <- trItemDiff(PROPS, anchor = nodif)
                  stats2 <- prov2$dist
                  if (max(abs(stats2)) <= Q) 
                    dif2 <- NULL
                  else dif2 <- (1:ncol(DATA))[abs(stats2) > Q]
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
            DIFitems <- (1:ncol(DATA))[abs(stats1) > Q]
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
        RES <- list(Dj = stats1, prop=PROPS,delta=prov1$delta,axisPar=prov1$pars, thr = Q, DIFitems = DIFitems, 
            purification = purify, nrPur = nrPur, difPur = difPur, 
            convergence = noLoop, names =itNames,anchor.names=dif.anchor,save.output=save.output,output=output)
    }
    class(RES) <- "TID"
    return(RES)
}
resToReturn<-internalTID()
if (save.output){
if (output[2]=="default") wd<-paste(getwd(),"/",sep="")
else wd<-output[2]
fileName<-paste(wd,output[1],".txt",sep="")
capture.output(resToReturn,file=fileName)
}
return(resToReturn)
}







# METHODS
plot.TID<-function (x, plot="dist", pch = 8, number = TRUE, col = "red", save.plot=FALSE,save.options=c("plot","default","pdf"),...) 
{
internalTID<-function(){
    res <- x
type<-switch(plot,dist=1,delta=2)
if (is.null(type)) stop("'plot' must be either 'dist' or 'delta'!",call.=FALSE)    
if (type==1){
yl<-c(min(c(res$Dj,-abs(res$thr)),na.rm=TRUE)-0.1,max(c(res$Dj,abs(res$thr)),na.rm=TRUE)+0.1)
if (!number){
plot(res$Dj,xlab = "Item", ylab = "Perpendicular distance",ylim = yl, pch = pch, main = "Transformed Item Difficulties")
if (!is.character(res$DIFitems)) points(res$DIFitems, res$Dj[res$DIFitems], pch = pch, col = col)
}
    else {
  plot(res$Dj, xlab = "Item", ylab = "Perpendicular distance", ylim = yl, 
            col = "white", main = "Transformed Item Difficulties")
        text(1:length(res$Dj), res$Dj, 1:length(res$Dj))
        if (!is.character(res$DIFitems)) text(res$DIFitems, res$Dj[res$DIFitems], res$DIFitems,  col = col)
}
abline(h = -abs(res$thr))
abline(h = abs(res$thr))
}
else{
xl<-yl<-c(max(c(min(res$delta),1),na.rm=TRUE)-1,min(c(max(res$delta,na.rm=TRUE),25))+1)
if (!number){
plot(res$delta[,1],res$delta[,2],xlab = "Reference group", ylab = "Focal group",xlim=xl,ylim = yl, pch = pch, main = "Delta plot")
if (!is.character(res$DIFitems)) points(res$delta[res$DIFitems,1],res$delta[res$DIFitems,2], pch = pch,col = col)
}
    else {
  plot(res$delta[,1],res$delta[,2],xlab = "Reference group", ylab = "Focal group",xlim=xl,ylim = yl, col = "white", main = "Delta plot")
        text(res$delta[,1],res$delta[,2], 1:length(res$Dj))
        if (!is.character(res$DIFitems)) text(res$delta[res$DIFitems,1],res$delta[res$DIFitems,2], res$DIFitems, col = col)
}
abline(res$axisPar[1],res$axisPar[2],lty=1)
abline(res$axisPar[1]+abs(res$thr)*sqrt(res$axisPar[2]^2+1),res$axisPar[2],lty=2)
abline(res$axisPar[1]-abs(res$thr)*sqrt(res$axisPar[2]^2+1),res$axisPar[2],lty=2)
}
}
internalTID()
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
internalTID()
}
dev.off()
}
if (plotype==2){
{
jpeg(filename=fileName)
internalTID()
}
dev.off()
}
cat("The plot was captured and saved into","\n"," '",fileName,"'","\n","\n",sep="")
}
}
else cat("The plot was not captured!","\n",sep="")
}



print.TID<-function (x, ...) 
{
    res <- x
    cat("\n")
    cat("Detection of Differential Item Functioning using Transformed Item Difficulties", 
        "\n")
      if (res$purification & is.null(res$anchor.names)) pur <- "with "
    else pur <- "without "
    cat("(TID) method, ", pur, "item purification","\n", "\n", sep = "")
    if (res$purification & is.null(res$anchor.names)) {
        if (res$nrPur <= 1) 
            word <- " iteration"
        else word <- " iterations"
        if (!res$convergence) {
            cat("WARNING: no item purification convergence after ", 
                res$nrPur, word, "\n", sep = "")
            loop <- NULL
            for (i in 1:res$nrPur) loop[i] <- sum(res$difPur[1, 
                ] == res$difPur[i + 1, ])
            if (max(loop) != length(res$Dj)) 
                cat("(Note: no loop detected in less than ", 
                  res$nrPur, word, ")", "\n", sep = "")
            else cat("(Note: loop of length ", min((1:res$nrPur)[loop == 
                length(res$Dj)]), " in the item purification process)", 
                "\n", sep = "")
            cat("WARNING: following results based on the last iteration of the purification", 
                "\n", "\n")
        }
        else cat("Convergence reached after ", res$nrPur, word, "\n", "\n", sep = "")
    }
if (is.null(res$anchor.names)) {
itk<-1:length(res$Dj)
cat("No set of anchor items was provided", "\n", "\n")
}
else {
itk<-(1:length(res$Dj))[!is.na(res$Dj)]
cat("Anchor items (provided by the user):", "\n")
if (is.numeric(res$anchor.names)) mm<-res$names[res$anchor.names]
else mm<-res$anchor.names
mm <- cbind(mm)
rownames(mm) <- rep("", nrow(mm))
colnames(mm) <- ""
print(mm, quote = FALSE)
cat("\n", "\n")
}
 cat("Summary statistics:", "\n", "\n")
 pval <- round(abs(res$Dj), 4)
    symb <- symnum(pval, c(0, 0.5, 1, 1.5,Inf), symbols = c("", 
        "*", "**", "***"))
    m1 <- cbind(round(res$prop[itk,],4),round(res$delta[itk,],4),round(res$Dj[itk], 4))
    m1 <- noquote(cbind(format(m1, justify = "right"), symb[itk]))
    if (!is.null(res$names)) rownames(m1) <- res$names[itk]
    else {
        rn <- NULL
        for (i in 1:nrow(m1)) rn[i] <- paste("Item", i, sep = "")
        rownames(m1) <- rn[itk]
    }
    colnames(m1) <- c("PropRef","PropFoc","DeltaRef","DeltaFoc","Dist.", "")
    print(m1)
    cat("\n")
    cat("Class. codes: 0 '' 0.5 '*' 1.0 '**' 1.5 '***' ", 
        "\n")
    cat(" (for absolute values of 'Dist.')", "\n", "\n")
    cat("Detection threshold: ", round(res$thr, 4), "\n", "\n", sep = "")
    if (is.character(res$DIFitems)) 
        cat("Items detected as DIF items:", res$DIFitems, "\n", "\n")
    else {
        cat("Items detected as DIF items:", "\n")
   if (!is.null(res$names)) m2 <- res$names
    else {
        rn <- NULL
        for (i in 1:length(res$Dj)) rn[i] <- paste("Item", i, sep = "")
        m2 <- rn
    }
        m2 <- cbind(m2[res$DIFitems])
        rownames(m2) <- rep("", nrow(m2))
        colnames(m2) <- ""
        print(m2, quote = FALSE)
        cat("\n", "\n")
    }
    if (!res$save.output) cat("Output was not captured!","\n")
    else {
if (res$output[2]=="default") wd<-paste(getwd(),"/",sep="")
else wd<-res$output[2]
fileName<-paste(wd,res$output[1],".txt",sep="")
cat("Output was captured and saved into file","\n"," '",fileName,"'","\n","\n",sep="")
}
}






