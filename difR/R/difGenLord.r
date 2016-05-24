difGenLord<-function(Data, group, focal.names, model, c=NULL, engine="ltm", discr=1, irtParam=NULL, nrFocal=2, same.scale=TRUE, anchor=NULL,alpha=0.05, purify=FALSE, nrIter=10,save.output=FALSE, output=c("out","default")) 
{
internalGLord<-function(){
if (!is.null(irtParam)){
nrItems<-nrow(irtParam)/(nrFocal+1)
dataName<-rownames(irtParam[1:nrItems,])
if (!is.null(anchor) & !same.scale){
dif.anchor<-anchor
if (is.numeric(anchor)) ANCHOR<-anchor
else{
ANCHOR<-NULL
for (i in 1:length(anchor)) ANCHOR[i]<-(1:length(dataName))[dataName==anchor[i]]
}
}
else {
ANCHOR<-1:nrItems
dif.anchor<-NULL
}
if (!same.scale){
prov<-vector("list",nrFocal+1)
for (i in 1:(nrFocal+1)) prov[[i]]<-irtParam[((i-1)*nrItems+1):(i*nrItems),]
irtParam<-prov[[1]]
for (gr in 1:nrFocal) irtParam<-rbind(irtParam,itemRescale(prov[[1]],prov[[gr+1]],items=ANCHOR))
}
mod<-as.character(ncol(irtParam))
model<-switch(mod,"2"="1PL","5"="2PL","6"="3PL","9"="3PL")
nPar<-switch(mod,"2"=1,"5"=2,"6"=2,"9"=3)
if (ncol(irtParam)!=6) Guess<-NULL
else {
Guess<-irtParam[1:nrItems,6]
if (length(unique(round(Guess,4)))==1) Guess<-unique(round(Guess,4))
}
Q<-qchisq(1-alpha,nPar*nrFocal)
itemParInit<-irtParam
estPar<-FALSE
}
else{
     if (length(group) == 1) {
           if (is.numeric(group)) {
              gr <- Data[, group]
              DATA <- Data[,(1:ncol(Data))!= group]
              colnames(DATA) <- colnames(Data)[(1:ncol(Data))!= group]
           }
           else {
              gr <- Data[, colnames(Data)==group]
              DATA <- Data[,colnames(Data)!= group]
              colnames(DATA) <- colnames(Data)[colnames(Data)!= group]
           }
    }
    else {
        gr <- group
        DATA <- Data
    }
Group <- rep(0, nrow(DATA))
nrFocal<-length(focal.names)
for (i in 1:nrFocal) Group[gr == focal.names[i]] <- i
nrItems<-ncol(DATA)
dataName<-colnames(DATA)
if (!is.null(anchor)){
dif.anchor<-anchor
if (is.numeric(anchor)) ANCHOR<-anchor
else{
ANCHOR<-NULL
for (i in 1:length(anchor)) ANCHOR[i]<-(1:ncol(DATA))[colnames(DATA)==anchor[i]]
}
}
else {
ANCHOR<-1:nrItems
dif.anchor<-NULL
}
irtParam<-NULL
GROUP<-0:nrFocal
for (indic in 1:length(GROUP)){
nr<-length(Group[Group==GROUP[indic]])
d0<-matrix(NA,nr,nrItems)
c0<-0
for (i in 1:length(Group)){
if (Group[i]==GROUP[indic]){
c0<-c0+1
d0[c0,]<-as.numeric(DATA[i,])
}
}
Guess<-c
if (is.null(Guess)) m0<-switch(model,"1PL"=itemParEst(d0,model="1PL",engine=engine,discr=discr),"2PL"=itemParEst(d0,model="2PL"),"3PL"=itemParEst(d0,model="3PL"))
else m0<-itemParEst(d0,model="3PL",c=Guess)
if (indic==1) irtParam<-m0
else irtParam<-rbind(irtParam,itemRescale(irtParam[1:nrItems,],m0,items=ANCHOR))
}
if (is.null(Guess)) nPar<-switch(model,"1PL"=1,"2PL"=2,"3PL"=3)
else nPar<-2
Q<-qchisq(1-alpha,nPar*nrFocal)
itemParInit<-irtParam
estPar<-TRUE
}
if (!purify | !is.null(anchor)) {
STATS<-genLordChi2(irtParam,nrFocal)
if ((max(STATS))<=Q) DIFitems<-"No DIF item detected"
else DIFitems<-(1:nrItems)[STATS>Q]
RES<-list(genLordChi=STATS,alpha=alpha,thr=Q,df=nPar*nrFocal,DIFitems=DIFitems,purification=purify,model=model,c=Guess,engine=engine,discr=discr,itemParInit=itemParInit,estPar=estPar,names=dataName,anchor.names=dif.anchor,focal.names=focal.names,save.output=save.output,output=output)
if (!is.null(anchor) & (RES$estPar | (!RES$estPar & !same.scale))) {
RES$genLordChi[ANCHOR]<-NA
for (i in 1:length(RES$DIFitems)){
if (sum(RES$DIFitems[i]==ANCHOR)==1) RES$DIFitems[i]<-NA
}
RES$DIFitems<-RES$DIFitems[!is.na(RES$DIFitems)]
}
}
else{
nrPur<-0
difPur<-NULL
noLoop<-FALSE
stats1<-genLordChi2(irtParam,nrFocal)
if (max(stats1)<=Q){
DIFitems<-"No DIF item detected"
noLoop<-TRUE
itemParFinal=irtParam
RES<-list(genLordChi=stats1,alpha=alpha,thr=Q,df=nPar*nrFocal,DIFitems=DIFitems,purification=purify,nrPur=nrPur,difPur=difPur,convergence=noLoop,model=model,c=Guess,engine=engine,discr=discr,itemParInit=itemParInit,itemParFinal=itemParFinal,estPar=estPar,names=dataName,anchor.names=NULL,focal.names=focal.names,save.output=save.output,output=output)
}
else{
dif<-(1:nrItems)[stats1>Q]
difPur<-rep(0,length(stats1))
difPur[dif]<-1
repeat{
if (nrPur>=nrIter) {
itemParFinal<-irtParam
break
}
else{
nrPur<-nrPur+1
nodif<-NULL
if (is.null(dif)) nodif<-1:nrItems
else{
for (i in 1:nrItems){
if (sum(i==dif)==0) nodif<-c(nodif,i)
}
}
prov<-vector("list",nrFocal+1)
for (i in 1:(nrFocal+1)) prov[[i]]<-irtParam[((i-1)*nrItems+1):(i*nrItems),]
irtParam<-prov[[1]]
for (gr in 1:nrFocal) irtParam<-rbind(irtParam,itemRescale(prov[[1]],prov[[gr+1]],items=nodif))
stats2<-genLordChi2(irtParam,nrFocal) 
if (max(stats2)<=Q) dif2<-NULL
else dif2 <- (1:nrItems)[stats2>Q]
difPur<-rbind(difPur,rep(0,nrItems))
difPur[nrPur+1,dif2]<-1
if (length(dif)!=length(dif2)) dif<-dif2
else{
dif<-sort(dif)
dif2<-sort(dif2)
if (sum(dif==dif2)==length(dif)) {
noLoop<-TRUE
itemParFinal<-irtParam
break
}
else dif<-dif2
}
}
}
if (!is.null(difPur)){
ro<-co<-NULL
for (ir in 1:nrow(difPur)) ro[ir]<-paste("Step",ir-1,sep="")
for (ic in 1:ncol(difPur)) co[ic]<-paste("Item",ic,sep="")
rownames(difPur)<-ro
colnames(difPur)<-co
}
RES<-list(genLordChi=stats2,alpha=alpha,thr=Q,df=nPar*nrFocal,DIFitems=dif2,purification=purify,nrPur=nrPur,difPur=difPur,convergence=noLoop,model=model,c=Guess,engine=engine,discr=discr,itemParInit=itemParInit,itemParFinal=itemParFinal,estPar=estPar,names=dataName,anchor.names=NULL,focal.names=focal.names,save.output=save.output,output=output)
}
}
class(RES)<-"GenLord"
return(RES)
}
resToReturn<-internalGLord()
if (save.output){
if (output[2]=="default") wd<-paste(getwd(),"/",sep="")
else wd<-output[2]
fileName<-paste(wd,output[1],".txt",sep="")
capture.output(resToReturn,file=fileName)
}
return(resToReturn)
}



# METHODS
plot.GenLord<-function (x, plot = "lordStat", item = 1, pch = 8, number = TRUE, 
    col = "red", colIC = rep("black", length(x$focal.names) + 
        1), ltyIC = 1:(length(x$focal.names) + 1), save.plot=FALSE,save.options=c("plot","default","pdf"),
    ref.name=NULL, ...) 
{
internalGLord<-function(){
    res <- x
    title <- expression(paste("Generalized Lord's ", chi^2))
    plotType <- switch(plot, lordStat = 1, itemCurve = 2)
    if (is.null(plotType)) 
        return("Error: misspecified 'type' argument")
    else {
        if (plotType == 1) {
    		if (!number) {
     			plot(res$genLordChi, xlab = "Item", ylab = expression(paste("Generalized Lord's ", 
          	 	chi^2, " statistic")), ylim = c(0, max(c(res$genLordChi, 
           		res$thr) + 1,na.rm=TRUE)), pch = pch, main = title)
       	if (!is.character(res$DIFitems)) 
           		points(res$DIFitems, res$genLordChi[res$DIFitems], 
                	pch = pch, col = col)
   	 	}
    	  	else {
       		 plot(res$genLordChi, xlab = "Item", ylab = expression(paste("Generalized Lord's ", 
          		 chi^2, " statistic")), ylim = c(0, max(c(res$genLordChi, 
          		 res$thr) + 1,na.rm=TRUE)), col = "white", main = title)
       		 text(1:length(res$genLordChi), res$genLordChi, 1:length(res$genLordChi))
       		 if (!is.character(res$DIFitems)) 
         	 	  text(res$DIFitems, res$genLordChi[res$DIFitems], 
              	  res$DIFitems, col = col)
    		}
    	  abline(h = res$thr)
	  }
	  else {
            it <- ifelse(is.character(item) | is.factor(item), 
                (1:length(res$names))[res$names == item], item)
if (is.na(res$genLordChi[it])) stop("Selected item is an anchor item!",call.=FALSE)
            J <- length(res$genLordChi)
		if (res$purification) matPar <- res$itemParFinal
            else matPar <- res$itemParInit
		nrFocal<-nrow(matPar)/J-1
		parItems<-matPar[it,]
		for (gr in 1:nrFocal) parItems<-rbind(parItems,matPar[it+gr*J,])
            nrpar <- ncol(matPar)
            nrpar <- paste("N", nrpar, sep = "")
           	parItem <- switch(nrpar, N2 = cbind(rep(1,nrFocal+1), parItems[, 1], 
                rep(0,nrFocal+1)), N5 = cbind(parItems[, 1:2], rep(0,nrFocal+1)),
		    N6 = parItems[,c(1, 2, 6)], N9 = parItems[, 1:3])
            seq <- seq(-4, 4, 0.1)
            mod <- function(t, s) t[3] + (1 - t[3]) * exp(t[1] * 
                (s - t[2]))/(1 + exp(t[1] * (s - t[2])))
            mainName <- ifelse(is.character(res$names[it]), res$names[it], 
                paste("Item ", it, sep = ""))
            plot(seq, mod(parItem[1,], seq), col = colIC[1], type = "l", 
                lty = ltyIC[1], ylim = c(0, 1), xlab = expression(theta), 
                ylab = "Probability", main = mainName)
		for (gr in 1:nrFocal) lines(seq, mod(parItem[gr+1,], seq),
			 col = colIC[gr+1], lty = ltyIC[gr+1])
            if (is.null(ref.name)) legnames <- "Reference"
            else legnames<-ref.name
                if (is.character(res$focal.names) | is.factor(res$focal.names)) 
                  legnames <- c(legnames, res$focal.names)
                else {
                  for (t in 1:length(res$focal.names)) legnames <- c(legnames, 
                    paste("Focal ", res$focal.names[t], sep = ""))
                }
                legend(-4, 1, legnames, col = colIC, lty = ltyIC, 
                  bty = "n")
          }
     }
}
internalGLord()
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
internalGLord()
}
dev.off()
}
if (plotype==2){
{
jpeg(filename=fileName)
internalGLord()
}
dev.off()
}
cat("The plot was captured and saved into","\n"," '",fileName,"'","\n","\n",sep="")
}
}
else cat("The plot was not captured!","\n",sep="")
}


print.GenLord<-function(x,...){
res<-x
cat("\n")
cat("Detection of Differential Item Functioning using generalized Lord's method","\n")
if (res$purification & is.null(res$anchor.names)) pur<-"with "
else pur<-"without "
if (is.null(res$c)) {
mod<-res$model
nrFocal<-res$df/switch(res$model,"1PL"=1,"2PL"=2,"3PL"=3)
}
else {
mod<-"constrained 3PL"
nrFocal<-res$df/2
}
cat("(",nrFocal," focal groups), with ",mod," model and ",pur, "item purification","\n","\n",sep="")
if (res$estPar){
if (res$model!="1PL" | res$engine=="ltm") cat("Engine 'ltm' for item parameter estimation","\n","\n")
else cat("Engine 'lme4' for item parameter estimation","\n","\n")
}
if (res$model=="1PL" & res$engine=="ltm") {
if (is.null(res$discr)) cat("Common discrimination parameter: estimated from 'ltm'","\n","\n")
else cat("Common discrimination parameter: fixed to ",res$discr,"\n","\n",sep="")
}
if (!is.null(res$c)){
if (length(res$c)==1) cat("Common pseudo-guessing value: ",res$c,"\n","\n",sep="")
else {
pg<-cbind(res$c)
rownames(pg)<-res$names
colnames(pg)<-"c"
cat("Common pseudo-guessing values:","\n","\n")
print(pg)
cat("\n")
}
}
if (res$purification & is.null(res$anchor.names)){
if (res$nrPur<=1) word<-" iteration"
else word<-" iterations"
if (!res$convergence) {
cat("WARNING: no item purification convergence after ",res$nrPur,word,"\n",sep="")
loop<-NULL
for (i in 1:res$nrPur) loop[i]<-sum(res$difPur[1,]==res$difPur[i+1,])
if (max(loop)!=length(res$genLordChi)) cat("(Note: no loop detected in less than ",res$nrPur,word,")","\n",sep="")
else cat("(Note: loop of length ",min((1:res$nrPur)[loop==length(res$genLordChi)])," in the item purification process)","\n",sep="")
cat("WARNING: following results based on the last iteration of the purification","\n","\n")
}
else cat("Convergence reached after ",res$nrPur,word,"\n","\n",sep="")
}
if (is.null(res$anchor.names)) {
itk<-1:length(res$genLordChi)
cat("No set of anchor items was provided", "\n", "\n")
}
else {
itk<-(1:length(res$genLordChi))[!is.na(res$genLordChi)]
cat("Anchor items (provided by the user):", "\n")
if (is.numeric(res$anchor.names)) mm<-res$names[res$anchor.names]
else mm<-res$anchor.names
mm <- cbind(mm)
rownames(mm) <- rep("", nrow(mm))
colnames(mm) <- ""
print(mm, quote = FALSE)
cat("\n", "\n")
}
cat("Generalized Lord's chi-square statistic:","\n","\n")
it<-rep("",length(res$genLordChi))
pval<-round(1-pchisq(res$genLordChi,res$df),4)
symb<-symnum(pval,c(0,0.001,0.01,0.05,0.1,1),symbols=c("***","**","*",".",""))
m1<-cbind(round(res$genLordChi[itk],4),pval[itk])
m1<-noquote(cbind(format(m1,justify="right"),symb[itk]))
if (!is.null(res$names)) rownames(m1)<-res$names[itk]
else{
rn<-NULL
for (i in 1:nrow(m1)) rn[i]<-paste("Item",i,sep="")
rownames(m1)<-rn[itk]
}
colnames(m1)<-c("Stat.","P-value","")
print(m1)
cat("\n")
cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ","\n")
cat("\n","Detection threshold: ",round(res$thr,4)," (significance level: ",res$alpha,")","\n","\n",sep="")
if (is.character(res$DIFitems)) cat("Items detected as DIF items:",res$DIFitems,"\n","\n")
else {
cat("Items detected as DIF items:","\n")
   if (!is.null(res$names)) m2 <- res$names
    else {
        rn <- NULL
        for (i in 1:length(res$genLordChi)) rn[i] <- paste("Item", i, sep = "")
        m2 <- rn
    }
        m2 <- cbind(m2[res$DIFitems])
rownames(m2)<-rep("",nrow(m2))
colnames(m2)<-""
print(m2,quote=FALSE)
cat("\n")
if (!x$save.output) cat("Output was not captured!","\n")
else {
if (x$output[2]=="default") wd<-paste(getwd(),"/",sep="")
else wd<-x$output[2]
fileName<-paste(wd,x$output[1],".txt",sep="")
cat("Output was captured and saved into file","\n"," '",fileName,"'","\n","\n",sep="")
}
}
}

