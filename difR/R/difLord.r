difLord<-function(Data, group, focal.name, model, c=NULL, engine="ltm", discr=1,irtParam=NULL, same.scale=TRUE, anchor=NULL, alpha=0.05, purify=FALSE, nrIter=10, save.output=FALSE, output=c("out","default"))
{
internalLord<-function(){
if (!is.null(irtParam)){
nrItems<-nrow(irtParam)/2
m0<-irtParam[1:nrItems,]
m1<-irtParam[(nrItems+1):(2*nrItems),]
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
if (same.scale) m1p<-m1
else m1p<-itemRescale(m0,m1,items=ANCHOR)
mod<-as.character(ncol(irtParam))
model<-switch(mod,"2"="1PL","5"="2PL","6"="3PL","9"="3PL")
if (ncol(irtParam)!=6) Guess<-NULL
else {
Guess<-irtParam[1:nrItems,6]
if (length(unique(round(Guess,4)))==1) Guess<-unique(round(Guess,4))
}
if (is.null(Guess)) Q<-switch(model,"1PL"=qchisq(1-alpha,1),"2PL"=qchisq(1-alpha,2),"3PL"=qchisq(1-alpha,3))
else Q<-qchisq(1-alpha,2)
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
    Group[gr == focal.name] <- 1

nr1<-sum(Group)
nr0<-length(Group)-nr1

d0<-matrix(NA,nr0,ncol(DATA))
d1<-matrix(NA,nr1,ncol(DATA))

c0<-c1<-0
for (i in 1:length(Group)){
if (Group[i]==0){
c0<-c0+1
d0[c0,]<-as.numeric(DATA[i,])
}
else {
c1<-c1+1
d1[c1,]<-as.numeric(DATA[i,])
}
}
Guess<-c
if (is.null(Guess)) {
Q<-switch(model,"1PL"=qchisq(1-alpha,1),"2PL"=qchisq(1-alpha,2),"3PL"=qchisq(1-alpha,3))
m0<-switch(model,"1PL"=itemParEst(d0,model="1PL",engine=engine,discr=discr),"2PL"=itemParEst(d0,model="2PL"),"3PL"=itemParEst(d0,model="3PL"))
m1<-switch(model,"1PL"=itemParEst(d1,model="1PL",engine=engine,discr=discr),"2PL"=itemParEst(d1,model="2PL"),"3PL"=itemParEst(d1,model="3PL"))
}
else{
Q<-qchisq(1-alpha,2)
m0<-itemParEst(d0,model="3PL",c=Guess)
m1<-itemParEst(d1,model="3PL",c=Guess)
}
nrItems<-ncol(DATA)
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
m1p<-itemRescale(m0,m1,items=ANCHOR)
irtParam<-rbind(m0,m1p)
same.scale<-TRUE
dataName<-colnames(DATA)
itemParInit<-rbind(m0,m1)
estPar<-TRUE
}
if (!purify | !is.null(anchor)) {
STATS<-LordChi2(m0,m1p)
if ((max(STATS))<=Q) DIFitems<-"No DIF item detected"
else DIFitems<-(1:nrItems)[STATS>Q]
RES<-list(LordChi=STATS,alpha=alpha,thr=Q,DIFitems=DIFitems,purification=purify,model=model,c=Guess,engine=engine,discr=discr,itemParInit=itemParInit,estPar=estPar,names=dataName,anchor.names=dif.anchor,save.output=save.output,output=output)
if (!is.null(anchor) & (RES$estPar | (!RES$estPar & !same.scale))) {
RES$LordChi[ANCHOR]<-NA
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
stats1<-LordChi2(m0,m1p)
if (max(stats1)<=Q){
DIFitems<-"No DIF item detected"
noLoop<-TRUE
itemParFinal=rbind(m0,m1p)
RES<-list(LordChi=stats1,alpha=alpha,thr=Q,DIFitems=DIFitems,purification=purify,nrPur=nrPur,difPur=difPur,convergence=noLoop,model=model,c=Guess,engine=engine,discr=discr,itemParInit=itemParInit,itemParFinal=itemParFinal,estPar=estPar,names=dataName,anchor.names=NULL, save.output=save.output,output=output)
}
else{
dif<-(1:nrItems)[stats1>Q]
difPur<-rep(0,length(stats1))
difPur[dif]<-1
repeat{
if (nrPur>=nrIter) {
itemParFinal<-rbind(m0,itemRescale(m0,m1,items=nodif))
break
}
else{
nrPur<-nrPur+1
nodif<-NULL
if (is.null(dif)==TRUE) nodif<-1:nrItems
else{
for (i in 1:nrItems){
if (sum(i==dif)==0) nodif<-c(nodif,i)
}
}
stats2<-LordChi2(m0,itemRescale(m0,m1,items=nodif))
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
itemParFinal<-rbind(m0,itemRescale(m0,m1,items=nodif))
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
RES<-list(LordChi=stats2,alpha=alpha,thr=Q,DIFitems=dif2,purification=purify,nrPur=nrPur,difPur=difPur,convergence=noLoop,model=model,c=Guess,engine=engine,discr=discr,itemParInit=itemParInit,itemParFinal=itemParFinal,estPar=estPar,names=dataName,anchor.names=NULL, save.output=save.output,output=output)
}
}
class(RES)<-"Lord"
return(RES)
}
resToReturn<-internalLord()
if (save.output){
if (output[2]=="default") wd<-paste(getwd(),"/",sep="")
else wd<-output[2]
fileName<-paste(wd,output[1],".txt",sep="")
capture.output(resToReturn,file=fileName)
}
return(resToReturn)
}


# METHODS

plot.Lord<-function (x, plot = "lordStat", item = 1, pch = 8, number = TRUE, 
    col = "red", colIC = rep("black", 2), ltyIC = c(1, 2), save.plot=FALSE,save.options=c("plot","default","pdf"), group.names=NULL, ...) 
{
internalLord<-function(){
    res <- x
    title <- expression(paste("Lord's ", chi^2))
    plotType <- switch(plot, lordStat = 1, itemCurve = 2)
    if (is.null(plotType)) 
        return("Error: misspecified 'type' argument")
    else {
       if (plotType == 1) {
    		if (!number) {
       	 plot(res$LordChi, xlab = "Item", ylab = expression(paste(chi^2, 
         	   " statistic")), ylim = c(0, max(c(res$LordChi, res$thr) + 
       	     1,na.rm=TRUE)), pch = pch, main = title)
      	  if (!is.character(res$DIFitems)) 
       	     points(res$DIFitems, res$LordChi[res$DIFitems], pch = pch, 
           	     col = col)
   	 	}
   		else {
     	   	plot(res$LordChi, xlab = "Item", ylab = expression(paste(chi^2, 
      	      " statistic")), ylim = c(0, max(c(res$LordChi, res$thr) + 
      	      1,na.rm=TRUE)), col = "white", main = title)
      	  text(1:length(res$LordChi), res$LordChi, 1:length(res$LordChi))
      	  if (!is.character(res$DIFitems)) 
      	      text(res$DIFitems, res$LordChi[res$DIFitems], res$DIFitems, 
       	         col = col)
       	}
   	 abline(h = res$thr)
	}
	else{
            it <- ifelse(is.character(item) | is.factor(item), 
                (1:length(res$names))[res$names == item], item)
if (is.na(res$LordChi[it])) stop("Selected item is an anchor item!",call.=FALSE)
		J<-length(res$LordChi)
		if (res$purification) matPar<-res$itemParFinal
		else matPar<-rbind(res$itemParInit[1:J,],
				   itemRescale(res$itemParInit[1:J,],
						res$itemParInit[(J+1):(2*J),]))
		nrpar<-ncol(matPar)
		nrpar<-paste("N",nrpar,sep="")
		parRef<-switch(nrpar,N2=c(1,matPar[it,1],0),
					   N5=c(matPar[it,1:2],0),
					   N6=matPar[it,c(1,2,6)],
					   N9=matPar[it,1:3])
		parFoc<-switch(nrpar,N2=c(1,matPar[J+it,1],0),
					   N5=c(matPar[J+it,1:2],0),
					   N6=matPar[J+it,c(1,2,6)],
					   N9=matPar[J+it,1:3])
            seq <- seq(-4,4, 0.1)
            mod <- function(t,s) t[3]+(1-t[3])*exp(t[1]*(s-t[2]))/(1+exp(t[1]*(s-t[2])))
            mainName <- ifelse(is.character(res$names[it]), res$names[it], 
                paste("Item ", it, sep = ""))
            plot(seq, mod(parRef,seq), col = colIC[1], 
                type = "l", lty = ltyIC[1], ylim = c(0, 1), xlab = expression(theta), 
                ylab = "Probability", main = mainName)
            lines(seq, mod(parFoc,seq), col = colIC[2], 
                  lty = ltyIC[2])
            if (is.null(group.names)) legnames<-c("Reference", "Focal")
            else legnames<-group.names
            legend(-4, 1, legnames, col = colIC, 
                  lty = ltyIC, bty = "n")
      }
   }
}
internalLord()
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
internalLord()
}
dev.off()
}
if (plotype==2){
{
jpeg(filename=fileName)
internalLord()
}
dev.off()
}
cat("The plot was captured and saved into","\n"," '",fileName,"'","\n","\n",sep="")
}
}
else cat("The plot was not captured!","\n",sep="")
}


print.Lord<-function(x,...){
res<-x
cat("\n")
cat("Detection of Differential Item Functioning using Lord's method","\n")
if (res$purification & is.null(res$anchor.names)) pur<-"with "
else pur<-"without "
if (is.null(res$c)) mod<-res$model
else mod<-"constrained 3PL"
cat("with ",mod," model and ",pur, "item purification","\n","\n",sep="")
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
if (max(loop)!=length(res$LordChi)) cat("(Note: no loop detected in less than ",res$nrPur,word,")","\n",sep="")
else cat("(Note: loop of length ",min((1:res$nrPur)[loop==length(res$LordChi)])," in the item purification process)","\n",sep="")
cat("WARNING: following results based on the last iteration of the purification","\n","\n")
}
else cat("Convergence reached after ",res$nrPur,word,"\n","\n",sep="")
}
if (is.null(res$anchor.names)) {
itk<-1:length(res$LordChi)
cat("No set of anchor items was provided", "\n", "\n")
}
else {
itk<-(1:length(res$LordChi))[!is.na(res$LordChi)]
cat("Anchor items (provided by the user):", "\n")
if (is.numeric(res$anchor.names)) mm<-res$names[res$anchor.names]
else mm<-res$anchor.names
mm <- cbind(mm)
rownames(mm) <- rep("", nrow(mm))
colnames(mm) <- ""
print(mm, quote = FALSE)
cat("\n", "\n")
}
cat("Lord's chi-square statistic:","\n","\n")
it<-rep("",length(res$LordChi))
df<-switch(res$model,"1PL"=1,"2PL"=2,"3PL"=3)
pval<-round(1-pchisq(res$LordChi,df),4)
symb<-symnum(pval,c(0,0.001,0.01,0.05,0.1,1),symbols=c("***","**","*",".",""))
m1<-cbind(round(res$LordChi[itk],4),pval[itk])
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
        for (i in 1:length(res$LordChi)) rn[i] <- paste("Item", i, sep = "")
        m2 <- rn
    }
        m2 <- cbind(m2[res$DIFitems])
rownames(m2)<-rep("",nrow(m2))
colnames(m2)<-""
print(m2,quote=FALSE)
cat("\n")
}
if (res$model=="1PL"){
cat("Effect size (ETS Delta scale):", "\n", "\n")
    cat("Effect size code:", "\n")
    cat(" 'A': negligible effect", "\n")
    cat(" 'B': moderate effect", "\n")
    cat(" 'C': large effect", "\n", "\n")
if (res$purification & is.null(res$anchor.names)) pars<-res$itemParFinal
else pars<-res$itemParInit
J<-nrow(pars)/2
mR<-pars[1:J,1]
mF<-itemRescale(pars[1:J,],pars[(J+1):(2*J),])[,1]
rr1<-round(mF-mR,4)
rr2<-round(-2.35*rr1,4)
    symb1 <- symnum(abs(rr2), c(0, 1, 1.5, Inf), symbols = c("A", 
        "B", "C"))
    matR2 <- cbind(rr1, rr2)[itk,]
    matR2 <- noquote(cbind(format(matR2, justify = "right"), 
        symb1[itk]))
    if (!is.null(res$names)) rownames(matR2) <- res$names[itk]
    else {
        rn <- NULL
        for (i in 1:nrow(matR2)) rn[i] <- paste("Item", i, sep = "")
        rownames(matR2) <- rn[itk]
    }
    colnames(matR2) <- c("mF-mR", "deltaLord", "")
    print(matR2)
    cat("\n")
    cat("Effect size codes: 0 'A' 1.0 'B' 1.5 'C'", "\n")
    cat(" (for absolute values of 'deltaLord')", "\n", "\n")
}
    if (!x$save.output) cat("Output was not captured!","\n")
    else {
if (x$output[2]=="default") wd<-paste(getwd(),"/",sep="")
else wd<-x$output[2]
fileName<-paste(wd,x$output[1],".txt",sep="")
cat("Output was captured and saved into file","\n"," '",fileName,"'","\n","\n",sep="")
}
}

