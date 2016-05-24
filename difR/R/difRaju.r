difRaju<-function(Data, group, focal.name, model, c=NULL, engine="ltm", discr=1,irtParam=NULL, same.scale=TRUE, anchor=NULL, alpha=0.05, signed=FALSE, purify=FALSE, nrIter=10,save.output=FALSE, output=c("out","default"))
{
internalRaju<-function(){
if (!is.null(irtParam)){
nrItems<-nrow(irtParam)/2
m0<-irtParam[1:nrItems,]
m1<-irtParam[(nrItems+1):(2*nrItems),]
dataName=rownames(m0)
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
model=switch(mod,"2"="1PL","5"="2PL","6"="3PL")
if (ncol(irtParam)!=6) Guess<-NULL
else{
Guess<-irtParam[1:nrItems,6]
if (length(unique(round(Guess,4)))==1) Guess<-unique(round(Guess,4))
}
itemParInit<-irtParam
estPar<-FALSE
}
else{
     if (length(group) == 1) {
           if (is.numeric(group)==TRUE) {
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
if (is.null(Guess)){
m0<-switch(model,"1PL"=itemParEst(d0,model="1PL",engine=engine,discr=discr),"2PL"=itemParEst(d0,model="2PL"))
m1<-switch(model,"1PL"=itemParEst(d1,model="1PL",engine=engine,discr=discr),"2PL"=itemParEst(d1,model="2PL"))
}
else{
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
dataName<-colnames(DATA)
itemParInit<-rbind(m0,m1)
estPar<-TRUE
}

if (!purify | !is.null(anchor)) {
STATS<-RajuZ(m0,m1p,signed=signed)$res[,3]
if (max(abs(STATS))<=qnorm(1-alpha/2)) DIFitems<-"No DIF item detected"
else DIFitems<-(1:nrItems)[abs(STATS)>qnorm(1-alpha/2)]
RES<-list(RajuZ=STATS,alpha=alpha,thr=qnorm(1-alpha/2),DIFitems=DIFitems,signed=signed,purification=purify,model=model,c=Guess,engine=engine,discr=discr,itemParInit=itemParInit,estPar=estPar,names=dataName,anchor.names=dif.anchor,save.output=save.output,output=output)
if (!is.null(anchor) & (RES$estPar | (!RES$estPar & !same.scale))) {
RES$RajuZ[ANCHOR]<-NA
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
stats1<-RajuZ(m0,m1p,signed=signed)$res[,3]
if (max(abs(stats1))<=qnorm(1-alpha/2)){
DIFitems<-"No DIF item detected"
noLoop<-TRUE
itemParFinal=rbind(m0,m1p)
RES<-list(RajuZ=stats1,alpha=alpha,thr=qnorm(1-alpha/2),DIFitems=DIFitems,signed=signed,purification=purify,nrPur=nrPur,difPur=difPur,convergence=noLoop,model=model,c=Guess,engine=engine,discr=discr,itemParInit=itemParInit,itemParFinal=itemParFinal,estPar=estPar,names=dataName,anchor.names=NULL, save.output=save.output,output=output)
}
else{ 
dif<-(1:nrItems)[abs(stats1)>qnorm(1-alpha/2)]
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
if (is.null(dif)) nodif<-1:nrItems
else{
for (i in 1:nrItems){
if (sum(i==dif)==0) nodif<-c(nodif,i)
}
}
stats2<-RajuZ(m0,itemRescale(m0,m1,items=nodif),signed=signed)$res[,3]
if (max(abs(stats2))<=qnorm(1-alpha/2)) dif2<-NULL
else dif2 <- (1:nrItems)[abs(stats2)>qnorm(1-alpha/2)]
difPur<-rbind(difPur,rep(0,nrItems))
difPur[nrPur+1,dif2]<-1
if (length(dif)!=length(dif2)) dif<-dif2
else{
dif<-sort(dif)
dif2<-sort(dif2)
if (sum(dif==dif2)==length(dif)){
noLoop<-TRUE
itemParFinal<-rbind(m0,itemRescale(m0,m1,items=nodif))
break
}
else dif<-dif2
}
}
}
if (is.null(difPur)==FALSE){
ro<-co<-NULL
for (ir in 1:nrow(difPur)) ro[ir]<-paste("Step",ir-1,sep="")
for (ic in 1:ncol(difPur)) co[ic]<-paste("Item",ic,sep="")
rownames(difPur)<-ro
colnames(difPur)<-co
}
RES<-list(RajuZ=stats2,alpha=alpha,thr=qnorm(1-alpha/2),DIFitems=dif2,signed=signed,purification=purify,nrPur=nrPur,difPur=difPur,convergence=noLoop,model=model,c=Guess,engine=engine,discr=discr,itemParInit=itemParInit,itemParFinal=itemParFinal,estPar=estPar,names=dataName,anchor.names=NULL, save.output=save.output,output=output)
}
}
class(RES)<-"Raj"
return(RES)
}
resToReturn<-internalRaju()
if (save.output){
if (output[2]=="default") wd<-paste(getwd(),"/",sep="")
else wd<-output[2]
fileName<-paste(wd,output[1],".txt",sep="")
capture.output(resToReturn,file=fileName)
}
return(resToReturn)
}


# METHODS
plot.Raj<-function(x,pch=8,number=TRUE,col="red", save.plot=FALSE,save.options=c("plot","default","pdf"),...) 
{
internalRaju<-function(){
res <- x
title<-paste("Raju's method (",res$model,")",sep="")
if (!number) {
plot(res$RajuZ,xlab="Item",ylab="Raju's statistic",ylim=c(min(c(res$RajuZ,-res$thr)-1,na.rm=TRUE),max(c(res$RajuZ,res$thr)+1,na.rm=TRUE)),pch=pch,main=title)
if (!is.character(res$DIFitems)) points(res$DIFitems,res$RajuZ[res$DIFitems],pch=pch,col=col)
}
else {
plot(res$RajuZ,xlab="Item",ylab="Raju's statistic",ylim=c(min(c(res$RajuZ,-res$thr)-1,na.rm=TRUE),max(c(res$RajuZ,res$thr)+1,na.rm=TRUE)),col="white",main=title)
text(1:length(res$RajuZ),res$RajuZ,1:length(res$RajuZ))
if (!is.character(res$DIFitems)) text(res$DIFitems,res$RajuZ[res$DIFitems],res$DIFitems,col=col)
}
abline(h=res$thr)
abline(h=-res$thr)
abline(h=0,lty=2)
}
internalRaju()
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
internalRaju()
}
dev.off()
}
if (plotype==2){
{
jpeg(filename=fileName)
internalRaju()
}
dev.off()
}
cat("The plot was captured and saved into","\n"," '",fileName,"'","\n","\n",sep="")
}
}
else cat("The plot was not captured!","\n",sep="")
}



print.Raj<-function(x,...){
res <- x
if (is.character(res)) cat("Error: unappropriate model specification!","\n")
else{
cat("\n")
cat("Detection of Differential Item Functioning using Raju's method","\n")
if (res$purification & is.null(res$anchor.names)) pur<-"with "
else pur<-"without "
cat("with ",res$model," model and ",pur, "item purification","\n","\n",sep="")
if (res$signed) cat("Type of Raju's Z statistic: based on signed area","\n","\n")
else cat("Type of Raju's Z statistic: based on unsigned area","\n","\n")
if (res$estPar){
if (res$model!="1PL" | res$engine=="ltm") cat("Engine 'ltm' for item parameter estimation","\n","\n")
else cat("Engine 'lme4' for item parameter estimation","\n","\n")
}
if (res$model=="1PL" & res$engine=="ltm") {
if (is.null(res$discr)) cat("Common discrimination parameter: estimated from 'ltm'","\n","\n")
else cat("Common discrimination parameter: fixed to ",res$discr,"\n","\n",sep="")
}
if (!is.null(res$c)){
if (length(res$c)==1) cat("Common pseudo-guessing value:",res$c,"\n","\n")
else 
{
mat<-cbind(res$c)
if (!is.null(res$names)) rownames(mat)<-res$names
else{
rn<-NULL
for (i in 1:nrow(mat)) rn[i]<-paste("Item",i,sep="")
rownames(mat)<-rn
}
colnames(mat)<-"c"
cat("Common pseudo-guessing values:","\n","\n")
print(mat)
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
if (max(loop)!=length(res$RajuZ)) cat("(Note: no loop detected in less than ",res$nrPur,word,")","\n",sep="")
else cat("(Note: loop of length ",min((1:res$nrPur)[loop==length(res$RajuZ)])," in the item purification process)","\n",sep="")
cat("WARNING: following results based on the last iteration of the purification","\n","\n")
}
else cat("Convergence reached after ",res$nrPur,word,"\n","\n",sep="")
}
if (is.null(res$anchor.names)) {
itk<-1:length(res$RajuZ)
cat("No set of anchor items was provided", "\n", "\n")
}
else {
itk<-(1:length(res$RajuZ))[!is.na(res$RajuZ)]
cat("Anchor items (provided by the user):", "\n")
if (is.numeric(res$anchor.names)) mm<-res$names[res$anchor.names]
else mm<-res$anchor.names
mm <- cbind(mm)
rownames(mm) <- rep("", nrow(mm))
colnames(mm) <- ""
print(mm, quote = FALSE)
cat("\n", "\n")
}
cat("Raju's statistic:","\n","\n")
pval<-round(2*(1-pnorm(abs(res$RajuZ))),4)
symb<-symnum(pval,c(0,0.001,0.01,0.05,0.1,1),symbols=c("***","**","*",".",""))
m1<-cbind(round(res$RajuZ[itk],4),pval[itk])
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
cat("\n","Detection thresholds: ",-round(res$thr,4)," and ",round(res$thr,4)," (significance level: ",res$alpha,")","\n","\n",sep="")
if (is.character(res$DIFitems)) cat("Items detected as DIF items:",res$DIFitems,"\n","\n")
else {
cat("Items detected as DIF items:","\n")
   if (!is.null(res$names)) m2 <- res$names
    else {
        rn <- NULL
        for (i in 1:length(res$RajuZ)) rn[i] <- paste("Item", i, sep = "")
        m2 <- rn
    }
        m2 <- cbind(m2[res$DIFitems])
rownames(m2)<-rep("",nrow(m2))
colnames(m2)<-""
print(m2,quote=FALSE)
cat("\n")
}
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
    colnames(matR2) <- c("mF-mR", "deltaRaju", "")
    print(matR2)
    cat("\n")
    cat("Effect size codes: 0 'A' 1.0 'B' 1.5 'C'", "\n")
    cat(" (for absolute values of 'deltaRaju')", "\n", "\n")
}
    if (!x$save.output) cat("Output was not captured!","\n")
    else {
if (x$output[2]=="default") wd<-paste(getwd(),"/",sep="")
else wd<-x$output[2]
fileName<-paste(wd,x$output[1],".txt",sep="")
cat("Output was captured and saved into file","\n"," '",fileName,"'","\n","\n",sep="")
}
}

