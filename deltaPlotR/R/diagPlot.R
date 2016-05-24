diagPlot<-function(x,pch=2,pch.mult=17,axis.draw=TRUE,thr.draw=FALSE,dif.draw=c(1,3),print.corr=FALSE,xlim=NULL,ylim=NULL,xlab=NULL,ylab=NULL,main=NULL,save.plot=FALSE,save.options=c("plot","default","pdf")){
internalDP<-function(){
ra<-range(x$Deltas)
rm<-ifelse(round(ra[1])<=0,0,round(ra[1])-1)
rM<-ifelse(round(ra[2]+2)>=26,26,round(ra[2])+2)
if (is.null(xlim)) xl<-c(rm,rM)
else xl<-xlim
if (is.null(ylim)) yl<-c(rm,rM)
else yl<-ylim
if (is.null(xlab)) xla<-"Reference group"
else xla<-xlab
if (is.null(ylab)) yla<-"Focal group"
else yla<-ylab
plot(x$Deltas[,1],x$Deltas[,2],pch=pch,xlim=xl,ylim=yl,xlab=xla,ylab=yla,main=main)
if (axis.draw){
if (is.null(dim(x$axis.par))) pars<-x$axis.par
else pars<-x$axis.par[nrow(x$axis.par),]
abline(pars[1],pars[2])
}
if (nrow(x$Deltas)!=nrow(unique(x$Deltas))){
N<-nrow(x$Deltas)
id<-rep(0,N)
for (i in 1:(N-1)){
for (j in (i+1):N){
if (sum(x$Deltas[i,]-x$Deltas[j,])==0) id[c(i,j)]<-id[c(i,j)]+1
}}
points(x$Deltas[id>0,1],x$Deltas[id>0,2],pch=pch.mult)
}
if (!is.character(x$DIFitems)) points(x$Deltas[x$DIFitems,1],x$Deltas[x$DIFitems,2],pch=dif.draw[1],cex=dif.draw[2])
if (thr.draw) {
if (is.null(dim(x$axis.par))) pars<-x$axis.par
else pars<-x$axis.par[nrow(x$axis.par),]
th<-x$thr[length(x$thr)]
abline(pars[1]+th*sqrt(pars[2]^2 + 1),pars[2],lty=2)
abline(pars[1]-th*sqrt(pars[2]^2 + 1),pars[2],lty=2)
}
if (print.corr){
rho<-cor(x$Deltas[,1],x$Deltas[,2])
legend(xl[1],yl[2],substitute(r[Delta] == x,list(x=round(rho,3))),bty="n")
}
 }
    internalDP()
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
            if (save.options[2] == "default") 
                wd <- file.path(getwd())
            else wd <- save.options[2]
nameF<-paste(save.options[1], switch(plotype, `1` = ".pdf", `2` = ".jpg"), sep = "")
            nameFile <- file.path(wd,nameF)
            if (plotype == 1) {
                {
                  pdf(file = nameFile)
                  internalDP()
                }
                dev.off()
            }
            if (plotype == 2) {
                {
                  jpeg(filename = nameFile)
                  internalDP()
                }
                dev.off()
            }
            cat("The plot was captured and saved into", "\n", 
                " '", nameFile, "'", "\n", "\n", sep = "")
        }
    }
    else cat("The plot was not captured!", "\n", sep = "")
}

