
clfs.strat <- function(data, stratf = NULL, maxx = NULL, com.est = TRUE, conf.int = FALSE, conf.int.level = NULL, no.iter = NULL, points = NULL, fig = TRUE, pvals = FALSE, pval.test = NULL)
{

LastContact.clfs <- data[,ncol(data)-1] - data[,1]

# compute the CLFS estimate for each level of the stratification factor:
strat.levels <- levels(as.factor(stratf[!is.na(data[,1])])) # levels of the stratification factor
pest.day <- matrix(0:maxx,maxx+1,1) # allocation of a matrix with point estimates (accompanied with confidence intervals) at each day
colnames(pest.day) <- "Day"
pest <- matrix(c(0,points),length(points)+1,1) # allocation of a matrix with point estimates (accompanied with confidence intervals) at the defined time points
colnames(pest) <- "Month"
no.risk <- matrix(c(0,points),length(points)+1,1) # allocation of a matrix with numbers at risk at the defined time points
colnames(no.risk) <- "Month"
sums <- NULL # allocation of matrix with data summary in each stratification group
for (i in 1:length(strat.levels)) {
  print(paste("Computation of the CLFS estimate for the level",i,"of the stratification factor."))
  inxi <- stratf==strat.levels[i]
  maxxi <- min(maxx,max(LastContact.clfs[inxi],na.rm=TRUE)) # maximum follow-up time in the data subsample or the maxx chosen by user
  res <- clfs.nostrat(data[inxi,], maxx=maxxi, com.est, conf.int, conf.int.level, no.iter, points=points[points<=floor(maxxi/(365/12))], fig=FALSE) # current (and common) leukaemia-free survival function estimates for the subsample
  colnames(res$pest.day) <- paste(strat.levels[i],colnames(res$pest.day),sep="-")
  colnames(res$pest) <- paste(strat.levels[i],colnames(res$pest),sep="-")
  colnames(res$no.risk) <- paste(strat.levels[i],colnames(res$no.risk),sep="-")
  pest.day <- cbind(pest.day,rbind(as.matrix(res$pest.day[1:(maxxi+1),2:ncol(res$pest.day)]),as.matrix(array(NA,c(maxx+1-(maxxi+1),ncol(res$pest.day)-1)))))
  pest <- cbind(pest,rbind(as.matrix(res$pest[,2:ncol(res$pest)]),as.matrix(array(NA,c(nrow(pest)-nrow(res$pest),ncol(res$pest)-1)))))
  no.risk <- cbind(no.risk,rbind(as.matrix(res$no.risk[,2:ncol(res$no.risk)]),as.matrix(array(0,c(nrow(no.risk)-nrow(res$no.risk),ncol(res$no.risk)-1)))))
  sums <- cbind(sums,as.numeric(as.matrix(res$summary[2])))
}
sums <- cbind(sums,c(rowSums(sums[1:3,]),floor(max(LastContact.clfs,na.rm=TRUE)/(365/12)),ceiling(max(rowSums(!is.na(data[,1:(ncol(data)-2)])))/2) ))
summary <- data.frame(cbind(c("The number of patients who achieved at least first disease remission","The number of patients with loss of the first disease remission","The number of patients who died after achievement of the first disease remission","Follow-up since achievement of the first disease remission (in months)","The maximum number of achieved disease remissions"),sums))
colnames(summary) <- c("",strat.levels,"Total")

# if conf.int=FALSE and com.est=FALSE, set column names of pest and pest.day:
if (!conf.int & !com.est) {
  colnames(pest.day) <- c("Day",paste(strat.levels,array("CLFS",length(strat.levels)),sep="-"))
  colnames(pest) <- c("Month",paste(strat.levels,array("CLFS",length(strat.levels)),sep="-"))
}

# if com.est=FALSE, set column names of no.risk:
if (!com.est) {
  colnames(no.risk) <- c("Month",paste(strat.levels,array("CLFS_Nrisk",length(strat.levels)),sep="-"))
}


# plot the CLFS and LFS estimates and their confidence intervals:
if (fig) {
  x=0:maxx
  yrs <- floor(maxx/365) # a number of years
  plot(0,1,pch='.',cex=0.01,xlab="Years after achievement of the first disease remission",ylab="Probability",axes=FALSE,xlim=c(0,maxx),ylim=c(0,1))
  axis(2,at=seq(0,1,0.2)) # set points in which tick-marks are drawn on the y-axis
  axis(1,at=seq(0,((yrs+1)*365),365),labels=seq(0,(yrs+1),1)) # set points in which tick-marks are drawn on the x-axis  
  color=c(1,2,4,5,8,6,3,7)
  if (com.est) {     
    if (conf.int) {
      for (i in 1:length(strat.levels)) {
        lines(x,pest.day[,1+i*6-5],type="S",lty=1,lwd=1,col=color[i]) # plot the lower confidence interval for the CLFS estimate
        lines(x,pest.day[,1+i*6-4],type="S",lty=1,lwd=2,col=color[i]) # plot the CLFS estimate
        lines(x,pest.day[,1+i*6-3],type="S",lty=1,lwd=1,col=color[i]) # plot the upper confidence interval for the CLFS estimate
        lines(x,pest.day[,1+i*6-2],type="S",lty=2,lwd=1,col=color[i]) # plot the lower confidence interval for the LFS estimate
        lines(x,pest.day[,1+i*6-1],type="S",lty=2,lwd=2,col=color[i]) # plot the LFS estimate
        lines(x,pest.day[,1+i*6],type="S",lty=2,lwd=1,col=color[i]) # plot the upper confidence interval for the LFS estimate
      }
      legend("bottomright",legend=c(paste(array("CLFS-",length(strat.levels)),strat.levels,array(paste(" with ",conf.int.level*100,"% conf. int.",sep=""),length(strat.levels)),sep=""),paste(array("LFS-",length(strat.levels)),strat.levels,array(paste(" with ",conf.int.level*100,"% conf. int.",sep=""),length(strat.levels)),sep="")),lwd=1,lty=c(array(1,length(strat.levels)),array(2,length(strat.levels))),col=c(color[1:length(strat.levels)],color[1:length(strat.levels)]),bty="n",cex=0.9)
    } else {
      for (i in 1:length(strat.levels)) {
        lines(x,pest.day[,1+i*2-1],type="S",lty=1,lwd=1,col=color[i]) # plot the CLFS estimate
        lines(x,pest.day[,1+i*2],type="S",lty=2,lwd=1,col=color[i]) # plot the LFS estimate
      }
      legend("bottomright",legend=c(paste(array("CLFS-",length(strat.levels)),strat.levels,sep=""),paste(array("LFS-",length(strat.levels)),strat.levels,sep="")),lwd=1,lty=c(array(1,length(strat.levels)),array(2,length(strat.levels))),col=c(color[1:length(strat.levels)],color[1:length(strat.levels)]),bty="n",cex=0.9)
    }
  } else {
    if (conf.int) {
      for (i in 1:length(strat.levels)) {
        lines(x,pest.day[,1+i*3-2],type="S",lty=1,lwd=1,col=color[i]) # plot the lower confidence interval for the CLFS estimate
        lines(x,pest.day[,1+i*3-1],type="S",lty=1,lwd=2,col=color[i]) # plot the CLFS estimate
        lines(x,pest.day[,1+i*3],type="S",lty=1,lwd=1,col=color[i]) # plot the upper confidence interval for the CLFS estimate
      }
      legend("bottomright",legend=paste(array("CLFS-",length(strat.levels)),strat.levels,array(paste(" with ",conf.int.level*100,"% conf. int.",sep=""),length(strat.levels)),sep=""),lwd=1,lty=1,col=c(color[1:length(strat.levels)]),bty="n",cex=0.9)
    } else {
      for (i in 1:length(strat.levels)) {
        lines(x,pest.day[,1+i],type="S",lty=1,lwd=1,col=color[i]) # plot the CLFS estimate
      }
      legend("bottomright",legend=strat.levels,lwd=1,lty=1,col=c(color[1:length(strat.levels)]),bty="n",cex=0.9)
    }
  }
}

# compute p-values for the comparison of the stratified curves at fixed points in time:
if (pvals) {
  if (length(strat.levels)==1) {
    pval <- NULL
  } else {   
    if (com.est) {
      pest.sel <- pest[,sort(c(1,seq(2,ncol(pest),6),seq(3,ncol(pest),6),seq(4,ncol(pest),6)))]
    } else {
      pest.sel <- pest
    }
    if (length(strat.levels)==2) {
      # compute p-values for the comparison of 2 stratified curves at fixed points in time:
      pval <- pvals.2cat(pest.sel,pval.test)
    } else {
      # compute p-values for the comparison of 3 or more stratified curves at fixed points in time:
      pval <- pvals.cat(pest.sel,pval.test)
    }
  }
} else {
  pval <- NULL
}

rownames(pest) <- NULL
rownames(pest.day) <- NULL

clfs.strat <- list(summary=summary,no.risk=no.risk,pest=pest,pest.day=pest.day,pval=pval)

}

