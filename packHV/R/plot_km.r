#' Kaplan-Meier plot with number of subjects at risk below
#'
#' Kaplan-Meier plot with number of subjects at risk below
#'
#' @param formula same formula than in \code{\link{survfit}} (\code{Surv(time,cens)~group} or \code{Surv(time,cens)~1}), where \code{cens} must equal to 0 (censorship) or 1 (failure)
#' @param data data frame with \code{time}, \code{cens} and \code{group}
#' @param test boolean, \code{TRUE} to compute and display the p-value of the log-rank test
#' @param xy.pvalue numeric vector of length 2, coordinates where to display the p-value of the log-rank test
#' @param conf.int boolean, \code{TRUE} to display the confidence interval of the curve(s)
#' @param times.print numeric vector, times at which to display the numbers of subjects at risk
#' @param legend character string (\code{"bottomright"} for example) or numeric vector (\code{c(x,y)}), where to place the legend of the curve(s)
#' @param nrisk.labels character vector to modify the levels of \code{group} in the table below the curve(s)
#' @param xlab character string, label of the time axis
#' @param ylab character string, label of the y axis
#' @param ylim numeric vector of length 2, minimum and maximum of the y-axis
#' @param left integer, size of left margin
#' @param bottom integer, number of lines in addition of the table below the graph
#' @param cex.mtext numeric, size of the numbers of subjects at risk
#' @param lwd width of the Kaplan-Meier curve(s)
#' @param lty type of the Kaplan-Meier curve(s)
#' @param col color(s) of the Kaplan-Meier curve(s)
#' @param \dots other arguments to be passed in \code{\link{plot.survfit}}
#' @return None
#' @author Hugo Varet
#' @examples
#' cgd$time=cgd$tstop-cgd$tstart
#' plot_km(Surv(time,status)~sex,data=cgd,col=c("blue","red"))

# last update: july 01, 2013

plot_km=function(formula,data,test=TRUE,xy.pvalue=NULL,conf.int=FALSE,times.print=NULL,nrisk.labels=NULL,legend=NULL,
                 xlab=NULL,ylab=NULL,ylim=c(0,1.02),left=4.5,bottom=5,cex.mtext=par("cex"),lwd=2,lty=1,col=NULL,...){
  # formula: Surv(temps,cens)~groupe, ces variables étant dans data
  # test: TRUE to perform a log-rank test (set to FALSE if only one level in groupe)
  # xy.pvalue: vector of length 2 (coordinates where to place the p-value of the test)
  # conf.int: TRUE to add the confidence interval(s) of the curves
  # times.print: temps auxquels les n at risk s'affichent
  # nrisk.labels: character vector to customize the levels of groupe
  # legend: where to display the legend ("topright","bottomright" or "bottomleft" or c(x,y))
  # left: nb de lignes dans la marge de gauche (controle donc la largeur)
  # bottom: nb de lignes en plus de celles du tableau dans la marge du bas (utile pour juxtaposer pls graph)
  # attention: dans formula, il ne faut que des noms de variables, i.e. pas de cut(x,...)
  formula <- deparse(substitute(formula))

  # extracting the variable names from formula
  formula <- paste(formula, collapse=" ")
  formula <- gsub("[[:space:]]+", " ", formula)
  varnames <- gsub("\\<Surv\\>|\\(|\\)| \\?\\(.*\\)", "", formula)
  varnames <- strsplit(varnames,  "[[:space:]]*(\\+|,|~)[[:space:]]*")
  varnames <- gsub("[[:space:]]+", "", unlist(varnames))

  if (any(!I(setdiff(varnames,"1") %in% names(data)))){
    stop(paste(paste(varnames,collapse=" and/or ")," not in ",deparse(substitute(data)),"\n",sep=""))
  }
  
  temps=data[,varnames[1]]
  cens=data[,varnames[2]]
  
  if (varnames[3]=="1"){
    groupe=factor(rep("# at risk",nrow(data)))
    name_at_risk=""
    test=FALSE; legend=NULL;
  } else{
    groupe=factor(data[,varnames[3]])
    if (length(levels(factor(data[,varnames[3]])))==1){
      stop(paste(varnames[3]," has only one level\n",sep=""))
    }
    name_at_risk="# at risk"
  }
    
  d=data.frame(temps,cens,groupe)
  if (any(is.na(d$temps) | is.na(d$cens) | is.na(d$groupe))){
    cat(paste(sum(is.na(d$temps) | is.na(d$cens) | is.na(d$groupe))," rows deleted due to missing values\n",sep=""))
  }
  d=d[I(!is.na(d$temps) & !is.na(d$cens) & !is.na(d$groupe)),]
  
  if (is.null(xlab)){xlab=varnames[1]}
  if (is.null(ylab)){ylab="Survival"}
  if (is.null(col)){col=1:nlevels(d$groupe)}
  
#  mar=c(bottom, left, top, right)
  par(mar=c(1+nlevels(groupe)+bottom, left, 4, 3) + 0.1,xaxs="i",yaxs="i")
  plot(survfit(Surv(temps,cens)~groupe,data=d),conf.int=conf.int,xlab=xlab,ylab=ylab,lwd=lwd,lty=lty,col=col,ylim=ylim,...)
  
  # affichage de la légende
  if (!is.null(legend)){
    legend(legend[1],if(length(legend)==2){legend[2]}else{NULL},legend=levels(d$groupe),col=col,lwd=lwd,lty=lty)
  }
  
  # calcul des n at risk
  if (is.null(times.print)){
    times.print=axis(1,labels=FALSE,tick=FALSE)
    times.print=times.print[times.print>=0]
  }
  n.risk=matrix(NA,nrow=1+nlevels(groupe),ncol=length(times.print),dimnames=list(c("temps",levels(groupe)),NULL))
  n.risk[1,]=times.print
  for (lev in levels(groupe)){
    tmp=d[d$groupe==lev,]
    tmp2=summary(survfit(Surv(temps,cens)~1,data=tmp),times.print)$n.risk
    if (length(tmp2)<length(times.print)){
      n.risk[lev,]=c(tmp2,rep(0,length(times.print)-length(tmp2)))
    } else{
      n.risk[lev,]=tmp2
    }
  }
  if (!is.null(nrisk.labels)){rownames(n.risk)[-1]=nrisk.labels}
    
  # affichage des n at risk
  range=range(axis(1,labels=FALSE,tick=FALSE))
  mtext(side=1,at=-0.065*(range[2]-range[1])+range[1],line=4,name_at_risk,cex=cex.mtext,adj=1)
  for (i in 2:nrow(n.risk)){
    mtext(side=1,at=times.print,line=i+3,n.risk[i,],cex=cex.mtext)
    mtext(side=1,at=-0.065*(range[2]-range[1])+range[1],line=i+3,rownames(n.risk)[i],cex=cex.mtext,adj=1)
  } 
  
  # Log-Rank test
  if (test){
    diff=survdiff(Surv(temps,cens)~groupe,data=d)
    p.value=1-pchisq(diff$chisq,df=length(levels(d$groupe))-1)
    p.value=paste("p",ifelse(p.value<0.001,"<0.001",paste("=",round(p.value,3),sep="")),sep="")
    if (!is.null(xy.pvalue)){
      text(xy.pvalue[1],xy.pvalue[2],p.value,adj=c(0,0))
    } else{
      text((range[2]-range[1])/20 + range[1],0.05,p.value,adj=c(0,0))
    }
  }  
}


#library(survival)
#cgd$time=cgd$tstop-cgd$tstart
#par(mfrow=c(2,2))
#plot_km(Surv(time,status)~sex,data=cgd,col=c("blue","red"),xaxt="n",xlim=c(0,200),main="xlim=c(0,200)")
#plot_km(Surv(time,status)~sex,data=cgd,col=c("blue","red"),xlim=c(0,200),main="xlim=c(0,200)")
#plot_km(Surv(time,status)~sex,data=cgd,col=c("blue","red"),xaxt="n",xlim=c(0,500),main="xlim=c(0,500)")
#plot_km(Surv(time,status)~sex,data=cgd,col=c("blue","red"),xlim=c(0,500),main="xlim=c(0,500)")     
#
#cgd$time=cgd$tstop-cgd$tstart
#plot_km(Surv(time,status)~sex,data=cgd,col=c("blue","red"),xaxt="n",xlim=c(0,510),main="xlim=c(0,510)")
#axis(1,at=c(0,100,200,300,400,500))
#
#par(mfrow=c(1,2))
#plot_km(Surv(time,status)~sex,data=cgd,col=c("blue","red"),xlim=c(200,500),main="xlim=c(200,500)")     
#plot_km(Surv(time,status)~sex,data=cgd,col=c("blue","red"),xlim=c(-200,500),main="xlim=c(-200,500)")     
#         
#plot_km(Surv(time,status)~sex,data=cgd,col=c("blue","red"),xlim=c(200,500),main="xlim=c(200,500)")     

