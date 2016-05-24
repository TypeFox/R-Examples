### Plot the monthly data

plot.monthglm<-function(x, alpha=0.05, ylim=NULL, ...){
  ## Checks
  if (class(x)!="monthglm"){stop("Object must be of class 'monthglm'")} 
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  ## y-axis limits
  if(is.null(ylim)==FALSE){this.y.lim=ylim}
  ## Create CIs
  refer=NA # reference level - only used for binomial or poisson
  z<-qnorm(1-(alpha/2))
  s<-summary(x$glm)
  type<-as.character(x$call$family)[1]
  out<-as.data.frame(matrix(data=NA,nrow=nrow(s$coef),ncol=3))
  names(out)<-c('mean','lower','upper')
  row.names(out)<-row.names(s$coef)
  out$mean<-s$coef[,1]
  out$lower<-s$coef[,1]-(z*s$coef[,2])
  out$upper<-s$coef[,1]+(z*s$coef[,2])

  ## Exponentiate the results if rate or odds ratio
  if (type=="poisson"|type=="binomial"){
    out$mean<-exp(out$mean)
    out$lower<-exp(out$lower)
    out$upper<-exp(out$upper)
    refer=1
  }
  index<-grep("months",row.names(out),ignore.case=TRUE,value=FALSE)
  toplot<-out[index,] # Select months

  ## Get month names
  compress<-gsub('months','',row.names(toplot))
  order<-vector(length=nrow(toplot),mode='numeric')
  for (i in 1:nrow(toplot)){
    order[i]<-sum(as.numeric(month.abb==compress[i])*(1:12))
  }

  ## plot
  ymin<-min(c(toplot$lower,refer),na.rm=T) # include reference
  ymax<-max(c(toplot$upper,refer),na.rm=T) # include reference
  if(is.null(ylim)){this.y.lim=c(ymin,ymax)}
  plot(order,toplot$mean,xaxt='n',xlab='',ylab='',
       ylim=this.y.lim,xlim=c(1,12),...)
  month.lab<-vector(mode='character',length=nrow(toplot))
  for (i in 1:nrow(toplot)){
    lines(c(order[i],order[i]),c(toplot$lower[i],toplot$upper[i]))
    months.num<-as.numeric(nochars(row.names(toplot)[i]))
    month.lab[i]<-month.abb[months.num]
  }
  month.lab<-substr(month.abb,1,1)
  axis(side=1,labels=month.lab,at=1:12)
  if (type=="poisson"|type=="binomial"){
    lines(c(1,12),c(1,1),lty=2) # reference list
    points(x$call$refmonth,1) # reference point
  }
  par(op) # restore graphic settings
}

