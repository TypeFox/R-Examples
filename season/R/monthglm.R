## monthglm.R
## fit a GLM using month as a factor
## option to add offset to control for uneven number of days
## Jan 2014

monthglm<-function(formula,data,family=gaussian(),refmonth=1,
                   monthvar='month',offsetmonth=FALSE,offsetpop=NULL){
  ## checks
  if (refmonth<1|refmonth>12){stop("Reference month must be between 1 and 12")}
  ## original call with defaults (see amer package)
  ans <- as.list(match.call())
  frmls <- formals(deparse(ans[[1]]))
  add <- which(!(names(frmls) %in% names(ans)))
  call<-as.call(c(ans, frmls[add]))

  monthvar=with(data,get(monthvar))
  cmonthvar=class(monthvar)
  ## If month is a character, create the numbers
  if(cmonthvar%in%c('factor','character')){
     if(cmonthvar=='character'){
        if(max(nchar(monthvar))==3){mlevels=substr(month.name,1,3)}else{mlevels=month.name}
        monthvar=factor(monthvar,levels=mlevels)
     }
     months=as.numeric(monthvar)
     data$month=months # add to data for flagleap
     months=as.factor(months)
     levels(months)[months]<-month.abb[months]
     months<-relevel(months.u,ref=month.abb[refmonth]) # set reference month
  }
  ## Transform month numbers to names
  if(cmonthvar%in%c('integer','numeric')){
    months.u<-as.factor(monthvar)  
    nums<-as.numeric(nochars(levels(months.u))) # Month numbers
    levels(months.u)[nums]<-month.abb[nums]
    months<-relevel(months.u,ref=month.abb[refmonth]) # set reference month
  }
  ## prepare data/formula
  parts<-paste(formula)
  f<-as.formula(paste(parts[2],parts[1],parts[3:length(formula)],'+months'))
  dep<-parts[2] # dependent variable
  days<-flagleap(data=data,report=FALSE,matchin=T) # get the number of days in each month
  l<-nrow(data)
  if(is.null(offsetpop)==FALSE){poff=with(data,eval(offsetpop))} else{poff=rep(1,l)} # 
  if(offsetmonth==TRUE){moff=days$ndaysmonth/(365.25/12)} else{moff=rep(1,l)} # days per month divided by average month length
###  data$off<-log(poff*moff)
  off<-log(poff*moff)  # 
  fit<-glm(formula=f,data=data,family=family,offset=off)
  ## return
  toret<-list()
  toret$call<-call
  toret$glm<-fit
  toret$fitted.values<-fitted(fit)
  toret$residuals<-residuals(fit)
  class(toret)<-'monthglm'
  return(toret)
}
