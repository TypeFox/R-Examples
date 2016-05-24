# cosinor.R
# cosinor function using a GLM
# available link functions = identity, log, logit, cloglog
# date = date for daily data, month for monthly data
# type =  monthly/weekly/daily/hourly
# phase results based on 1 cycle per year
# Aug 2014

cosinor<-function(formula, date, data,family=gaussian(), alpha=0.05,
                  cycles=1, rescheck=FALSE, type='daily', offsetmonth=FALSE,
                  offsetpop=NULL,text=TRUE){
  ## checks
  classes = lapply(data, class) # classes of all variables
  this.class=as.character(classes[which(names(data)==date)]) # also used later
  if (!is.logical(offsetmonth)){
    stop("Error: 'offsetmonth' must be of type logical")}
  if (type!='daily'&type!='weekly'&type!='monthly'&type!='hourly'){stop("type must be daily, weekly, monthly or hourly")}
  if (type=='hourly'&length(grep('POSIXct',this.class))==0){
    stop("date variable must be of class POSIXct when type='hourly'")}
  if (type=='daily'&this.class!='Date'){
    stop("date variable must be of class Date when type='daily'")}
  if (alpha<=0|alpha>=1){stop("alpha must be between 0 and 1")}
  if (type=='hourly'&offsetmonth==TRUE){
    stop("do not use monthly offset for hourly data")}
  
  ## original call with defaults (see amer package)
  link<-family$link
  ans <- as.list(match.call())
  frmls <- formals(deparse(ans[[1]]))
  add <- which(!(names(frmls) %in% names(ans)))
  call<-as.call(c(ans, frmls[add],link=link))

  ## make the formula
  parts<-paste(formula)
  f<-as.formula(paste(parts[2],parts[1],parts[3:length(formula)],'+cosw+sinw'))

  ## get the year/hour fraction
  to.frac=subset(data,select=date)[,1]
  if(type=='hourly'){
    number = as.numeric(format(to.frac,'%H'))+(as.numeric(format(to.frac,'%M'))/60)+(as.numeric(format(to.frac,'%S'))/60*60)
    frac = number/24
  }
  if(type!='hourly'){
    class(to.frac)=this.class # return to class (needed for date class)
    frac = yrfraction(to.frac,type=type) # 
  }
  data$cosw<-cos(frac*2*pi*cycles)
  data$sinw<-sin(frac*2*pi*cycles)
  newdata<-data.frame(cosw=data$cosw,sinw=data$sinw) # used later
  poff=rep(1,nrow(data))
  if(is.null(offsetpop)==FALSE){poff=offsetpop}
  moff=rep(1,nrow(data))
  if(offsetmonth==TRUE){
    days = flagleap(data=data, report=FALSE, matchin=TRUE) # get the number of days in each month
    moff = days$ndaysmonth/(365.25/12) # days per month divided by average month length
  }
  offset <- log(poff*moff)
  # generalized linear model
  model<-glm(f, data=data, family=family, offset=offset)
  s<-summary(model)
  res<-residuals(model)

  ## create predicted data (intercept + sinusoid)
  cnames<-row.names(s$coefficients)
  cindex<-sum(as.numeric(cnames=='cosw')*(1:length(cnames)))
  sindex<-sum(as.numeric(cnames=='sinw')*(1:length(cnames)))
  fitted=fitted(model) # standard fitted values
  pred<-s$coefficients[1,1]+(s$coefficients[cindex,1]*newdata$cosw)+
    (s$coefficients[sindex,1]*newdata$sinw)
  # back-transform:
  if (s$family$link=='log'){pred<-exp(pred)}
  if (s$family$link=='logit'){pred<-exp(pred)/(1+exp(pred))}
  if (s$family$link=='cloglog'){pred<-1-exp(-exp(pred))}
  # return:
  toret<-list()
  toret$call<-call
  toret$glm<-model # changed to model rather than summary
  toret$fitted.plus<-fitted
  toret$fitted.values<-pred
  toret$residuals<-res
  toret$date<-date
  class(toret)<-'Cosinor'
  return(toret)
}
