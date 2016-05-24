## casecross.R
## time-stratified case-crossover
## Oct 2011
## assumes date variable is called 'date'
## quicker version

casecross<-function(formula, data, exclusion=2, stratalength=28,
                    matchdow=FALSE, usefinalwindow=FALSE, matchconf='',
                    confrange=0,stratamonth=FALSE){
  outcome <- dow <- case <-timex <- dow.x <- dow.y <- matchday.x <- matchday.y <- windownum.x <- windownum.y <- NULL # Setting some variables to NULL first (for R CMD check)
  thisdata<-data
  ## Checks
  if (class(thisdata$date)!="Date"){
    stop("date variable must be in date format, see ?Dates")} 
  if (exclusion<0){stop("Minimum value for exclusion is zero")} 
  parts<-paste(formula)
  dep<-parts[2] # dependent variable
  indep<-parts[3] # dependent variable
  if (length(formula)<=2){stop("Must be at least one independent variable")} 
  ## original call with defaults (see amer package)
  ans <- as.list(match.call())
  frmls <- formals(deparse(ans[[1]]))
  add <- which(!(names(frmls) %in% names(ans)))
  call<-as.call(c(ans, frmls[add]))
  thisdata$dow<-as.numeric(format(thisdata$date,'%w'));
  ## Slim down the data 
  f<-as.formula(paste(parts[2],parts[1],parts[3],'+date+dow'))
  if (substr(matchconf,1,1)!=""){
    f<-as.formula(paste(dep,"~",indep,'+date+dow+',matchconf))
  }
  datatouse<-model.frame(f,data=thisdata,na.action=na.omit) # remove cases with missing covariates
  ## Check for irregularly spaced data
  if(any(diff(datatouse$date)>1)){
     cat('Note, irregularly spaced data...\n')
     cat('...check your data for missing days\n')
  }
  datediff<-as.numeric(datatouse$date)-min(as.numeric(thisdata$date)) # use minimum data in entire sample
  time<-as.numeric(datediff)+1 # used as strata number

  ## Create strata
  if (stratamonth==TRUE){
    month<-as.numeric(format(datatouse$date,'%m'));
    year<-as.numeric(format(datatouse$date,'%Y'));
    matchday<-as.numeric(format(datatouse$date,'%d'));
    yrdiff<-year-min(year);
    windownum<-(yrdiff*12)+month;
  }
  if (stratamonth==FALSE){
    ## Get the earliest time and difference all dates from this time
    ## Increase strata windows in jumps of 'stratalength'
    windownum<-floor(datediff/stratalength)+1
    nwindows<-floor(nrow(thisdata)/stratalength)+1
    matchday<-datediff-((windownum-1)*stratalength)+1 # Day number in strata
    ## Exclude the last window if it is less than 'stratalength'
    lastwindow<-datatouse[datatouse$windownum==nwindows,]
    if (nrow(lastwindow)>0){ # only apply to data sets with some data in the final window
      lastlength<-max(time[windownum==nwindows])-
        min(time[windownum==nwindows])+1
      if (lastlength<stratalength&usefinalwindow==FALSE) datatouse <-
        datatouse[windownum<nwindows,]
    }
  }
  ## Create the case data
  n<-nrow(datatouse)
  cases<-datatouse
  cases$case<-1 # binary indicator of case
  cases$timex<-1 # Needed for conditional logistic regression
  cases$windownum<-windownum
  cases$time<-time
  cases$diffdays<-NA
  cases$matchday<-matchday
  posout<-sum(as.numeric(names(datatouse)==as.character(f[2]))*
              (1:ncol(datatouse))) # get the position of the dependent variable
  cases$outcome<-datatouse[,c(posout)]
  # October 2011, removed nonzerocases
  # Create a case number for matching
  if (substr(matchconf,1,1)==""){
    cases.tomerge<-subset(cases,select=c(matchday,time,outcome,windownum,dow))}
  if (substr(matchconf,1,1)!=""){
     also<-sum(as.numeric(names(cases)==matchconf)*(1:length(names(cases))))
     cases.tomerge<-subset(cases,
                           select=c(matchday,time,outcome,windownum,dow,also))
  }
  ncases<-nrow(cases)
  cases.tomerge$casenum<-1:ncases
  # Duplicate case series to make controls
  maxwindows<-max(cases$windownum)
  rowstorep<-NA
  casenum<-NA
  # Fix for missing windows (thanks to Yuming)
  windowrange<-as.numeric(levels(as.factor(windownum)))
  for (k in windowrange){
    small=min(cases$time[cases$windownum==k])
    large=max(cases$time[cases$windownum==k])
    these<-rep(small:large,large-small+1)
    rowstorep<-c(rowstorep,these)
    casenum<-c(casenum,these[order(these)])
  }
  controls<-cases[rowstorep[2:length(rowstorep)],] # can fall over if there's missing data
  controls<-subset(controls,select=c(-case,-timex,-time,-outcome))
  # Replace case number
  controls$casenum<-casenum[2:length(rowstorep)]
  # Merge cases with controls by case number
  controls<-merge(controls,cases.tomerge,by='casenum')
  controls<-controls[controls$windownum.x==controls$windownum.y,] # must be in same stratum window
  controls$case<-0 # binary indicator of case
  controls$timex<-2 # Needed for conditional logistic regression
  controls$diffdays<-abs(controls$matchday.x-controls$matchday.y)
  controls<-controls[controls$diffdays>exclusion,] # remove the exclusion window
  # match on day of the week
  if (matchdow==TRUE){controls<-controls[controls$dow.x==controls$dow.y,]} 
  # match on a confounder
  if (substr(matchconf,1,1)!=""){
     one<- paste(matchconf,'.x',sep='')
     two<- paste(matchconf,'.y',sep='')
     find1<-grep(one,names(controls))
     find2<-grep(two,names(controls))
     matchdiff<-abs(controls[,find1]-controls[,find2])
     controls<-controls[matchdiff<=confrange,] 
     controls<-subset(controls,select=c(-casenum,-dow.x,-dow.y,-matchday.x,-matchday.y,-windownum.x,-windownum.y,-find1,-find2))
     findc<-sum(as.numeric(names(cases)==matchconf)*(1:length(names(cases))))
     final.cases<-subset(cases,select=c(-dow,-matchday,-windownum,-findc))
  }
  if (substr(matchconf,1,1)==""){
    controls<-subset(controls,select=c(-casenum,-dow.x,-dow.y,-matchday.x,-matchday.y,-windownum.x,-windownum.y))
    final.cases<-subset(cases,select=c(-dow,-matchday,-windownum))
  }
  finished<-rbind(final.cases,controls)
  ## Remove empty controls
  finished<-finished[finished$outcome>0,]
  ## Count the number of control days without a case day, and the total number of cases
  onlycntl<-finished[finished$case==0,]
  ncases<-nrow(table(onlycntl$time))
  which.times=unique(onlycntl$time)
  extra.only=final.cases[final.cases$time%in%which.times,]
  ncontrols<-round(mean(as.numeric(table(onlycntl$time))),1)
  ## Run the conditional logistic regression
  finalformula<-as.formula(paste('Surv(timex,case)~',indep,'+strata(time)'))
  c.model<-coxph(finalformula,
                 weights=outcome,
                 data=finished,method=c("breslow"))
  toret<-list()
  toret$call<-call
  toret$c.model<-c.model
  class(toret$c.model)<-"coxph"
  toret$ncases<-sum(extra.only$outcome)
  toret$ncasedays<-ncases
  toret$ncontroldays<-ncontrols
  class(toret)<-'casecross'
  return(toret)
}
