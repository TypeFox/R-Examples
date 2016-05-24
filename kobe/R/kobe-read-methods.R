utils::globalVariables(c("sims"))

#' @importFrom reshape melt cast
#' @importFrom plyr mlply llply

if (!isGeneric('kobe'))
  setGeneric('kobe', function(object,method,...) standardGeneric('kobe'))

setMethod('kobe',  signature(object='character',method="character"), 
          function(object,
                   method=c("aspic","adapt","bsp","mfcl","ss","sam","vpa"),
                   what  =c("sims","trks","pts","smry","wrms")[1],
                   dir   ="",
                   prob  =c(0.75,0.5,0.25),
                   year  =NULL,
                   nwrms =10,...) {
   
    method=tolower(method)
    if (any("2box" == method)) method["2box" == method]="adapt"   
    switch(substr(method[1],1,2),
           ad=kobe2box( object,what=what,prob=prob,year=year,nwrms=nwrms,...),
           as=kobeAspic(object,dir=dir,what=what,prob=prob,year=year,nwrms=nwrms,...),
           bs=kobeBsp(  object,dir=dir,what=what,prob=prob,year=year,nwrms=nwrms,...),
           mf=kobeMFCL( object,dir=dir,what=what,prob=prob,year=year,nwrms=nwrms,...),
           ss=kobeSS(   object,what=what,prob=prob,year=year,nwrms=nwrms,...))
    })


setMethod('kobe',  signature(object="data.frame",method="missing"),  
          function(object,method,what=c("sims","trks","pts","smry","wrms")[1],dir="",
                   prob=c(0.75,0.5,.25),year=NULL,nwrms=10){
  
  res=llply(object, function(x,what=what,prob=prob,year=year,nwrms=nwrms)
    kobeFn(object,what=what,prob=prob,year=year,nwrms=nwrms),
            what=what,prob=prob,year=year,nwrms=nwrms)
  
  res=list(trks=ldply(res, function(x) x$trks),
           pts =ldply(res, function(x) x$pts),
           smry=ldply(res, function(x) x$smry),
           wrms=ldply(res, function(x) x$wrms),
           sims=ldply(res, function(x) x$sims))
  
  if (length(what)==1)
    return(res[[what]])
  else
    return(res[what]) })

