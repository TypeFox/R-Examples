utils::globalVariables(c("variable","value"))

########################################################################################
#### aspic stuff for Kobe ##############################################################
########################################################################################

#### exported function
kobeBsp=function(f,b,dir="",what=c("sims","trks","pts","smry","wrms")[1],
                 prob=c(0.75,0.5,.25),year=NULL,nwrms=10){
  
  res=ioBsp(f,b,what=what,prob=prob,nwrms=nwrms,year=year)
  
#   if (length(prb)>1){
#     res=mlply(prb, function(x,bio,prob=prob,nwrms=nwrms,what=what)
#       ioAspic(bio=bio,prb=x,prob=prob,nwrms=nwrms,what=what),
#               bio=bio,prob=prob,nwrms=nwrms,what=what)
#     
#     res=list(trks=ldply(res, function(x) x$trks),
#              pts =ldply(res, function(x) x$pts),
#              smry=ldply(res, function(x) x$smry),
#              wrms=ldply(res, function(x) x$wrms),
#              sims=ldply(res, function(x) x$sims))
#  }
  
  if (length(what)==1)
    return(res[[what]])
  else
    return(res[what]) }

## Heavy lifting functions 
ioBsp=function(f,b,prob=c(0.75,0.5,.25),
               what=c("sims","trks","pts","smry","wrms")[1],
               year=NULL,nwrms=10,firstyear=0){
  
  if (!all(what %in% c("trks","pts","smry","wrms","sims"))) stop("what not in valid options")
  
  res=getBSP(f,b,firstyear)
  
  if (is.null(year)) pts=max(res$year)
  
  trks. =NULL
  pts.  =NULL
  smry. =NULL
  wrms. =NULL
  sims. =NULL
  
  if ("trks" %in% what){ 
    stock  =ddply(res,.(year),function(x) quantile(x$stock,    prob, na.rm=TRUE))
    harvest=ddply(res,.(year),function(x) quantile(x$harvest,  prob, na.rm=TRUE))
    trks.=data.frame(melt(stock,id.vars="year"),"harvest"=melt(harvest,id.vars="year")[,3])
    names(trks.)[c(2,3)]=c("Percentile","stock")}
  
  if ("pts" %in% what)
    pts.=res[res$year %in% pts,]
  
  if ("sims" %in% what)
    sims.=res
  
  if ("smry" %in% what)
    smry. =ddply(data.frame(res,kobeP(res$stock,res$harvest)),
                 .(year), function(x) data.frame(stock      =median(x$stock,       na.rm=TRUE),
                                                 harvest    =median(x$harvest,     na.rm=TRUE),
                                                 red        =mean(  x$red,         na.rm=TRUE),
                                                 yellow     =mean(  x$yellow,      na.rm=TRUE),
                                                 green      =mean(  x$green,       na.rm=TRUE),
                                                 overFished =mean(  x$overFished,  na.rm=TRUE),
                                                 overFishing=mean(  x$overFishing, na.rm=TRUE)))
  
  if ("wrms" %in% what)
    wrms.=res[res$iter %in% sample(unique(res$iter),nwrms),c("iter","year","stock","harvest")]
  
  res=list(trks=trks.,pts=pts.,smry=smry.,wrms=wrms.,sims=sims.)
  
      #if (length(what)==1) return(res[[what]])
  
  return(list(trks=trks.,pts=pts.,smry=smry.,wrms=wrms.,sims=sims.))}

#fileF="/home/laurie/Desktop/Dropbox/collaboration/alb-wg/Inputs/BSP/for Kobe plots/kobeFalb2013SAf11.csv"
#fileB="/home/laurie/Desktop/Dropbox/collaboration/alb-wg/Inputs/BSP/for Kobe plots/kobeBalb2013SAf11.csv"
#getBSP(fileF,fileB,1928)
  
getBSP=function(fileF,fileB,firstyear=0){
  resF=transform(melt(read.csv(fileF,header=F)),year=as.numeric(variable)+firstyear,
                        harvest=value)[,-(1:2)]
  resB=transform(melt(read.csv(fileB,header=F)),year=as.numeric(variable)+firstyear,
                        stock=value)[,-(1:2)]
  res=cbind(resF,stock=resB$stock)
  return(res)}

