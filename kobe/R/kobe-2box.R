# #######################################################################################
# ### SS stuff for Kobe #################################################################
# #######################################################################################

##note need to think about years
utils::globalVariables(c("stock","harvest","year","tac"))

nmsRef <- c("iter", 
            "fmsy",   "ymsy",   "yrmsy",    "srmsy",   "sprmsy",  "ssbmsy",
            "fmax",             "yrmax",    "srmax",   "sprmax",  "ssbmax",
            "f0.1",             "yr0.1",    "sr0.1",   "spr0.1",  "ssb0.1",
            "f20",              "yr20",     "sr20",               "ssb20",
            "f30",              "yr30",     "sr30",               "ssb30",
            "f40",              "yr40",     "sr40",               "ssb40",
            "f90max", "y90max", "yr90max",  "sr90max",            "ssb90max",
            "f75max", "y75max", "yr75max",  "sr75max",            "ssb75max")

 
#setMethod('kobe2box', signature(object='character'),
kobe2box=function(object,proxy=c("fmsy","fmax","f0.1","f20","f30","f40","f90max","f75max")[3], 
                         what=c("sims","trks","pts","smry","wrms")[1],
                         prob=c(0.75,0.5,0.25),year=NULL,nwrms=10){
            
    if (length(object)==1)
       res=io2box(object,proxy=proxy,what=what,prob=prob,nwrms=nwrms,year=year)
            
    if (length(object)>1){
       res=mlply(object, function(x,proxy=proxy,what=what,prob=prob,nwrms=nwrms,year=year)
                           io2box(x,proxy=proxy,what=what,prob=prob,nwrms=nwrms,year=year),
                                    proxy=proxy,what=what,prob=prob,nwrms=nwrms,year=year)
              
              res=list(trks=ldply(res, function(x) x$trks),
                       pts =ldply(res, function(x) x$pts),
                       smry=ldply(res, function(x) x$smry),
                       wrms=ldply(res, function(x) x$wrms),
                       sims=ldply(res, function(x) x$sims))
       }
    
    if (length(what)==1) 
      return(res[[what]]) 
    else 
      return(res[what])}


## Heavy lifting functions ##############################################################
readKobe2box=function(dir,proxy=c("fmsy","fmax","f0.1","f20","f30","f40","f90max","f75max")[3]){

  rfp=read.table(paste(dir,"BENCH-1.OUT",sep="/"),header=F,skip=1,col.names=nmsRef)[,c("iter",proxy,paste("ssb",substr(proxy,2,nchar(proxy)),sep=""))]
 
  bio=read.table(paste(dir,"SSBIO-1.OUT",sep="/"),header=F,skip=0)
  names(bio)=c("tac","iter",seq(dim(bio)[2]-2))
  bio=merge(bio,rfp)
  bio=melt(bio,id.vars=c("tac",names(rfp)),variable_name="year")
  names(bio)[6]="stock"
  
  hvt=read.table(paste(dir,"Fapex-1.OUT",sep="/"),header=F,skip=0)
  names(hvt)=c("tac","iter",seq(dim(hvt)[2]-2))
  hvt=merge(hvt,rfp)
  hvt=melt(hvt,id.vars=c("tac",names(rfp)),variable_name="year")
  names(hvt)[6]="harvest"
  
  res=merge(bio,hvt[,c("tac","iter","year","harvest")])
  res=transform(res,stock  =stock/res[,paste("ssb",substr(proxy,2,nchar(proxy)),sep="")],
                    harvest=harvest/res[,proxy])
  
  res=res[do.call(order,res[,c("iter","tac","year")]),]
  
  dimnames(res)[[1]]=sort(as.numeric(dimnames(res)[[1]]))
    
  return(res)}

io2box=function(x,proxy=c("fmsy","fmax","f0.1","f20","f30","f40","f90max","f75max")[3],
                prob=c(0.75,0.5,.25),what=c("sims","trks","pts","smry","wrms")[1],nwrms=10,year=NULL){
  
  if (!all(what %in% c("trks","pts","smry","wrms","sims"))) stop("what not in valid options")
  
  res=readKobe2box(x,proxy)

  sims.=NULL
  trks.=NULL
  pts. =NULL
  smry.=NULL
  wrms.=NULL
   
  if ("sims" %in% what)
      sims.=res

   if ("trks" %in% what){ 
       ssb =ddply(res,.(year,tac),function(x) quantile(x$ssb,    prob))
       hvt =ddply(res,.(year,tac),function(x) quantile(x$harvest,prob))
       trks.=data.frame(melt(ssb,id.vars=c("year","tac")),harvest=melt(hvt,id.vars=c("year","tac"))[,4])
       names(trks.)[3:4]=c("Percentile","ssb")}

   if ("pts" %in% what & !is.null(year))
        pts.=res[res$year %in% year,]
           
   if ("smry" %in% what)
       smry.   =ddply(res,  .(year), function(x) data.frame(ssb        =median(x$ssb,       na.rm=TRUE),
                                                            harvest    =median(x$harvest,   na.rm=TRUE),
                                                            red        =mean(x$red,         na.rm=TRUE),
                                                            yellow     =mean(x$yellow,      na.rm=TRUE),
                                                            green      =mean(x$green,       na.rm=TRUE),
                                                            overFished =mean(x$overFished,  na.rm=TRUE),
                                                            overFishing=mean(x$overFishing, na.rm=TRUE)))
    
   if ("wrms" %in% what)
       wrms.=res[res$iter %in% sample(unique(res$iter),nwrms),c("iter","year","ssb","harvest")]
 
  
  return(list(trks=trks.,pts=pts.,smry=smry.,wrms=wrms.,sims=sims.))}

