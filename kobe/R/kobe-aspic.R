########################################################################################
#### aspic stuff for Kobe ##############################################################
########################################################################################

#### exported function
kobeAspic=function(object,prb=NULL,dir="",
                   what=c("sims","trks","pts","smry","wrms")[1],
                   prob=c(0.75,0.5,.25),year=NULL,nwrms=10){
    
    if (dir!=""){
      bio=paste(dir,object,sep="/")
      prb=paste(dir,prb,   sep="/")}
    else bio=object
    
    if (length(prb) %in% 0:1)
       res=ioAspic(bio,prb,what=what,prob=prob,nwrms=nwrms,year=year)
    
    if (length(prb)>1){
       res=mlply(prb, function(x,bio,prob=prob,nwrms=nwrms,what=what)
                                   ioAspic(bio=bio,prb=x,prob=prob,nwrms=nwrms,what=what),
                      bio=bio,prob=prob,nwrms=nwrms,what=what)
                 
      res=list(trks=ldply(res, function(x) x$trks),
               pts =ldply(res, function(x) x$pts),
               smry=ldply(res, function(x) x$smry),
               wrms=ldply(res, function(x) x$wrms),
               sims=ldply(res, function(x) x$sims))
      }
    
    if (length(what)==1)
       return(res[[what]])
    else
       return(res[what]) }

#### read ASPIC bootstrapped assessment file 
aspicBio=function(file){
  
  rtn=NULL
  for (iFile in file){
    t.  <-scan(iFile,skip=4)
    nits<-scan(iFile,skip=1,nmax=1)
    yrs <-scan(iFile,skip=2,nmax=2)
    nyrs<-diff(yrs)
    nval<-nyrs*2+3
    
    yrs <-yrs[1]:yrs[2]
    
    b.  <-data.frame(stock  =t.[unlist(tapply(((1:nits)-1)*nval+2,     1:nits,function(x,y=nyrs+1)  x:(x+y-1)))],year=yrs,               iter=rep(1:nits,each=length(yrs)))
    f.  <-data.frame(harvest=t.[unlist(tapply(((1:nits)-1)*nval+nyrs+4,1:nits,function(x,y=nyrs)    x:(x+y-1)))],year=yrs[-length(yrs)], iter=rep(1:nits,each=(length(yrs)-1)))
    
    bmsy<-data.frame(bmsy=t.[unlist(tapply(((1:nits)-1)*nval+1,     1:nits,function(x,y=1)      x:(x+y-1)))],iter=1:nits)
    fmsy<-data.frame(fmsy=t.[unlist(tapply(((1:nits)-1)*nval+nyrs+3,1:nits,function(x,y=1)      x:(x+y-1)))],iter=1:nits)
    
    sh   =merge(b.,  f.,  by=c("iter","year"))
    rf   =merge(bmsy,fmsy,by=c("iter"))
    res  =merge(sh,rf, by=c("iter"))
    
    res$stock  =res$stock/res$bmsy
    res$harvest=res$harvest/res$fmsy
   
    rtn=rbind(rtn,cbind(run=iFile,res))
    }

  return(rtn[do.call("order",rtn[,c("run","iter","year")]),c("run","iter","year","stock","harvest","bmsy","fmsy")])}

#### reads ASPIC projection file
aspicPrb=function(file){
  ## Stuff
  nits<-scan(file,skip=1,nmax=1)
  yrs <-scan(file,skip=2,nmax=2)
  t.  <-scan(file,skip=4)
  ncol<-yrs[2]-yrs[1]+2
  
  ## stock
  first<-rep((1:nits-1)*ncol*2,each=yrs[2]-yrs[1]+1)+(1:(ncol-1))+1
  b.   <-data.frame("stock"=t.[first],year=yrs[1]:yrs[2],iter=rep(1:nits,each=length(yrs[1]:yrs[2])))
  
  first<-((1:nits-1)*ncol*2)+1
  bmsy <-data.frame(bmsy=t.[first],iter=1:nits)
  b.   <-merge(b.,bmsy,by="iter")
  
  ## F
  first<-rep((1:nits-1)*ncol*2+ncol,each=yrs[2]-yrs[1]+1)+(1:(ncol-1))+1
  f.   <-data.frame("harvest"=t.[first],year=yrs[1]:yrs[2],iter=rep(1:nits,each=length(yrs[1]:yrs[2])))
  
  first<-((1:nits-1)*ncol*2)+ncol+1
  fmsy <-data.frame(fmsy=t.[first],iter=1:nits)
  f.   <-merge(f.,fmsy,by="iter")
  
  res=merge(b.,f.,by=c("iter","year"))
  
  res$stock  =res$stock/res$bmsy
  res$harvest=res$harvest/res$fmsy
  
  return(res[do.call("order",res[,c("iter","year")]),c("iter","year","stock","harvest","bmsy","fmsy")])}


## Heavy lifting functions 
ioAspic=function(bio,prb,prob=c(0.75,0.5,.25),
                 what=c("sims","trks","pts","smry","wrms")[1],
                 year=NULL,nwrms=10){
    
    if (!all(what %in% c("trks","pts","smry","wrms","sims"))) stop("what not in valid options")
    
    if (all(tolower(getExt(bio)) %in% "bio")) 
       bio.=aspicBio(bio) else stop("Arg not a .bio file")
    
    if (!is.null(prb)){
      if (tolower(getExt(prb)) %in% "prb") 
          prb.=aspicPrb(prb) else stop("First  arg not a .prb file")
      res=rbind(bio.,prb.)
      res=res[!is.na(res$stock) & !is.na(res$harvest),]
      }
    else
      res=bio.
    
    if (is.null(year)) pts=max(bio.$year)
      
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
