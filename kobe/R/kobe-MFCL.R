########################################################################################
#### SS stuff for Kobe #################################################################
########################################################################################
  
utils::globalVariables(c("bio"))

#setMethod('kobeMFCL', signature(object='character'),
kobeMFCL=function(object,dir="",what=c("sims","trks","pts","smry","wrms")[1],
                                prob=c(0.75,0.5,0.25),year=NULL,nwrms=10){
                     
            if (length(object)==1)
              res=ioMFCL(object,prob=prob,nwrms=nwrms)
            
            if (length(object)>1){
              res=mlply(object, function(x,prob=prob,nwrms=nwrms,what=what)
                ioMFCL(x,what=what,prob=prob,year=year,nwrms=nwrms),
                        bio=bio,what=what,year=year,prob=prob,nwrms=nwrms)
              
              res=list(trks=ldply(res, function(x) x$trks),
                       pts =ldply(res, function(x) x$pts),
                       smry=ldply(res, function(x) x$smry),
                       wrms=ldply(res, function(x) x$wrms),
                       sims=ldply(res, function(x) x$sims))
            }
            
            if (length(what)==1) return(res[[what]]) else return(res[what])}

getplotdat1<-function (h = "", plotrepfile, skip = 1) {
  dat <- readLines(plotrepfile)
  recnum <- grep(h, dat)
  scanText(dat[recnum + skip], what = 0)}

scanText<-function(string, what = character(0), ...){
  ## Like scan() but reading from a vector of character strings
  tc <- textConnection(string)
  result <- scan(tc, what = what, quiet = TRUE, ...)
  close(tc)
  return(result)}

getB2Bmsy=
  function(plotrepfile="plot.rep"){
    ##==============================================================
    ## Biomass at MSY
    ##==============================================================
    getplotdat1(plotrepfile,h="# Total biomass over total biomass at MSY")
  }

getF2Fmsy =
  function(plotrepfile="plot.rep"){
    ##==============================================================
    ## Biomass at MSY
    ##==============================================================
    getplotdat1(plotrepfile,h="# Aggregate F over F at MSY")
  }

getSB2SBmsy=
  function(plotrepfile="plot.rep"){
    ##==============================================================
    ## Biomass at MSY
    ##==============================================================
    getplotdat1(plotrepfile,h="# Adult biomass over adult biomass at MSY")
  }

getyrs=
  function(plotrepfile="plot.rep"){
    #==============================================================
    # Vector of years
    #==============================================================
    yr <- range(floor(unlist(getrtimes(plotrepfile))))
    seq(yr[1],yr[2])
  }

getrtimes<-
  function(plotrepfile="plot.rep"){
    ##==============================================================
    ## Time of each realization by fishery (down)
    ##==============================================================
    dat <- getplotdat4("# Time of each realization by fishery",plotrepfile)
    nreal <- getnreal(plotrepfile)
    nfish <- getnfish(plotrepfile)
    splitter <- rep(seq(nfish), nreal)
    split(dat, splitter)
  }

getnfish <-
  function(plotrepfile="plot.rep"){
    ##==============================================================
    ## Number of fisheries
    ##==============================================================
    getplotdat1(plotrepfile,h="# Number of fisheries")
  }

getnreal <-
  function(plotrepfile="plot.rep"){
    ###==============================================================
    ### Number of realizations per fishery
    ###==============================================================
    getplotdat1(plotrepfile,h="# Number of realizations per fishery")
  }

getplotdat4 <- function(h="",plotrepfile) {
  ##=================================================
  ## Start listing after header h.  Quit if encounter
  ##  "^#"
  ##=================================================
  dat <- readLines(plotrepfile)
  rec1 <- grep(h, dat)
  if(length(rec1) <= 0)
    stop(paste('"',h,'"',"not found in",plotrepfile," Die yuppie scum!"))
  recnum <- rec1+1
  tt <- numeric(0)
  for(i in recnum:length(dat)) {
    if (regexpr("^#", dat[i]) != -1) break
    tt <- c(tt, scanText(dat[i], what = 0))
  }
  tt
}


object="/home/laurie/Desktop/Dropbox/collaboration/Shelton/ALBN/4B/plot-09.par.rep"
  
## Heavy lifting functions ##############################################################
ioMFCL=function(object,what=c("sims","trks","pts","smry","wrms")[1],prob=c(0.75,0.5,0.25),year=NULL,nwrms=10){

  if (!all(what %in% c("trks","pts","smry","wrms","sims"))) stop("what not in valid options")
  
  res=data.frame(year=getyrs(object),biomass=getB2Bmsy(object),harvest=getF2Fmsy(object),stock=getSB2SBmsy(object),ssb=getSB2SBmsy(object))
  
  trks. =NULL
  pts.  =NULL
  smry. =NULL
  wrms. =NULL
  sims. =NULL
  
  if ("trks" %in% what){ 
    stock  =ddply(res,.(year),function(x) quantile(x$stock,    prob, na.rm=TRUE))
    harvest=ddply(res,.(year),function(x) quantile(x$harvest,  prob, na.rm=TRUE))
    trks.=data.frame(melt(stock,id.vars="year"),harvest=melt(harvest,id.vars="year")[,3])
    names(trks.)[c(2,3)]=c("Percentile","stock")}
  
  if ("pts" %in% what & !is.null(year))
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
  
  #if ("wrms" %in% what)
  #  wrms.=res[res$iter %in% sample(unique(res$iter),nwrms),c("iter","year","ssb","harvest")]
  
  res=list(trks=trks.,pts=pts.,smry=smry.,wrms=wrms.,sims=sims.)
  
  return(res)}  

