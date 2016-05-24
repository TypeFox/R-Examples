## general interface for finding probabilities

## rpel = function(string, envir=.GlobalEnv) {
##   eval(parse(text=string), envir=envir)
## }


####
dpqrfuncs = function(type=c("p","d","q","r"), container=NULL) {

  randomFamilies = list(
    "norm" = c("mean","sd",0,1),
    "t" = c("df","ncp","",0),
    "f" = c("df1","df2","",""),
    "chisq" = c("df","ncp","",0),
    "unif" = c("min","max",0,1),
    "exp" =  c("rate",1),
    "weibull" = c("shape","scale","",1),
    "gamma" = c("shape","rate","",1),
    "beta" = c("shape1","shape2","",""),
    "cauchy" = c("location","scale","",""),
    "logis" = c("location","scale",0,1),
    "lnorm" = c("meanlog","sdlog",0,1),
    "pois" = c("lambda","",1,""),
    "binom"=c("size","prob",1,.5),
    "geom" = c("prob",""),
    "nbinom" = c("size","prob","","")
    )
  helpPage = c(
    "norm" = "Normal",
    "t" = "TDist",
    "f" = "FDist",
    "chisq" = "Chisquare",
    "unif" = "Uniform",
    "exp" =  "Exponential",
    "weibull" = "Weibull",
    "gamma" = "Gamma",
    "beta" = "Beta",
    "cauchy" = "Cauchy",
    "logis" = "Logistic",
    "lnorm" = "Lognormal",
    "pois" = "Poisson",
    "binom"= "Binomial",
    "geom" = "Geometric",
    "nbinom" = "NegBinomial"
    )
  
  theFirstArg = c("p"="q","d"="x","q"="p","r"="n")
  doWhat = c("p"="Probability","d"="p.d.f","q"="Quantile","r"="Random sample")
  
  
  type = match.arg(type)
    
  gp = ggroup(horizontal=FALSE, container=container)
  glabel(paste("Find a",doWhat[type]), container=gp)

  ## widgets
  distSelector = gdroplist(paste(type,names(randomFamilies),sep=""))
  sampleSize = gedit(1) #, coerce.with=rpel)


  param1label = glabel("mean")
  param1 = gedit("0") #, coerce.with=rpel)
  param2label = glabel("sd")
  param2 = gedit("1") #, coerce.with=rpel)
  doLog = gcheckbox("")
  doLowerTail = gcheckbox("",checked=TRUE)
  doLogP = gcheckbox("")
  
  saveAs = gedit("")

  output = gtext("", font.attr=c(style="monospace"))

  ## layout
  tbl = glayout(container=gp)
  tbl[1,1] = glabel("Choose a distribution")
  tbl[1,2] = distSelector
  tbl[2,1] = glabel(theFirstArg[type])
  tbl[2,2] = sampleSize
  tbl[3,1] = param1label
  tbl[3,2] = param1
  tbl[4,1] = param2label
  tbl[4,2] = param2
  
  i = 5
  if(type == "d") {
    tbl[5,1] = glabel("log"); tbl[5,2] = doLog
    i = 6
  } else if (type %in% c("p","q")) {
    tbl[5,1] = glabel("lower.tail");tbl[5,2] = doLowerTail
    tbl[6,1] = glabel("log.p"); tbl[6,2] = doLogP
    i = 7
  }
  
  tbl[i,1] = glabel("Save output as:")
  tbl[i,2] = saveAs

  bgp = ggroup(container=gp)
  tbl[i+1,2] <- bgp
  addSpring(bgp)
  findSample = gbutton("ok",container=bgp)
  helpButton = gbutton("help", container=bgp)

  visible(tbl) <- TRUE

  ## replace, using pmg.cli now
  ## add(gp, output, expand = TRUE)  

  ### actions
  ## select distr, update parameters
  addhandlerchanged(distSelector, handler = function(h,...) {

    theFunc = svalue(distSelector)
    theDist = substr(theFunc,2,stop=nchar(theFunc))
    theParams = randomFamilies[[theDist]]

    ## enabled?
    if(length(theParams) == 2) {
      svalue(param2) <- ""; enabled(param2) <- FALSE
    } else {
      ## ensure they are on
      enabled(param2label) <- TRUE
      enabled(param2) <- TRUE
    }

    ## param1
    svalue(param1label) <- theParams[1]
    svalue(param1) <- theParams[3]
    if(length(theParams) == 2) {
      svalue(param2label) <- ""
      svalue(param2) <- ""
    } else {
      svalue(param2label) <- theParams[2]
      svalue(param2) <- theParams[4]
    }

    ## clear output
    svalue(output) <- ""
  })
    

  ## click OK
  addhandlerchanged(findSample, handler = function(h,...) {

    theFunc = svalue(distSelector)
    theDist = substr(theFunc,2,stop=nchar(theFunc))
    theParams = randomFamilies[[theDist]]


    ## We paste together a command to use with pmg.cli, not a do.
##     theArgs = list()
##     theArgs[[theFirstArg[type]]]=svalue(sampleSize)
##     theArgs[[theParams[1]]] <- svalue(param1)
##     if(length(theParams) > 2)
##       theArgs[[theParams[2]]] <- svalue(param2)
##     if(type == "d") {
##       if(svalue(doLog))
##         theArgs[["log"]]=TRUE
##     } else if (type %in% c("p","q")) {
##       if(svalue(doLowerTail) == FALSE)
##         theArgs[["lower.tail"]] = FALSE
##       if(svalue(doLogP))
##         theArgs[["log.p"]] = TRUE
##     }

##     res = do.call(theFunc, theArgs)

##     ## save 
##     if(svalue(saveAs) != "")
##       assign(svalue(saveAs),res, envir=.GlobalEnv)
##     ## update output
##     oldWidth = getOption("width"); options("width"=60)
##     svalue(output) <- capture.output(print(res))
##     options(width=oldWidth)

    cmd = paste(theFunc,"(",sep="")
    cmd = paste(cmd, theFirstArg[type], "=", svalue(sampleSize), sep="")
    cmd = paste(cmd, ",", theParams[1], "=",  svalue(param1), sep="")
    if(length(theParams) >2)
      cmd = paste(cmd, ",", theParams[2], "=",  svalue(param2), sep="")
    if(type == "d") {
      if(svalue(doLog))
        cmd = paste(cmd, ", log=TRUE", sep="")
    }
    if(type %in% c("p","q")) {
      if(svalue(doLowerTail) == FALSE)
        cmd = paste(cmd, ", lower.tail=FALSE", sep="")
      if(svalue(doLogP))
        cmd = paste(cmd, ", log.p=TRUE", sep="")
    }
    cmd = paste(cmd, ")", sep="")

    ## do we save as something?
    if(svalue(saveAs) != "")
      names(cmd) <- make.names(svalue(saveAs))

    svalue(pmg.cli) <- cmd
      
      

    
  })
  ## help
  addhandlerchanged(helpButton, handler = function(h,...) {
    theFunc = svalue(distSelector)
    theDist = substr(theFunc,2,stop=nchar(theFunc))    
    ghelp(helpPage[theDist], container=pmgWC$new(paste("Help on",theFunc)))
  })

  ##
  return(gp)
}

