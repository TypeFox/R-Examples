## functions by Yvonnick Noel yvonnick.noel@uhb.fr

## TODO
## * histogramAndDensity -- get density drawn right
## probabilityCalculator -- get Type of calculation spelled correctly
## * speed up the drawing
## constructionOFNormal -- add handler to no. of variables., ... m and n
histogramAndDensity = function(container=gwindow("Window 1: idea of a density", visible=FALSE), ...) {

  availDists = c(Uniform = "unif", Normal = "norm", Gamma = "gamma")
  theParams = list(
    "norm" = c("mean","sd",0,1),
    "t" = c("df","ncp","",0),
    "chisq" = c("df","ncp","",0),
    "f" = c("df1","df2","",""),
    "binom"=c("size","prob",1,.5),
    "pois" = c("lambda","",1,""),
    "gamma" = c("shape","rate","",1),
    "beta" = c("shape1","shape2","",""),
    "unif" = c("min","max",0,1)
    )
  
  
  updatePlot = function(h, ...) {
    ### makeEmptyPlot - windows issue
    plot.new()
    plot.window(xlim=c(0,1),ylim=c(0,1))

    
    
    rfunc = paste("r",availDists[svalue(distribution)],sep="")
    dfunc = paste("d",availDists[svalue(distribution)],sep="")
    x <- do.call(rfunc, list(svalue(sampleSize),svalue(param1),svalue(param2)))
    breaks <- svalue(cutpoints)
    if(is.null(breaks))
        breaks <- "sturges"

    hh=hist(x,breaks=breaks,plot=FALSE)
    ## now plot histogram
    if(svalue(displayWhat) == "Counts") {
      plot(hh,
           main = paste("Distribution",svalue(distribution)),
           ylab=svalue(displayWhat)
           ) 
    } else if(svalue(displayWhat) == "Frequencies") {
      hh$counts = hh$counts / svalue(sampleSize)
      plot(hh,
           main = paste("Distribution",svalue(distribution)),
           ylab=svalue(displayWhat)
           ) 
    } else {
      hist(x,probability=TRUE,
           main = paste("Distribution",svalue(distribution)),
           ylab=svalue(displayWhat)
           ) 
      addDensity()
    }
  }

  addDensity = function(...) {
    dfunc = paste("d",availDists[svalue(distribution)],sep="")
    dFunc = get(dfunc)
    f <- function(x) dFunc(x, svalue(param1), svalue(param2))
    curve(f, lwd=2, col="red",add=T)
  }

  ## the distribution used
  distribution = gdroplist(names(availDists),horizontal=FALSE,handler=updatePlot)
  addhandlerchanged(distribution, handler = function(h,...) {
    theDist = availDists[svalue(h$obj)]
    svalue(param1label) = theParams[[theDist]][1]
    svalue(param1) = theParams[[theDist]][3]
    
    if(theParams[[theDist]][2] == "") {
      svalue(param2label) = theParams[[theDist]][2]
      svalue(param2) = ""      
      enabled(param2) <- FALSE
    } else {
      enabled(param2) <- TRUE
      svalue(param2label) = theParams[[theDist]][2]
      svalue(param2) = theParams[[theDist]][4]
    }
  })
  


  sampleSize = gradio(c(500, 5000, 50000),handler=updatePlot)
  displayWhat = gradio(c("Counts","Frequencies","Density of frequencies"),handler=updatePlot)
  displayFunc = gcheckbox("Theoretical density",handler=addDensity)

  param1label = glabel(theParams[[availDists[1]]][1])
  param2label = glabel(theParams[[availDists[1]]][2])
  param1 = gedit(theParams[[availDists[1]]][3],width=5,coerce.with=rpel)
  param2 = gedit(theParams[[availDists[1]]][4],width=5,coerce.with=rpel)

  cutpoints = gedit("",coerce.with=rpel)


  BigGroup = ggroup(container= container, ...)
  group = ggroup(horizontal = FALSE, container = BigGroup)
  
  
  tmp = gframe("Distribution", container = group)
  distribGroup = glayout(container=tmp)
  distribGroup[1,1]=glabel("Law")
  distribGroup[1,2]=distribution
  distribGroup[2,1]=param1label
  distribGroup[2,2]=param1
  distribGroup[3,1]=param2label
  distribGroup[3,2]=param2
  visible(distribGroup)=TRUE
  
  tmp = gframe("Sample size", container = group)
  add(tmp, sampleSize)
  
  tmp = gframe("Display", container = group)
  add(tmp,displayWhat)
                                        #add(tmp,displayFunc)
  
  tmp = gframe("Cutpoints", container = group)
  add(tmp,cutpoints,expand=TRUE)
  
  addSpring(group)
  
  buttonGroup=ggroup(container=group)
  
  if(missing(container))
    gbutton("cancel", container=buttonGroup, handler = function(h,...) dispose(container))
  
  addSpring(buttonGroup)
  gbutton("display",container=buttonGroup, handler=updatePlot)
  
  add(BigGroup, ggraphics())

  visible(container) <- TRUE
  
  invisible(BigGroup)
  
}



## probability calculator
probabilityCalculator = function(container=gwindow("Probability caculator")) {

    availDists = c(Normal="norm",Student="t","Chi-2"="chisq",
        Fisher="f",Binomial="binom",Poisson="pois",Gamma="gamma",Beta="beta")

    theParams = list(
        "norm" = c("mean","sd",0,1),
        "t" = c("df","ncp","",0),
        "chisq" = c("df","ncp","",0),
        "f" = c("df1","df2","",""),
        "binom"=c("size","prob",1,.5),
        "pois" = c("lambda","",1,""),
    "gamma" = c("shape","rate","",1),
        "beta" = c("shape1","shape2","",""),
        "unif" = c("min","max",0,1)
        )
    
    initOptions = function(h, ...) {
        
        r2s.distrib = svalue(distribution)
        r2s.is1P = r2s.distrib %in% c("Student","Chi-2","Poisson")
        if(r2s.is1P) svalue(param2)=""
        svalue(result)=""
        svalue(value)=""
        
    }

    updatePlot = function(h, ...) {
  
        r2s.distrib = svalue(distribution)
        r2s.param1 = svalue(param1)
        r2s.param2 = svalue(param2)
        r2s.value = svalue(value)
        

        r2s.p = svalue(calcWhat) == "Find quantile"
        r2s.right = svalue(side)=="to right"

        r2s.isDiscrete = r2s.distrib %in% c("Binomial","Poisson")
        r2s.is1P = r2s.distrib %in% c("Student","Chi-2","Poisson")
        r2s.is01 = function(x) (x>=0)&&(x<=1)
        r2s.isInteger = function(x) abs(x)==round(x)
        r2s.probf = availDists

        r2s.dfunction = eval(parse(text=paste("d",r2s.probf[r2s.distrib],sep="")))
        r2s.pfunction = eval(parse(text=paste("p",r2s.probf[r2s.distrib],sep="")))
        r2s.qfunction = eval(parse(text=paste("q",r2s.probf[r2s.distrib],sep="")))
        r2s.rfunction = eval(parse(text=paste("r",r2s.probf[r2s.distrib],sep="")))

        ## Chosen distribution has two parameters
        if(!r2s.is1P) {
            ## Check parameter values
            if(r2s.distrib=="Binomial") {
                stopifnot(r2s.isInteger(r2s.param1) && r2s.is01(r2s.param2)) }
            if(r2s.distrib=="Fisher") {
                stopifnot(r2s.isInteger(r2s.param1) && r2s.isInteger(r2s.param2))
            }

            if(r2s.p) {
                r2s.prob = r2s.value
                if(!r2s.isDiscrete) {
                    if(r2s.right) r2s.value = r2s.qfunction(1-r2s.prob,r2s.param1,r2s.param2)
                    else          r2s.value = r2s.qfunction(r2s.prob,r2s.param1,r2s.param2)
                } else {
                    if(r2s.right) r2s.value = r2s.qfunction(1-r2s.prob,r2s.param1,r2s.param2)
                    else          r2s.value = r2s.qfunction(r2s.prob,r2s.param1,r2s.param2)
                }} else {
                    if(!r2s.isDiscrete) {
                        if(r2s.right) { r2s.prob = 1-r2s.pfunction(r2s.value,r2s.param1,r2s.param2) } else {
                            r2s.prob = r2s.pfunction(r2s.value,r2s.param1,r2s.param2)}} else {
                                if(r2s.right) { r2s.prob = 1-r2s.pfunction(r2s.value-1,r2s.param1,r2s.param2) } else {
                                    r2s.prob = r2s.pfunction(r2s.value,r2s.param1,r2s.param2)}}}
            
            r2s.dens = r2s.dfunction(r2s.value,r2s.param1,r2s.param2)
                                        # Chosen distribution has only one parameter
        } else {
            
            svalue(param2)=""
            
            if(r2s.distrib=="Student") {
                stopifnot(r2s.isInteger(r2s.param1)) }
            if(r2s.distrib=="Chi-2") {
                stopifnot(r2s.isInteger(r2s.param1)) }
            
            if(r2s.p) {
                r2s.prob = r2s.value
                if(!r2s.isDiscrete) {
                    if(r2s.right) r2s.value = r2s.qfunction(1-r2s.prob,r2s.param1)
                    else          r2s.value = r2s.qfunction(r2s.prob,r2s.param1)
     } else {
         if(r2s.right) r2s.value = r2s.qfunction(1-r2s.prob,r2s.param1)
         else          r2s.value = r2s.qfunction(r2s.prob,r2s.param1)
     }} else {
       if(!r2s.isDiscrete) {
           if(r2s.right) { r2s.prob = 1-r2s.pfunction(r2s.value,r2s.param1) } else {
                         r2s.prob = r2s.pfunction(r2s.value,r2s.param1)}} else {
         if(r2s.right) { r2s.prob = 1-r2s.pfunction(r2s.value-1,r2s.param1) } else {
                         r2s.prob = r2s.pfunction(r2s.value,r2s.param1)}}}

   r2s.dens = r2s.dfunction(r2s.value,r2s.param1)
 }

 # Result
 svalue(result)=ifelse(r2s.p,
         paste("x =",format(r2s.value,digits=4, nsmall=4)),
         paste("p =",format(r2s.prob,digits=4, nsmall=4)))

 # Affichage
 r2s.xlab="X"
 r2s.title = paste("Distribution :",r2s.distrib)
 r2s.ylab = expression(f(X==x))
 from = 0

 if(!r2s.is1P) {
   if(!r2s.isDiscrete) {
     from = ifelse(r2s.distrib=="Normal",r2s.param1-4*r2s.param2,0)
     to = ifelse(r2s.distrib=="Normal",r2s.param1+4*r2s.param2,max(r2s.rfunction(1000,r2s.param1,r2s.param2)))
     curve(function(x) r2s.dfunction(x,r2s.param1,r2s.param2),n=1000,from=from,to=to,lwd=2,main=r2s.title,xlab=r2s.xlab,ylab=r2s.ylab)
     if(!r2s.right) {
       r2s.z = seq(from,r2s.value,len=1000)
       for(i in r2s.z) lines(rbind(c(i,0),c(i,r2s.dfunction(i,r2s.param1,r2s.param2))),lwd=2,col="red")
     } else {
         r2s.z = seq(r2s.value,to,len=1000)
         for(i in r2s.z) lines(rbind(c(i,0),c(i,r2s.dfunction(i,r2s.param1,r2s.param2))),lwd=2,col="red") }
     r2s.dum=curve(function(x) r2s.dfunction(x,r2s.param1,r2s.param2),add=TRUE,n=1000,from=from,to=to,lwd=2,main=r2s.title,xlab=r2s.xlab,ylab=r2s.ylab)
   } else {
     from = 0
     to = ifelse(r2s.distrib=="Binomial",r2s.param1,max(r2s.rfunction(1000,r2s.param1,r2s.param2)))
     r2s.z = 0:to
     plot(r2s.z,r2s.dfunction(r2s.z,r2s.param1,r2s.param2),type="h",lwd=2,main=r2s.title,xlab=r2s.xlab,ylab=r2s.ylab)
     if(!r2s.right) {
       for(i in 0:(r2s.value-1)) { lines(rbind(c(i,0),c(i,r2s.dfunction(i,r2s.param1,r2s.param2))),lwd=2,col="red") }
       lines(rbind(c(r2s.value,0),c(r2s.value,r2s.prob-r2s.pfunction(r2s.value-1,r2s.param1,r2s.param2))),lwd=2,col="red")
     } else {
       for(i in r2s.param1:(r2s.value+1)) {
       lines(rbind(c(i,0),c(i,r2s.dfunction(i,r2s.param1,r2s.param2))),lwd=2,col="red") }
       lines(rbind(c(r2s.value,0),c(r2s.value,r2s.prob-1+r2s.pfunction(r2s.value,r2s.param1,r2s.param2))),lwd=2,col="red")}}
 # One parameter distributions
 } else {
   if(!r2s.isDiscrete) {
     from = ifelse(r2s.distrib=="Student",min(r2s.rfunction(1000,r2s.param1)),0)
     to = max(r2s.rfunction(1000,r2s.param1))
     curve(function(x) r2s.dfunction(x,r2s.param1),n=1000,from=from,to=to,lwd=2,main=r2s.title,xlab=r2s.xlab,ylab=r2s.ylab)
     if(!r2s.right) {
       r2s.z = seq(from,r2s.value,len=1000)
       for(i in r2s.z) lines(rbind(c(i,0),c(i,r2s.dfunction(i,r2s.param1))),lwd=2,col="red")
     } else {
         r2s.z = seq(r2s.value,to,len=1000)
         for(i in r2s.z) lines(rbind(c(i,0),c(i,r2s.dfunction(i,r2s.param1))),lwd=2,col="red") }
     r2s.dum=curve(function(x) r2s.dfunction(x,r2s.param1),add=TRUE,n=1000,from=from,to=to,lwd=2,main=r2s.title,xlab=r2s.xlab,ylab=r2s.ylab)
   } else {
     from = 0
     to = max(r2s.rfunction(1000,r2s.param1))
     r2s.z = 0:to
     plot(r2s.z,r2s.dfunction(r2s.z,r2s.param1),type="h",lwd=2,main=r2s.title,xlab=r2s.xlab,ylab=r2s.ylab)
     if(!r2s.right) {
       for(i in 0:(r2s.value-1)) { lines(rbind(c(i,0),c(i,r2s.dfunction(i,r2s.param1))),lwd=2,col="red") }
       lines(rbind(c(r2s.value,0),c(r2s.value,r2s.prob-r2s.pfunction(r2s.value-1,r2s.param1))),lwd=2,col="red")
     } else {
       for(i in to:(r2s.value+1)) {
       lines(rbind(c(i,0),c(i,r2s.dfunction(i,r2s.param1))),lwd=2,col="red") }
       lines(rbind(c(r2s.value,0),c(r2s.value,r2s.prob-1+r2s.pfunction(r2s.value,r2s.param1))),lwd=2,col="red")}}
 }

}

  distribution = gdroplist(names(availDists),horizontal=FALSE,handler=initOptions)
  addhandlerchanged(distribution,handler = function(h,...) {
    theDist = availDists[svalue(h$obj)]

    svalue(param1label) = theParams[[theDist]][1]
    svalue(param1) = theParams[[theDist]][3]

    if(theParams[[theDist]][2] == "") {
      svalue(param2label) = theParams[[theDist]][2]
      svalue(param2) = ""      
      enabled(param2) <- FALSE
    } else {
      enabled(param2) <- TRUE
      svalue(param2label) = theParams[[theDist]][2]
      svalue(param2) = theParams[[theDist]][4]
    }
  })


  calcWhat = gradio(c("Find probability","Find quantile"))
side = gradio(c("to left","to right"))
  param1label = glabel(theParams[[availDists[1]]][1])
  param2label = glabel(theParams[[availDists[1]]][2])
  param1 = gedit(theParams[[availDists[1]]][3],width=5,coerce.with=rpel)
  param2 = gedit(theParams[[availDists[1]]][4],width=5,coerce.with=rpel)

  value  = gedit(width=15, handler = updatePlot,coerce.with=rpel)
  result = glabel("")

BigGroup = ggroup(container= container)
group = ggroup(horizontal = FALSE, container = BigGroup)

tmp = gframe("Distribution", container = group)
distribGroup = glayout(container=tmp)
distribGroup[1,1]=glabel("Law")
distribGroup[1,2]=distribution
distribGroup[2,1]=param1label
distribGroup[2,2]=param1
distribGroup[3,1]=param2label
distribGroup[3,2]=param2
visible(distribGroup)=TRUE

tmp = gframe("Type of calculation", container = group)
add(tmp,calcWhat)

tmp = gframe("Cumulative", container = group)
add(tmp,side)

tmp = gframe("Result", container = group)
  resultGroup = glayout()
  resultGroup[1,1]=glabel("Value")
  resultGroup[1,2]=value
  resultGroup[2,1]=glabel("Result")
  resultGroup[2,2]=result
  visible(resultGroup) <- TRUE
  add(tmp,resultGroup)

addSpring(group)

buttonGroup=ggroup(container=group)
  if(missing(container))
    gbutton("cancel", container=buttonGroup, handler = function(h,...) dispose(container))
addSpring(buttonGroup)
gbutton("update",container=buttonGroup, handler=updatePlot)

add(BigGroup, ggraphics())
  invisible(BigGroup)
}

constructionOfNormal = function(container = gwindow("Construction of normal from X1 + X2 + ... + Xn", visible=FALSE), ...) {
  
  
  availDists = c(Uniform = "unif", Binomial = "binom", Normal = "norm")
  theParams = list(
    "norm" = c("mean","sd",0,1),
    "t" = c("df","ncp","",0),
    "chisq" = c("df","ncp","",0),
    "f" = c("df1","df2","",""),
    "binom"=c("size","prob",1,.5),
    "pois" = c("lambda","",1,""),
    "gamma" = c("shape","rate","",1),
    "beta" = c("shape1","shape2","",""),
    "unif" = c("min","max",0,1)
    )
  
  updatePlot = function(h, ...) {
      
      rfunc = paste("r",availDists[svalue(distribution)],sep="")
      y = do.call(rfunc, list(svalue(sampleSize)*svalue(nvar),
        svalue(param1),svalue(param2)))
      z = rowSums(matrix(y,svalue(sampleSize),svalue(nvar)))
      
      xlab="Values of the variable"
      title = paste("Distribution :",svalue(distribution))
      ylab = "Densities/Probabilities"
      if(svalue(distribution)!="Binomial") {
          hist(z,freq=FALSE,main=title,xlab=xlab,ylab=ylab)
      } else {
          res = plot(table(z)/svalue(sampleSize),main=title,xlab=xlab,ylab=ylab)
      }
      if(svalue(displayFunc)) curve(function(x) dnorm(x,mean(z),sd(z)),
                                    from=min(z),to=max(z),add=TRUE,lwd=2,col=2)
  }
  
  distribution = gdroplist(names(availDists),horizontal=FALSE)
  addhandlerchanged(distribution,handler = function(h,...) {
      theDist = availDists[svalue(h$obj)]

      svalue(param1label) = theParams[[theDist]][1]
      svalue(param1) = theParams[[theDist]][3]

      if(theParams[[theDist]][2] == "") {
          svalue(param2label) = theParams[[theDist]][2]
          svalue(param2) = ""      
          enabled(param2) <- FALSE
      } else {
          enabled(param2) <- TRUE
          svalue(param2label) = theParams[[theDist]][2]
          svalue(param2) = theParams[[theDist]][4]
      }
  })
  


  param1label = glabel(theParams[[availDists[1]]][1])
  param2label = glabel(theParams[[availDists[1]]][2])
  param1 = gedit(theParams[[availDists[1]]][3],width=5,coerce.with=rpel)
  param2 = gedit(theParams[[availDists[1]]][4],width=5,coerce.with=rpel)


  sampleSize  =  gedit("500",width=5,coerce.with=rpel)
  nvar  =  gedit(1, width=5, handler=updatePlot,coerce.with=rpel)

  displayFunc = gcheckbox("Show normal law",handler=updatePlot)


  BigGroup = ggroup(container= container, ...)
  group = ggroup(horizontal = FALSE, container = BigGroup)


  tmp = gframe("Distribution", container = group)
  distribGroup = glayout(container=tmp)
  distribGroup[1,1]=glabel("Law")
  distribGroup[1,2]=distribution
  distribGroup[2,1]= param1label
  distribGroup[2,2]=param1
  distribGroup[3,1]= param2label
  distribGroup[3,2]=param2
  visible(distribGroup)=TRUE

  tmp = gframe("Score", container = group)
  distribSample = glayout(container=tmp)
  distribSample[1,1]=glabel("n for: X1 + ... + Xn")
  distribSample[1,2]=nvar
  distribSample[2,1]=glabel("Number of simulations")
  distribSample[2,2]=sampleSize
  visible(distribSample)=TRUE
  
  tmp = gframe("update", container = group)
  add(tmp,displayFunc)

  addSpring(group)

  buttonGroup=ggroup(container=group)
  if(missing(container))
      gbutton("cancel", container=buttonGroup,
              handler = function(h,...) dispose(container))
  addSpring(buttonGroup)
  gbutton("update",container=buttonGroup, handler=updatePlot)

  add(BigGroup, ggraphics())

  visible(container) <- TRUE
  invisible(BigGroup)
}
