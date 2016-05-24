### Demos:
## histogram -- select bins
## binomial
## normal sampling -- clt
## confidence intervals
## power
## probability calculator
## 
#### IMPLEMENT ME:
## bootstrap
## robustness of regression line

## some in YN.R from Yvonnick Noel yvonnick.noel@uhb.fr


pmgDemo = function(lst, container=NULL) {

  showHelp = function(txt) {
    gmessage(txt)
  }
  
  tmpfile = tempfile()
  
  ## lst is a list with components
  ## widgets
  ## function to draw to a file
  ##

  group = ggroup(container=container)
  lgroup = ggroup(horizontal=FALSE)
  rgroup = ggroup()
  add(group, gpanedgroup(lgroup, rgroup), expand=TRUE)

  ##
  graphic = gimage(); add(rgroup, graphic, expand=TRUE)
  size(graphic) <- c(480,480)
  ## lgroup
  tbl = glayout()
  ctr = 1
  for(i in names(lst$widgets)) {
    tbl[ctr,1] <- glabel(i)
    tbl[ctr,2] <- lst$widgets[[i]]
    ctr <- ctr + 1
  }
  tbl[ctr,2] <- gbutton("help",handler = function(h,...) {
    showHelp(lst$help)
  })
  add(lgroup, tbl)
  visible(tbl) <- TRUE


  ## add the callback  
  theCallback = function(...) {
    file = do.call(lst$f, list(file = tmpfile, lapply(lst$widgets,svalue)))
    svalue(graphic) <- file
  }
  sapply(lst$widgets, function(i) addhandlerchanged(i, handler = theCallback))
  
  ## call callback
  theCallback(file, lapply(lst$widgets,svalue))
  return(group)
}

##################################################



templateForDemo = list(
  widgets = list(

    ),
  f = function(file, valLst) {
    png(file, width=480, height=480)
    on.exit(dev.off())

    return(file)
  },
  help = paste(
    "",
    collapse = "\n")
  )





## open a window showing the possible demos
pmg.teachingDemos = function(...) {

### This is really annoying, but to avoid calling these when defining the list, I put them inside this function
histogramDemo = list(
  widgets = list(
    "data set" = gedit(""),
    "% no. bins" = gslider(from=0.01,to=1, by = 0.01, value=.5)
    ),
  f = function(file, valList) {
    dataSetName = valList[[1]]
    dataSet = svalue(dataSetName)
    if(is.na(dataSetName) || dataSetName == "")
      return(file)

    noBins = as.integer(length(dataSet)*as.numeric(valList[[2]]))
    noBins = max(1, noBins)
    
    png(file, width=480, height=480)
    on.exit(dev.off())

    theTitle = paste(noBins,"bins and ", length(dataSet),"observations" )
    
    print(histogram(dataSet, nint=noBins, xlab=theTitle))
    
    return(file)
  },
  help = paste("Histogram:\n",
    "Drag or type in variable name containing data set",
    " adjust slider to adjust number of bins as a proportion of sample size."
    )
  )
### confidence intervals
confIntervalDemo = list(
  widgets = list(
    mu = gedit("0", coerce.with=as.numeric),
    sigma = gedit("1", coerce.with=as.numeric),
    n = gedit("10", coerce.with=as.numeric),
    "Conf. level" = gradio(c(.80,0.90,0.95), selected=3),
    "No intervals" = gradio(c(10,25,50,100)),
    "Resample" = gbutton("again")
    ),
  f = function(file, valList) {
    mu = valList[[1]]
    sigma = valList[[2]]
    n = valList[[3]]
    confLevel = valList[[4]]
    noIntervals = valList[[5]]

    m = noIntervals
    
    res = matrix(NA,nrow=2,ncol=m)
    for(i in 1:m) res[,i]= t.test(rnorm(n,mu,sigma),conf.level=confLevel)$conf.int
    
    missed = 1 + (res[1,] >= mu | res[2,] <= mu)

    png(file, width = 480, height = 480)
    on.exit(dev.off())

    
    ## make plot
    matplot(res,rbind(1:m,1:m)/m,
          type="l",
          col=c("black","red")[missed],
          lwd=missed,
          lty=missed,
          yaxt="n",
          xlab="",ylab="",
          main=paste(m," ",100*confLevel,"% confidence intervals")  
           )
    abline(v=mu)

    return(file)
  },
  help = paste(
    "Demo of confidence intervals.",
    "Adjust parameters or click 'new sample' to draw a new sample.",
    "The graphic illustrates several simulated confidence intervals. Those",
    "that do not contain the true parameter are flagged.",
    collapse = "\n"
    )

  )

##################################################
CLTDemo = list(
  widgets = list(
    Distribution = gradio(c("Normal","Long tailed","Skewed")),
    n = gslider(from=5,to=50,by=5, value=10),
    m = gslider(from=5,to=25,by=2, value=10),
    Statistic = gradio(c("mean","median")),
    Resample = gbutton("again")
    ),
  f = function(file, valLst) {
    png(file, width=480, height=480)
    on.exit(dev.off())

    Distribution = valLst$Distribution
    n = valLst$n
    m = valLst$m
    Statistic = valLst$Statistic

    
    switch(Distribution,
           "Normal" = {
             curve(function(x) dnorm(x), from=-3, to = 3, lwd=2)
             abline(v=0)
           },
           "Long tailed" = {
             curve(function(x) dt(x,df=3), from = -3.5, to = 3.5, lwd=2)
             abline(v=0)
           },
           "Skewed" = {
             curve(function(x) dexp(x), from = 0, to = 3, lwd=2)
             if(Statistic == "mean")
               abline(v=1)
             else
               abline(v=0.69)
           })

    ySize = switch(Distribution,
      "Normal" = .4,
      "Long tailed" = .4,
      "Skewed" = 1
      )


    ## Now draw the samples
    tmp = numeric(m)
    for(i in 1:m) {
      theSample = switch(Distribution,
        "Normal" = rnorm(n),
        "Long tailed" = rt(n,df=3),
        "Skewed" = rexp(n)
        )
      y = i*ySize/m
      points(theSample, rep(y, n), col=gray(.6))
      tmp[i] = do.call(Statistic, list(theSample))
      points(tmp[i], y, pch=17, cex=2,col=gray(.8)) # lighter
    }

    ## plot sample at bottom with rug
    rug(tmp)
    dens = density(tmp)
    lines(dens$x, dens$y * ySize / max(dens$y), col="red",lwd=2)
    

    return(file)
  },
  help = paste(
    "The Central Limit Theorem.",
    "The graphic illustrates m different saples of size n",
    "illustrated by open circles at the same level.",
    "Each sample is summarized by a statistic marked with a triangle.",
    "These statistics are shown on the x axis using rug() and summarized",
    "with a density plot (not to scale) in red."
    )
  )

## start with an assignment
.binomialDemo.results <- c()
#assign(".binomialDemo.results",c(),envir=.GlobalEnv)
binomialDemo = list(
  widgets = list(
    "Size of sample, n" = gslider(from=5,to=50,by=5,value=10),
    "Success probablity, p" = gslider(from=.05,to=.95,by=.05,value=0.5),
    "Number at a time" = gradio(c(1,10,25)),
    "Sample one" = gbutton("click"),
    "Clear results" = gcheckbox("")
    ),
  f = function(file, valLst) {
    png(file, width=480, height=480)
    on.exit(dev.off())

    n = valLst[[1]]
    p = valLst[[2]]
    m = valLst[[3]]
    clear = valLst[[5]]

    if(!exists(".binomialDemo.results"))
        assign_global(".binomialDemo.results",c())
    if(clear)
      .binomialDemo.results <<- c()

    ## show sample
    
    par(fig=c(0,1,.8,1))
    par(mai=c(0,0,0,0))

    x = rbinom(n, size=1, prob=p) 
    plot.new()
    plot.window(xlim=c(1,n),ylim=c(0,1))
    
    metsCols = c("orange","blue")
    points(1:n,rep(.5,n),pch=16,cex=2,col=metsCols[1 + x])

    ## do I need to add more?
    .binomialDemo.results <<- c(.binomialDemo.results,sum(x))
    if(m > 1)
      .binomialDemo.results <<- c(.binomialDemo.results, rbinom(m-1, n , p))

    ## show results
    par(fig=c(0,1, 0 , .8), new=TRUE)
    par(mar=c(3,2,2,1))
    ## could spice ths up
    if(length(.binomialDemo.results) > 1 ) {
      ## draw a special histogram of the results
      y = .binomialDemo.results
      tbl = table(y)
      plot.new();
      plot.window(xlim=range(y) + c(-1/2,1/2), ylim = c(0,max(tbl)))
      axis(1); axis(2)
      
      cols = c(gray(.9),"blue")[1+as.numeric(names(tbl) == sum(x))]

      names(cols) = names(tbl)
      ## for each val in tbl
      sapply(names(tbl), function(i) {
        x = as.numeric(i)
        polygon(c(x-1/2,x-1/2,x+1/2,x+1/2),c(0,tbl[i],tbl[i],0),col=cols[i])
      })
    }

    return(file)
  },
  help = paste(
    "The binomial distribution. The top figure shows n independent trials ",
    "with success probability p. A success is colored blue.",
    "The number of successes follows a binomial distribution.",
    "The sequence of the number of successess is plotted with a histogram.",
    "To clear out the results, click 'Clear results' twice.",
    collapse="\n")
 )


##################################################
## power

## http://www.amstat.org/publications/jse/v11n3/anderson-cook.html
powerDemo = list(
  widgets = list(
    "Alternative" = gradio(c("less","greater","two.sided")),
    "mu2 - mu1" = gslider(from=-5,to=5, by=.1, value=0),
    "n"  = gslider(from=1, to=25, by=1, value=1),
    alpha = gradio(c(.01, .05, .10), selected=2)
    ),
  f = function(file, valLst) {
    png(file, width=480, height=480)
    on.exit(dev.off())

    HA = valLst[[1]]
    mu2 = valLst[[2]]
    n = valLst[[3]]
    alpha = valLst[[4]]

    x = seq(-5,5, length=1000)
    y1 = dnorm(x, 0, 1/sqrt(n))
    y2 = dnorm(x, mu2, 1/sqrt(n))

    par(mai=c(0,0,0,0))
    plot.new()
    plot.window(xlim=c(-5,5), ylim = c(0, max(y1) + max(y2)))
    axis(1)

    plotPoly = function(x, y, z, toLeft = TRUE, col="blue") {
      i = min(which(x > z))
      n = length(x)
      if(toLeft) {
        polygon(c(x[1],x[1:i],x[i]),c(min(y),y[1:i],min(y)), col=col)
        polygon(c(x[i+1],x[(i+1):n],x[n]),c(min(y),y[(i+1):n],min(y)))
      } else {
        polygon(c(x[1],x[1:i],x[i]),c(min(y),y[1:i],min(y)))
        polygon(c(x[i+1],x[(i+1):n],x[n]),c(min(y),y[(i+1):n],min(y)), col=col)
      }
    }

    Beta = 0
    Beta = switch(HA,
           "less" = {
             z = qnorm(alpha,0,1/sqrt(n))
             abline(v=z)
             plotPoly(x,y1,z, toLeft=TRUE)
             plotPoly(x,y2 + max(y1),z, toLeft=FALSE, col="red")
             ## beta =
             1 - pnorm(z, mu2, 1/sqrt(n))
           },
           "greater" = {
             z = qnorm(1-alpha,0,1/sqrt(n))
             abline(v=z)
             plotPoly(x,y1,z, toLeft=FALSE)
             plotPoly(x,y2 + max(y1),z, toLeft=TRUE, col="red")
             ## return beta
             pnorm(z, mu2, 1/sqrt(n))
           },
           "two.sided" = {
             z = qnorm(1-alpha/2,0,1/sqrt(n))
             abline(v=z)
             abline(v=-z)
             plotPoly(x,y1,-z, toLeft=TRUE)
             plotPoly(x,y1, z, toLeft=FALSE)
             plotPoly(x,y2 + max(y1),z, toLeft=TRUE, col="red")
             plotPoly(x,y2 + max(y1),-z, toLeft=TRUE, col="white")
             ## return beta
              pnorm(z,mu2,1/sqrt(n)) - pnorm(-z,mu2,1/sqrt(n))
           })
    Beta = floor(Beta*100)/100
    text(4,0,label="Null", pos=3)
    text(4,max(y1),label="Alternative", pos=3)
    text(4,max(y1) + 1/2*max(y2), label=paste("Power=",1-Beta))
    return(file)
  },
  help = paste(
    "Power demo. The power of a statistical test indicates the probability of a p-value larger than alpha when the alternative hypothesis is correct.  The power for the z-test depends on the mean postulated for the altervative hypothesis. In this case, the difference mu2 - mu1 is used. The graphic shows in blue the area associated to alpha. and in red the probabliity the NULL is accepted  The power is this latter quantity subtracted from 1.",
    collapse = "\n")
  )

  
  win = pmgWC$new("P M G teaching demos", width=450,height=350)
  group = ggroup(horizontal = FALSE, container=win, expand=TRUE)
  mb = list()
  mb$Demo$"Histogram bin selection"$handler =
    function(...) {
      add(nb,pmgDemo(histogramDemo),label="Histogram")
      svalue(sb) <- ""
    }
  mb$Demo$"histogram and density"$handler = 
    function(...) {
        histogramAndDensity(nb, label="View histogram and density")
#      add(nb, histogramAndDensity(nb), label = "View histogram and density")
      svalue(sb) <- ""
    }
  mb$Demo$"Binomial distribution"$handler =
    function(...) {
      add(nb,pmgDemo(binomialDemo), label="Binomial")
      svalue(sb) <- ""
    }
  mb$Demo$"construction of normal"$handler = 
    function(...) {
        constructionOfNormal(nb, label = "See normal from sum")
        svalue(sb) <-  ""
    }
  mb$Demo$"Central Limit Theorem"$handler =
    function(...) {
      add(nb,pmgDemo(CLTDemo), label="CLT")
      svalue(sb) <-  ""
    }
  mb$Demo$"probability calculator"$handler = 
    function(...) {
      add(nb, probabilityCalculator(NULL), label = "probability calculator")
      svalue(sb) <- ""
    }
  mb$Demo$"Confidence Intervals"$handler =
    function(...) {
      add(nb, pmgDemo(confIntervalDemo), label = "Confidence intervals")
      svalue(sb) <- ""
    }
  mb$Demo$"Power of test"$handler =
    function(...) {
      add(nb, pmgDemo(powerDemo), label = "Power of a test")
      svalue(sb) <-  ""
    }


  add(group, gmenu(mb))
  nb = gnotebook(closebuttons=TRUE)
  add(group, nb, expand=TRUE)

sb <- gstatusbar(gettext("Select a demo from the menubar above"), container=group)
}
                
