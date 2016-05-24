## This is a practice for making better widgets for beginning students
## It seems like a good student widget should have:
## * a graphic
## * dynamic response to change of variables
## * visible output, students don't like quiet GUIs

## in misc.R
## is.gdataframecolumn = function(obj) {
##   if(class(obj)[1] == "GtkTreeViewColumn")
##     return(TRUE)
##   else
##     return(FALSE)
## }

## make a graphic save to png file, return file name
dp = function(x,probs=NULL,
  cols = gray(seq(.25, .75, length=length(probs))),
  f = tempfile()
  ) {
  width = 200; height= 200
  png(f, pointsize=4, width=width, height=height)
  par(mai=c(0,0,0,0))
  
  if(missing(x)) {
    ## make an empty graphic
    plot.new(); plot.window(xlim=c(0,1),ylim=c(0,1))
    text(.5,.5,"Add variable")
  } else {
    ## make a pretty graphic

    d = density(x, na.rm=TRUE)
    

    plot(d, xlab="",ylab="",main="")
    rug(x)

    cumSums = cumsum(d$y[-1]*diff(d$x))
    cumSums[1] = 0; cumSums[512] = 1
    if(!is.null(probs)) {
      n = length(probs)
      probs = sort(probs)
      for(i in n:1) {
        j = max(which(cumSums <= probs[i]))
        polygon(c(d$x[1:j],rev(d$x[1:j])),c(d$y[1:j],rep(0,j)),col=cols[i])
      }
    } 
  }

  ## return filename of graphic
  dev.off()
  return(f)
}


## Write a widget to display quantiles
## then abstract to make general widget

quantileWidget = function(container=NULL) {

  ## Define layout and primary widgets:
  ## f
  ## plotGraphic
  ## xVar
  ## probsVar
  ## outputArea
  group = ggroup(horizontal=FALSE, container=container)
  obj=group                                # will return this, can use to add things to with
                                        # tag()
  f = dp()                              # default value
  tag(obj,"graphicFilename") <- f
  plotGraphic = gimage(f, container=group)
  gseparator(horizontal=TRUE, container=group)

  tbl = glayout()
  tbl[1,1] = glabel("x:")
  xVar = glabel("Add variable here", editable=TRUE)
  font(xVar)<-list(style="bold")

  tag(obj,"xVarData") <- NULL
  tbl[1,2] = xVar

  tbl[2,1] <- glabel("probs:")
  probsVar = gdroplist(c("c()","c(.25,.5,.75)","seq(.2,.8,by=.2)",
    "seq(.1,.9,by=.1)","seq(.05,.95,by=.05)"))
  tag(obj,"probsVarData") <- NULL
  tbl[2,2] = probsVar

  add(group,tbl)
  visible(tbl)<-TRUE

  outputArea = gtext("")
  add(group, outputArea, expand=TRUE)

  ## Now add actions
  addhandlerchanged(xVar,
                    handler = function(h,...) {
                      ids = tag(obj,"dropHandlers")
                      if(length(ids) > 0) {
                        removehandler(obj,ids)
                        tag(obj,"dropHandlers") <- list()
                      }
                      tag(obj, "xVarData") <- svalue(h$obj)
                      ## put popup on 1
                      ## svalue(tag(obj,"actionPopup"),index=TRUE) <- 1
                      
                      update()
                    })
  adddroptarget(xVar,
                handler=function(h, ...) {
                  tag(obj,"xVarData") <- h$dropdata
                  svalue(xVar) <- id(h$dropdata)
                  ## put popup on 1
                  ## svalue("actionPopup",index=TRUE) <- 1
                  
                  update()
                  
                  ## now bind to be dynamic *if* a treeviewcolumn
                  if(is.gdataframecolumn(h$dropdata)) {
                    view.col = h$dropdata
                    id = addhandlerchanged(view.col,
                      signal = "edited",
                      handler=function(h,...) update()
                      )
                    dropHandlers = tag(obj,"dropHandlers")
                    dropHandlers[[length(dropHandlers)+1]] = list(
                                  view.col = view.col,
                                  id = id
                                  )
                    tag(obj,"dropHandlers") <- dropHandlers
                  }
                })

  ## for probs we have to parse text to make sure it is valid R code
  addhandlerchanged(probsVar, function(h,...) {
    errMsg = "The value of Probs: should evaluate to a vector of probabilities, e.g. c(.25,.5,.75)"

    val = svalue(probsVar)              # command like "c(1,2,3)"
    vals = try(eval(parse(text=val), envir=.GlobalEnv), silent=TRUE)
    if(inherits(vals, "try-error")) {
      cat(errMsg); cat("\n")
      return()
    } else {
      if(any(vals >1) || any(vals < 0)) 
        cat("Probs: are probabilities between 0 and 1\n")
      vals = vals[vals <=1 && vals >= 0] # these are probabilities
      vals = sort(vals)                 # sort if funny
      if(length(vals) >= 1) {
        tag(obj,"probsVarData") <- vals
        update()
      }
    }
  })

  update = function() {
    ## update the graph, update the function
    x = svalue(tag(obj,"xVarData"))
    if(!is.numeric(x)) {
      cat("the x variable is not numeric\n")
      return()
    }
    
    probs = tag(obj,"probsVarData")
    f = tag(obj,"graphicFilename")
    f = dp(x, probs, f=f)
    svalue(plotGraphic) <- f            # update graphic

    ## get numeric summaries, trim down digits
    oldDigits = getOption("digits")
    options("digits"=2)
    if(length(probs) >= 1)
      outValue = capture.output(quantile(x, probs=probs, na.rm=TRUE))
    else
      outValue = "specify a value for 'probs' (a vector of probabilities)"
    options("digits"=oldDigits)

    dispose(outputArea)
    add(outputArea, "Quantiles:\n",font.attr=c(style="monospace",color="blue"))
    add(outputArea, outValue, font.attr=c(style="monospace"))

    svalue(outputArea) <- c("Quantiles:\n",outValue)
  }

  ## return the widget
  return(obj)
}
                                
