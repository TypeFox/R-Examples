# make.args code courtesy of Matthew Lundberg
make.args <- function(..., PRE.ARGS=list(), POST.ARGS=list()) {
  a <- list()
  l <- c(PRE.ARGS, list(...), POST.ARGS)
  for (name in unique(names(l))) {
    a[[name]] <- l[[name]] # First occurrence will be found.
  }
  return(a)
}

`plot.randomLCA` <-
  function(x,...,graphtype=ifelse(x$random,"marginal","conditional"),conditionalp=0.5,classhorizontal=TRUE) {
    if (!inherits(x, "randomLCA"))
      stop("Use only with 'randomLCA' objects.\n")
    if (missing(graphtype)) graphtype <- ifelse(x$random,"marginal","conditional")
    # calculate everything to be plotted
    graphdata <- NULL
    if (graphtype=="marginal") {
      graphdata <- calcMargProb(x)
      graphdata <- cbind(perc=rep(0,dim(graphdata)[1]),graphdata)
    }
    if (graphtype=="conditional")  graphdata <- calcCondProb(x,conditionalp)
    if (graphtype=="conditional2") graphdata <- calcCond2Prob(x,conditionalp)
    # set up the x axis labels
    #    browser()
    if (x$level2) blocksize <- x$level2size
    else blocksize <- x$blocksize
    noblocks <- dim(x$patterns)[2] %/% blocksize
    if (x$level2 | (blocksize!=dim(x$patterns)[2])) {
      thenames <- names(x$patterns)
      thenames <- strsplit(thenames,"\\.")
      temp <- NULL
      for (i in 1:noblocks) {
        temp <- c(temp,thenames[[1+(i-1)*blocksize]][2])
      }
      thenames <- temp
    } else thenames <- names(x$patterns)
    # first decide if there are multiple blocks, otherwise plot all classes on one graph
    if (blocksize==dim(x$patterns)[2]) {
      if (graphtype=="marginal") {
        arg.list <- make.args(...,
                              x = as.formula("outcomep~outcome"),
                              group=graphdata$class,
                              data=graphdata,
                              ylim=c(-0.05,1.05),
                              xlab="Outcome",
                              ylab="Outcome Prob.",
                              scales=list(x=list(at=1:length(thenames),labels=thenames))
        )
        print(do.call(xyplot, arg.list))
      }
      if ((graphtype=="conditional") || (graphtype=="conditional2")) {
        if (length(conditionalp)>1) {
          arg.list <- make.args(...,
                                x = as.formula("outcomep~outcome|perc"),
                                group=graphdata$class,
                                data=graphdata,
                                ylim=c(-0.05,1.05),
                                xlab="Outcome",
                                ylab="Outcome Prob.",
                                scales=list(x=list(at=1:length(thenames),labels=thenames),alternating=FALSE)
          )
          print(do.call(xyplot, arg.list))
        }
        else {
          arg.list <- make.args(...,
                                x = as.formula("outcomep~outcome"),
                                group=graphdata$class,
                                data=graphdata,
                                ylim=c(-0.05,1.05),
                                xlab="Outcome",
                                ylab="Outcome Prob.",
                                scales=list(x=list(at=1:length(thenames),labels=thenames))
          )
          print(do.call(xyplot, arg.list))
        }
      }
    }
    else {
      if (classhorizontal) {
        if (graphtype=="marginal"){
          arg.list <- make.args(...,
                                x = as.formula("outcomep~block|class"),
                                group=graphdata$outcome,
                                data=graphdata,
                                ylim=c(-0.05,1.05),
                                xlab="Block",
                                ylab="Outcome Prob.",
                                scales=list(x=list(at=1:length(thenames),labels=thenames),alternating=FALSE)
          )
          print(do.call(xyplot, arg.list))
        }
        
        if ((graphtype=="conditional") || (graphtype=="conditional2")) {
          if (length(conditionalp)>1) {
            arg.list <- make.args(...,
                                  x = as.formula("outcomep~block|class*perc"),
                                  group=graphdata$outcome,
                                  data=graphdata,
                                  ylim=c(-0.05,1.05),
                                  xlab="Block",
                                  ylab="Outcome Prob.",
                                  scales=list(x=list(at=1:length(thenames),labels=thenames),alternating=FALSE)
            )
            print(do.call(xyplot, arg.list))
          }		else {
            arg.list <- make.args(...,
                                  x = as.formula("outcomep~block|class"),
                                  group=graphdata$outcome,
                                  data=graphdata,
                                  ylim=c(-0.05,1.05),
                                  xlab="Block",
                                  ylab="Outcome Prob.",
                                  scales=list(x=list(at=1:length(thenames),labels=thenames),alternating=FALSE)
            )
            print(do.call(xyplot, arg.list))
          }
        }
      }
      else  {
        if (graphtype=="marginal") {
          arg.list <- make.args(...,
                                x = as.formula("outcomep~block|class"),
                                group=graphdata$outcome,
                                data=graphdata,
                                ylim=c(-0.05,1.05),
                                xlab="Block",
                                ylab="Outcome Prob.",
                                scales=list(x=list(at=1:length(thenames),labels=thenames),alternating=FALSE)
          )
          print(do.call(xyplot, arg.list))
        }
        if ((graphtype=="conditional") || (graphtype=="conditional2")) {
          if (length(conditionalp)>1) 
            arg.list <- make.args(...,
                                  x = as.formula("outcomep~block|perc*class"),
                                  group=graphdata$outcome,
                                  data=graphdata,
                                  ylim=c(-0.05,1.05),
                                  xlab="Block",
                                  ylab="Outcome Prob.",
                                  scales=list(x=list(at=1:length(thenames),labels=thenames),alternating=FALSE)
            )
          print(do.call(xyplot, arg.list))
        }
        else {
          if (length(conditionalp)>1) 
            arg.list <- make.args(...,
                                  x = as.formula("outcomep~block|perc*class"),
                                  group=graphdata$outcome,
                                  data=graphdata,
                                  ylim=c(-0.05,1.05),
                                  xlab="Block",
                                  ylab="Outcome Prob.",
                                  scales=list(x=list(at=1:length(thenames),labels=thenames),alternating=FALSE)
            )
          print(do.call(xyplot, arg.list))
        }
      }
    }
  }

