"ci.plot" <-
function(lm.object, ...)
  UseMethod("ci.plot")

"ci.plot.lm" <-
function(lm.object,
         xlim=range(data[, x.name]),
         newdata,
         conf.level=.95,
         data=model.frame(lm.object),
         newfit,
         ylim,
         pch=19,
         lty=c(1,3,4,2),
         lwd=2,
         main.cex=1,
         main=list(paste(100*conf.level,
           "% confidence and prediction intervals for ",
           substitute(lm.object), sep=""), cex=main.cex), ...
         ) {
  formula.lm <- formula(lm.object)
  x.name <- as.character(formula.lm[[3]])
  y.name <- as.character(formula.lm[[2]])
  missing.xlim <- missing(xlim)       ## R needs this
  missing.ylim <- missing(ylim)       ## R needs this
  missing.newdata <- missing(newdata) ## R needs this
  default.newdata <- data.frame(seq(xlim[1], xlim[2], length=51))
  names(default.newdata) <- x.name
  if (missing.xlim) xlim <- xlim + diff(xlim)*c(-.02,.02) ## needed
  if (missing.newdata) {
    newdata <- default.newdata
    newdata.x <- numeric()
  }
  else {
    if (is.na(match(x.name, names(newdata))))
      stop(paste("'newdata' must be a data.frame containing a column named '",
                 x.name, "'", sep=""))
    if (missing.xlim)
      xlim=range(xlim, newdata[[x.name]])
    newdata.x <- as.data.frame(newdata)[,x.name]
    newdata <- rbind(as.data.frame(newdata)[,x.name, drop=FALSE],
                     default.newdata)
    newdata <- newdata[order(newdata[,x.name]), , drop=FALSE]
  }
  if (missing.xlim) xlim <- xlim + diff(xlim)*c(-.02,.02) ## repeat is needed
  if (missing(newfit)) {
           new.p <-
             predict(lm.object, newdata=newdata,
                     se.fit=TRUE, level=conf.level,
                     interval = "prediction")
           new.c <-
             predict(lm.object, newdata=newdata,
                     se.fit=TRUE, level=conf.level,
                     interval = "confidence")
           newfit <- new.p
           newfit$ci.fit <- new.c$fit[,c("lwr","upr"), drop=FALSE]
           dimnames(newfit$ci.fit)[[2]] <- c("lower","upper")
           attr(newfit$ci.fit,"conf.level") <- conf.level
           newfit$pi.fit <- new.p$fit[,c("lwr","upr"), drop=FALSE]
           dimnames(newfit$pi.fit)[[2]] <- c("lower","upper")
           attr(newfit$pi.fit,"conf.level") <- conf.level
           newfit$fit <- newfit$fit[,"fit", drop=FALSE]
         }
  tpgsl <- trellis.par.get("superpose.line")
  tpgsl$lty <- lty
  tpgsl$lwd <- lwd
  tpgsl$col[1] <- 0
  tpgsl <- Rows(tpgsl, 1:4)
  if (missing.ylim) {
    ylim <- range(newfit$pi.fit, data[,y.name])
    ylim <- ylim + diff(ylim)*c(-.02,.02) ## needed
  }
  xyplot(formula.lm, data=data, newdata=newdata, newfit=newfit,
         newdata.x=newdata.x,
         xlim=xlim, ylim=ylim, pch=pch,
         par.settings=list(superpose.line=tpgsl),
         panel=function(..., newdata.x) {
           panel.ci.plot(...)
           if (length(newdata.x) > 0)
             panel.rug(x=newdata.x)
         },
         main=main,
         key=list(border=TRUE,
           space="right",
           text=list(c("observed","fit","conf int","pred int")),
           points=list(
             pch=c(pch,32,32,32),
             col=c(trellis.par.get("plot.symbol")$col, tpgsl$col[2:4])
             ),
           lines=tpgsl),
         ...)
}

## source("~/HH-R.package/HH/R/ci.plot.R")
