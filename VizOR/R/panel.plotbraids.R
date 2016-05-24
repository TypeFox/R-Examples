##' Panel function for function 'plotbraids'
##'
##' This panel function manages details of braid geometry.  Key data
##' preparation steps are performed by the high-level 'plotbraids'
##' function.
##' @title panel.plotbraids
##' @param x 
##' @param y 
##' @param subscripts 
##' @param ... 
##' @param formula 
##' @param data 
##' @param idvar 
##' @param stratify 
##' @param steps 
##' @return None
##' @author David C. Norris
##' @seealso \code{\link{plotbraids}}
##' @keywords internal hplot
##' @export panel.plotbraids
panel.plotbraids <- function(x, y, subscripts, ..., formula, data, idvar, stratify, steps, outside){
  statevar <- as.character(formula[[2]]) # Let our formula be of the form
  timevar <- as.character(formula[[3]])  # 'statevar ~ timevar [| condvar]',
  condvar <- tryCatch(as.character(formula[[3]][[3]]), # with optional conditioning variable.
                      error=function(e) NULL)
  if(length(timevar)>1)
    timevar <- timevar[[2]]
  data.long <- data[subscripts,c(idvar, statevar, timevar)]
  ## Reshape the longitudinal 'data' to wide-form
  data <- reshape(data.long, v.names=statevar, timevar=timevar, idvar=idvar, direction="wide")
  if(any(is.na(data))){
    ## If data[subscripts,] turns out to have been unbalanced,
    ## then append a 'CENSORED' level to the statevar variable,
    ## redo the reshape operation, and replace the resulting NA
    ## values with 'CENSORED':
    data.long[[statevar]] <- colored(data.long[[statevar]],
                                     c(CENSORED="white", key(data.long[[statevar]])))
    data <- reshape(data.long, v.names=statevar, timevar=timevar, idvar=idvar, direction="wide")
    statecolumns <- grep(paste(statevar,".",sep=""), names(data))
    NAs <- is.na(data[,statecolumns])
    data[,statecolumns][NAs] <- 'CENSORED'
  }
  tx.key <- key(data.long[[statevar]])
  tx.key <- tx.key[!is.na(names(tx.key))] # Assumption: user wants only explicit factor levels plotted (not NAs)
  ## Start data frame for lane markers
  lanes <- data.frame(trt=names(tx.key))
  ## Aggregate the data to get a frequency column in place of patnum
  N <- dim(data)[1]
  data <- aggregate(list(x=rep(1/N,N)), # TODO: Might it not be clearer to call this 'y', not 'x'?
                    by=data[,paste(statevar,steps,sep=".")],
                    FUN=sum)
  last <- dim(data)[1]
  for(i in steps){
    col <- paste(statevar, i, sep=".")
    data <- data[order(data[[col]]),]
    top <- paste(col, "top", sep=".")
    bot <- paste(col, "bot", sep=".")
    data[[top]] <- cumsum(c(0, data$x[-last]))
    data[[bot]] <- cumsum(data$x)
    ## Add the ith period to 'lanes':
    x.i <- aggregate(data$x, by=data.frame(trt=data[[col]]), FUN=sum)
    names(x.i)[names(x.i)=='x'] <- paste('x', i, sep=".")
    lanes <- merge(lanes, x.i, by='trt', all=TRUE)
    nas <- is.na(lanes[,length(lanes)])
    lanes[nas,length(lanes)] <- 0
  }
  names(lanes)[1] <- statevar # This name change enables merges with 'data'
  ## We calculate the lane markers next, since the 'strata'
  ## concept are defined by straightening these out.
  lanes$x.max <- apply(lanes[,-1], 1, max)
  ## Set a minimum (relative) height for the lanes,
  ## to ensure that the stratum labels fit.
  ##lanes$x.max <- pmax(lanes$x.max, 0.05)
  lanes <- lanes[match(names(tx.key), lanes[,1]),]
  lanes[,-1] <- cumsum(lanes[,-1])
  ## The x.max values now give the tops of each of the strata,
  ## and may (optionally) be used at this point to adjust the
  ## x.n.(top|bot) columns.
  if(stratify){
    ## Convert lanes$x.n to whitespace adjustments
    lanes[,-c(1,length(lanes))] <- lanes$x.max - lanes[,-c(1,length(lanes))]
    ## TODO: Consider doing this calculation with reshaped (long-form) frames
    for(i in steps){
      trt.i <- paste(statevar,i,sep=".")
      x.i <- paste("x",i,sep=".")
      data <- merge(data, lanes[,c(statevar, x.i)], by.x=trt.i, by.y=statevar)
    }
    trts.bot <- paste(statevar, steps, "bot", sep=".")
    trts.top <- paste(statevar, steps, "top", sep=".")
    adjs <- paste("x", steps, sep=".")
    data[,trts.bot] <- data[,trts.bot] + data[,adjs]
    data[,trts.top] <- data[,trts.top] + data[,adjs]
    ## TODO: Consider transforming 'lanes' in parallel as check on correctness
    lanes[,-c(1,length(lanes))] <- lanes$x.max
    lanes[,-1] <- lanes[,-1]/max(lanes$x.max)
  }
  data <- data[order(data[[paste(statevar, steps[1], sep=".")]]),]
  numeric.columns <- -seq(length(steps))
  data[,numeric.columns] <- data[,numeric.columns]/max(data[,numeric.columns])
  ## Plot!
  snake <- function(y, r=0.2, S=51){ # Smooth sinusoidal transitions between points
    curve.t <- seq(-r, r, length=S)
    curve.y <- (1 + sin(seq(-pi/2, pi/2, length=S)))/2
    ans <- list(t=0, y=y[1])
    for(i in 2:length(y)){
      ans$t <- c(ans$t, curve.t + (i-1))
      ans$y <- c(ans$y, rep(y[i-1], S) + (y[i] - y[i-1])*curve.y)
    }
    ans$t <- c(ans$t, i)
    ans$y <- c(ans$y, y[length(y)])
    ans
  }
  ## Draw the braids
  for(k in seq(nrow(data))){
    y.top <- unlist(data[k,grep(paste(statevar,"[1-9]","top",sep="."),names(data))])
    y.bot <- unlist(data[k,grep(paste(statevar,"[1-9]","bot",sep="."),names(data))])
    S.top <- snake(y.top)
    S.bot <- snake(y.bot)
    state.1 <- paste(statevar,steps[1],sep=".")
    color.key <- key(data[[state.1]])
    grid.polygon(x=c(S.top$t, rev(S.bot$t)),
                 y=c(S.top$y, rev(S.bot$y)),
                 default.units="native",
                 gp=gpar(col=rgb(0,0,0,alpha=0),
                   fill=color.key[[data[k,state.1]]],
                   lwd=1))
  }
  ## Draw the lanes
  for(k in seq(nrow(lanes))){
    y <- unlist(lanes[k,grep("x[.][1-9]",names(lanes))])
    S <- snake(y)
    grid.lines(x=S$t, y=S$y,
               default.units="native",
               gp=gpar(col='black', lwd=1))
  }
  ## Position and draw the y axis labels
  ylabels <- as.character(lanes[[statevar]])
  yat <- lanes[[paste("x",steps[1],sep=".")]]
  yat <- (c(0,yat[-length(yat)]) + yat)/2
  ## Next 2 lines show how to get height of label text, should it prove necessary to implement fine adjustments:
  ##label.height <- grobHeight(textGrob("X", gp=do.call("gpar", trellis.par.get("axis.text"))))
  ##label.height <- convertUnit(label.height, 'npc', valueOnly=TRUE) # TODO: Understand why 'native' is wrong
  pushViewport(viewport(clip="off")) # prevent clipping of labels for thin strata at top or bottom panel edges
  panel.axis(side="left",
             labels=c(ylabels,ylabels), # Strangely, doubling these vectors
             at=c(yat,yat),             # seems to be necessary.
             draw.labels=TRUE, outside=outside, ticks=FALSE)
  popViewport()
}
