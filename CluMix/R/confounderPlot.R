confounderPlot <-
function(data, S, x, y, labels, returnS=FALSE, ...){
# data: data frame with variables of interest
# S: similarity matrix; if missing it will be calculated from data
# x: predictor of main interest, for which confounders / collinearities shall be detected
# y: outcome variable
# labels: variable names used for plotting - have to be in corresponding order with columns of data!
# returnS: shall similarity matrix be returned
  
  # similarity matrix
  if(missing(S))
    S <- similarity.variables(data)
  
  if(!identical(names(data), colnames(S)))
    stop("names of 'data' and 'S' must coincide")
  
  if(missing(labels))
    labels <- names(data)
  names(labels) <- names(data)
  
  # color: categorical/continuous variables
  dc <- sapply(data, data.class)
  col <- ifelse(dc == "numeric", 1, 4)
  
  # highlight x and y (bold)
  font <- ifelse(names(data) %in% c(x,y), 2, 1)
  
  # plot
  par(mar=c(5,4,4,7))
  #plot(0:1, 0:1, type="n", xlab=paste("similarity to", labels[x]), ylab=paste("similarity to", labels[y]), ...)
  #text(x=S[x,], y=S[y,], labels=labels, col=col, font=font, xpd=T)
  plot(S[x,], S[y,], ylim=c(0,1.05), xlim=c(0,1.05), frame=FALSE, pch=16, col=col, xlab=paste("similarity to", labels[x]), ylab=paste("similarity to", labels[y]), ...)
  text(x=S[x,], y=S[y,] + 0.03, labels=labels, col=col, font=font, xpd=T)
  rect(0, 0, 1, 1, lty=3)
  legend(1.05, 0.9, c("continuous", "categorical"), col=c(1,4), text.col=c(1,4), xpd=T)
  
  if(returnS)
    return(S)
}
