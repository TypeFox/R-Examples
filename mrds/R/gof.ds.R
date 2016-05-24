#' Compute chi-square goodness-of-fit test for ds models
#'
#' @param model \code{ddf} model object
#' @param breaks distance cut points
#' @param nc number of distance classes
#'
#' @return list with chi-square value, df and p-value
#' @seealso ddf.gof
#' @author Jeff Laake
#' @importFrom stats pchisq
gof.ds <- function(model,breaks=NULL,nc=NULL){
  # Functions used: predict.ds

  width <- model$meta.data$width
  left <- model$meta.data$left
  xmat <- model$ds$aux$ddfobj$xmat
  n <- dim(xmat)[1]

  # If breaks are given use those; otherwise create them
  if(is.null(breaks)){
    breaks <- model$ds$aux$breaks
  }
  if(is.null(nc)){
    nc <- round( sqrt(length(xmat$distance)))
  }
  if(is.null(breaks)){
    breaks <- left + ((width-left)/nc)*(0:nc)
  }else{
    nc <- length(breaks)-1
  }

  # Get predicted values for ds component
  expected.1 <- rep(0,nc)
  for(j in 1:nc){
    expected.1[j] <- sum(predict(model,compute=TRUE,integrate=TRUE,
                                 int.range=matrix(c(breaks[j],breaks[j+1]),
                                 nrow=1))$fitted/model$fitted,
                         na.rm=TRUE)
  }

  # Compute observed values of distance bins
  observed.count.1 <- table(cut(xmat$distance,breaks,include.lowest=TRUE))
  chisq.1 <- sum((observed.count.1-expected.1)^2/expected.1,na.rm=TRUE)
  df.1 <- nc-1-length(model$par)
  if(df.1 <= 0){
    df.1 <- NA
    p.1 <- NA
  }else{
    p.1 <- 1-pchisq(chisq.1,df.1)
  }

  # build and return list
  return(list(chi1=list(observed=observed.count.1,
              expected=expected.1,
              chisq=chisq.1,
              p=p.1,
              df=df.1)))
}
