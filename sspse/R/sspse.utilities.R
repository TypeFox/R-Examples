#' @keywords internal
mode.density <- function(x,lbound=min(x,na.rm=TRUE),ubound=max(x,na.rm=TRUE), smooth=0.35, h=0){
#     posdensN <- locfit(~ x, alpha=c(2*smooth,0.3))
#     posdensN <- locfit(~ x)
#     xp <- floor(lbound):ceiling(ubound)
      x <- x[!is.na(x)]
      if(length(x)==0){return(NA)}
      if(length(x[x>=lbound&x<=ubound])==0){return(NA)}
      if(length(x[x>=lbound&x<=ubound])==1){return(x[x>=lbound&x<=ubound])}
#     posdensN <- try(locfit( ~ lp(x, nn=2*smooth, h=h, maxk=500)),silent=TRUE)
#     if(inherits(posdensN,"try-error")){
       posdensN <- density(x, from=lbound, to=ubound)
#     }else{
#      xp <- seq(lbound, ubound, length=10000)
#      posdensN <- try(predict(posdensN, newdata=xp),silent=TRUE)
#      if(inherits(posdensN,"try-error")){
#       posdensN <- density(x, from=lbound, to=ubound)
#      }else{
#       posdensN <- list(x=xp,y=posdensN)
#      }
#     }
      posdensN$x[which.max(posdensN$y)]
}
#' @keywords internal
ll.density <- function(x,lbound=min(x,na.rm=TRUE),ubound=max(x,na.rm=TRUE), smooth=0.35, h=0){
#     posdensN <- locfit(~ x, alpha=c(2*smooth,0.3))
#     posdensN <- locfit(~ x)
#     xp <- floor(lbound):ceiling(ubound)
      x <- x[!is.na(x)]
      if(length(x)==0){return(NA)}
      xp <- seq(lbound, ubound, length=10000)
#     posdensN <- try(locfit( ~ lp(x, nn=2*smooth, h=h, maxk=500)),silent=TRUE)
      posdensN <- .catchToList(bgk_kde(x,n=2^(ceiling(log(ubound-lbound)/log(2))),MIN=lbound,MAX=ubound))
      if(!is.null(posdensN$error)){
       posdensN <- density(x, from=lbound, to=ubound)
       list(x=xp,y=posdensN)
      }else{
#      posdensN <- try(predict(posdensN, newdata=xp),silent=TRUE)
       posdensN <- .catchToList(spline(x=posdensN$value[1,],y=posdensN$value[2,],xout=xp)$y)
       if(!is.null(posdensN$error)){
        posdensN <- density(x, from=lbound, to=ubound)
        list(x=xp,y=posdensN)
       }else{
        list(x=xp,y=posdensN$value)
       }
      }
}
#' @keywords internal
HPD.density <- function(x,lbound=min(x,na.rm=TRUE),ubound=max(x,na.rm=TRUE)){
      x <- x[!is.na(x)]
      if(length(x)==0){return(c(NA,NA))}
      posdensN <- density(x, from=lbound, to=ubound)
      cy <- cumsum(posdensN$y/sum(posdensN$y))
      cy <- c(posdensN$x[which.max(cy>0.025)],
              posdensN$x[which.max(cy>0.975)])
      if(is.na(cy[1])) cy[1] <- lbound
      if(is.na(cy[2])) cy[2] <- ubound
      cy
}

.catchToList <- function(expr) {
  val <- NULL
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, w$message)
    invokeRestart("muffleWarning")
  }
  myError <- NULL
  eHandler <- function(e) {
    myError <<- e$message
    NULL
  }
  val <- tryCatch(withCallingHandlers(expr, warning = wHandler), error = eHandler)
  list(value = val, warnings = myWarnings, error=myError)
} 
