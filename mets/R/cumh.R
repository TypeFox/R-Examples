
cumh <- function(formula,data,...,time,
                 timestrata=quantile(data[,time],c(0.25,0.5,0.75,1)),
                 cens.formula=NULL,cens.model="aalen",
                 cumulative=FALSE, 
                 silent=FALSE) {

  time.  <- substitute(time)
  if (!is.character(time.)) time. <- deparse(time.)
  time <- time.
  if (!is.null(cens.formula)) {    
    m <- match.call(expand.dots = TRUE)[1:3]
    Terms <- terms(cens.formula, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    censMod <- eval(m, parent.frame())
    censtime <- model.extract(m, "response")
    status <- censtime[,2]
###  data[,"_status"] <- 
  }

  res <- list(); i <- 0
  ht <- c()
  outcome <- as.character(terms(formula)[[2]])
  coefs <- c()
  y0 <- data[,outcome]
  for (i in seq(length(timestrata))) {
    t <- timestrata[i]
    data[,outcome] <- y0
    newdata <- data
    if (!cumulative) {
      if (i==1) {
        idx <- data[,time]<t
      } else {
        idx <- (timestrata[i-1]<=data[,time] & data[,time]<t)
      }
    } else {
      newdata <- data
      newdata[,]
      newdata[,outcome] <- data[,outcome]*(data[,time]<t)
      if (!is.null(cens.formula)) {
        
      }      
    }
    if (!silent) {
      message(t)
    }
    if (!cumulative)
      res[[i]] <- summary(twinlm(formula,data=data[idx,],...))
    else {
      
      res[[i]] <- summary(bptwin(formula,data=data,...))
    }

    coefs <- c(coefs, list(res[[i]]$all))
    ht <- rbind(ht,c(t,res[[i]]$heritability[1,]))
  }
  
  coeftype <- c()
  for (i in seq(nrow(coefs[[1]]))) {
    rr <- matrix(unlist(lapply(coefs,function(z) z[i,])),ncol=3,byrow=TRUE)
    colnames(rr) <- colnames(coefs[[1]])
    rr <- cbind(time=timestrata,rr)
    coeftype <- c(coeftype,list(rr)); names(coeftype)[length(coeftype)] <- rownames(coefs[[1]])[i]
  }
  
  rownames(ht) <- timestrata
  colnames(ht) <- c("time","Heritability","Std.Err","2.5%","97.5%")
  res <- (list(time=timestrata, ht=ht,models=res,coef=coefs,
               coeftype=coeftype, timevar=time))
  class(res) <- "cumh"
  res
}


##' @export
summary.cumh <- function(object,...) object 

##' @export
print.cumh <- function(x,type=seq(nrow(x$coef[[1]])),...) {
  for (i in type) {    
    cat(i, ": ", names(x$coeftype)[i], "\n",sep="")
    rr <- matrix(unlist(lapply(x$coef,function(z) z[i,])),ncol=3,byrow=TRUE)
    colnames(rr) <- colnames(x$coef[[1]])
    rr <- cbind(time=x$time,rr)
    print(rr)
  }
  invisible(x)
}


##' @export
plot.cumh <- function(x,...,type=1,lwd=2,col,fillcol,alpha=0.2,ylim=c(0,1),xlab=x$timevar,ylab="Heritability",idx=seq(nrow(x$ht)),legend=TRUE,legendpos="topleft") {

  add <- FALSE
  if (missing(col)) col <- seq(length(type))
  if (alpha>0 & missing(fillcol)) fillcol <- Col(col,alpha)
  count <- 0
  for (tt in type) {
    count <- count+1
    zz <- x$coeftype[[tt]][idx,,drop=FALSE]
    if (!add) {    
      plot(zz[,1:2,drop=FALSE],type="l",ylim=ylim,lwd=lwd,
           ylab=ylab,xlab=xlab,col=col[count],...)
    }
    add <- TRUE
    xx <- with(x, c(zz[,1],rev(zz[,1])))
    yy <- with(x, c(zz[,3],rev(zz[,4]))) 
    polygon(xx,yy,col=fillcol[count])
    lines(zz[,1:2,drop=FALSE],lwd=lwd,col=col[count],...)
  }
  if (legend) graphics::legend(legendpos,names(x$coeftype)[type],col=col,lwd=lwd,lty=1)
  invisible(x)    
}


