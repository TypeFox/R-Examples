##' @export
plotcr <- function(x,col,lty,legend=TRUE,which=1:2,cause=1:2,
                   ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...) {
    if (inherits(try(find.package("prodlim"),silent=TRUE),"try-error")) {
        stop("Needs prodlim")
    }

  dots <- list(...)
  if (is.null(dots$xlab)) dots$xlab <- "Time"
  if ((!is.data.frame(x) | !is.matrix(x)) && ncol(x)<2) stop("Wrong type of data")  
  if (ncol(x)==2) {
    if (is.null(dots$curvlab)) {
      causes <- setdiff(unique(x[,2]),0)
      dots$curvlab <- seq(length(causes))
    }
    colnames(x)[1:2] <- c("t","status")
    co <- prodlim::prodlim(Hist(t,status)~1,data=data.frame(x))    
    ##    co <- cuminc(x[,1],x[,2])
    ##    if (any(x[,1]<0)) {
    ##      do.call(co,p)
    ## for (i in seq(length(co))) {
    ##   co[[i]]$time <- co[[i]]$time[-1]
    ##   co[[i]]$est <- co[[i]]$est[-1]
    ## }
    ##if (is.null(dots$xlim)) dots$xlim <- range(x[,1])
    ##    }
    ##do.call("plot", c(list(x=co), dots))
    if (is.null(dots$lwd)) dots$lwd <- 1
    if (missing(lty)) lty <- seq_len(length(co$cuminc))
    if (missing(col)) col <- rep(1,length(co$cuminc))
    if (is.null(dots$lwd)) dots$lwd <- 1
    dots$lty <- lty[1]
    dots$col <- col[1]
    do.call("plot", c(list(x=co), dots))
    dots$add <- TRUE
    for (i in seq_len(length(co$cuminc)-1)+1) {
      dots$lty <- lty[i]
      dots$col <- col[i]
      dots$cause <- i
      do.call("plot", c(list(x=co), dots))
    }
    if (legend)
      legend("topleft",names(co$cuminc),col=col,lty=lty,pch=-1)
    return(invisible(co))
  }

  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  t <- as.matrix(x[,1:2]); cause0 <- as.matrix(x[,3:4])
  causes <- sort(setdiff(unique(cause0),0))
  if (is.null(dots$curvlab))
    dots$curvlab <- seq(1:length(unique(causes)))
  if (missing(col)) col <- c("seagreen","darkred","darkblue","goldenrod","mediumpurple")
  if (1%in%which) {
    plot(t,type="n",...)
    count <- 1
    for (i in causes) {
      points(t[cause0[,1]==causes[i],],col=Col(col[count],0.5),pch=2)
      points(t[cause0[,2]==causes[i],],col=Col(col[count],0.5),pch=6)
      count <- count+1
    }
    points(t[cause0[,1]==0 & cause0[,2]==0,],col=Col("black",0.2),pch=1)  
    if (legend)
    legend("topleft", c("Subj 1, Cause 1", "Subj 2, Cause 1",
                        "Subj 1, Cause 2", "Subj 2, Cause 2",
                        "Double Censoring"), pch=c(2,6,2,6,1), col=c(rep(col[1:length(causes)],each=2),"black"))
  }
  if (2%in%which) {
    dots$curvlab <- NULL
    colnames(x)[1:4] <- c("t1","t2","cause1","cause2")
    co1 <- prodlim::prodlim(Hist(t1,cause1)~1,data.frame(x))
    co2 <- prodlim::prodlim(Hist(t2,cause2)~1,data.frame(x))
    if (is.null(dots$lwd)) dots$lwd <- 1
    if (missing(lty)) lty <- seq_len(length(co1$cuminc))
    if (missing(col)) col <- rep(1,length(co1$cuminc))
    if (is.null(dots$lwd)) dots$lwd <- 1
    dots$lty <- lty[1]
    dots$col <- col[1]
    dots$cause <- cause[1]
    do.call("plot", c(list(x=co1), dots))    
    dots$add <- TRUE
    do.call("plot", c(list(x=co2), dots))
##    for (i in seq_len(length(co1$cuminc)-1)+1) {
    for (i in seq_len(length(cause)-1)+1) {
        dots$lty <- lty[i]
        dots$col <- col[i]
        dots$cause <- cause[i]
        do.call("plot", c(list(x=co1), dots))
        do.call("plot", c(list(x=co2), dots))  
    }
    if (legend && length(cause)>1)
      legend("topleft",names(co1$cuminc)[cause],col=col,lty=lty,pch=-1,bg="white")
  }
}
