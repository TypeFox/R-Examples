plotscore <-
function(param=c(2,.5), fam="pow", bounds, reverse=FALSE, legend=TRUE, ...){
    ## Plot scoring rules (for two-alternative rules only)
    if(length(param) > 2) stop("plotscore is only for two-alternative rules.\n")

    ## For deprecated scaling argument
    dots <- list(...)
    
    if(exists("dots$scaling")){
      if(dots$scaling) bounds <- c(0,1)
    }
    
    p <- seq(.01,.99,.01)
    
    if(missing(bounds)) bounds <- NULL
    sc1 <- calcscore(p, rep(1,length(p)), fam, param, bounds=bounds, reverse=reverse)
    sc0 <- calcscore(p, rep(0,length(p)), fam, param, bounds=bounds, reverse=reverse)

    ymin <- min(sc1,sc0)
    ymax <- max(sc1,sc0)
    
    yl <- c(ymin - .05*(ymax - ymin), ymax + .05*(ymax - ymin))

    main.arg <- list(x=p, y=sc1)
    supplied <- list(...)
    default <- list(type="l", ylim=yl, xlab="Forecast", ylab="Score")
    nomatch <- setdiff(c("type","xlab","ylab","ylim"), names(supplied))
    plot.args <- c(main.arg, supplied, default[nomatch])

    do.call(plot, plot.args)
    lines(p, sc0, lty=2)
    if(legend) legend(.8, yl[2] - .1, c("d=1","d=0"), lty=c(1,2))
}
