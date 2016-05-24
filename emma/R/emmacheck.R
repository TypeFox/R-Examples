emmacheck  <- function(x, graph, fn1 = NULL, fn2 = NULL, fn3 = NULL, fn4 = NULL, nresp=1)
{
  x$yspace <- estimateModel(x,graph,nresp=nresp)
  x$yspace[x$tested,] <- as.data.frame(x$ypop)

  d <- distance(x$xpop,x$xspace,x$yspace,x$weight,x$opt)

  Gb.pred <- which.min(d$fit)

  if(! Gb.pred %in% x$tested){
    x$add <- 1
    print("PERFORM THE FOLLOWING ADDITIONAL EXPERIMENT:")
    print(x$xspace[Gb.pred,])
    x$xpop <- rbind(x$xpop,x$xspace[Gb.pred,])
    x$tested <- c(x$tested,Gb.pred)
    x$ypop <- rbind(x$ypop,rep(0,ncol(x$ypop)))

    if (!is.null(fn1)) x$ypop[nrow(x$ypop),1] <- fn1(x$xspace[Gb.pred,])
    if (!is.null(fn2)) x$ypop[nrow(x$ypop),2] <- fn2(x$xspace[Gb.pred,])
    if (!is.null(fn3)) x$ypop[nrow(x$ypop),3] <- fn3(x$xspace[Gb.pred,])
    if (!is.null(fn4)) x$ypop[nrow(x$ypop),4] <- fn4(x$xspace[Gb.pred,])

    if (is.null(fn1)) { 
      ypop <- NULL
      print("ADD THE MEASURED RESPONSE VALUES TO ypop")
    }
  }

return(x)
}

