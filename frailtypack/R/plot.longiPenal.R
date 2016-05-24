"plot.longiPenal" <- function (x, type.plot="Hazard", conf.bands=TRUE, pos.legend="topright", cex.legend=0.7, main, color=2, ...)
{

  plot.type <- charmatch(type.plot, c("Hazard", "Survival"),nomatch = 0)
  if (plot.type == 0) {
    stop("estimator must be Hazard or Survival")
  }

  if(missing(main))
    main<-""

  if(plot.type==1){ # hazard

    if(conf.bands){
      matplot(x$xD[-1,1], x$lamD[-1,,1], col=color, type="l", lty=c(1,2,2), xlab="Time",ylab="Hazard function", main=main, ...)
    }else{
      plot(x$xD[-1,1], x$lamD[-1,1,1], col=color, type="l", lty=1, xlab="Time",ylab="Hazard function", main=main, ...)
    }

  }else{ # survival

    if (x$typeof == 0){
      if (conf.bands){
        matplot(x$xD[,1], x$survD[,,1], col=color, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", main=main, ...)
     }else{
        plot(x$xD[,1], x$survD[,1,1], col=color, type="l", lty=1, xlab="Time",ylab="Baseline survival function", main=main, ...)
       }
    }else{
      if (conf.bands){
        matplot(x$xSuD[,1], x$survD[,,1], col=color, type="l", lty=c(1,2,2), xlab="Time",ylab="Baseline survival function", main=main, ...)
       }else{
        plot(x$xSuD[,1], x$survD[,1,1], col=color, type="l", lty=1, xlab="Time",ylab="Baseline survival function", main=main, ...)
      }
    }

  }

  legend(pos.legend, c("event"), lty=1, col=color, cex=cex.legend, ...)

  return(invisible())
}
