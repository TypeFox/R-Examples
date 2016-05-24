"drift" <-
function(za=-zb,zb,t,t2=t,pow=NULL,drft=NULL,conf=NULL,zval=zb[length(zb)]){
   if ({length(za) != length(zb)}|{length(t) != length(zb)}|{length(t) != length(t2)}){
      stop("Analysis times, upper bounds, and lower bounds must all be vectors of the same length.")
    }
   if ({min(t) <= 0}|{max(t) > 1}|{min(t2) <= 0}){
      stop("Analysis times must be in (0,1].  Second time scale values must be positive.")
    }
   if (sum({t-c(0,t[-length(t)]) <= 0}|{t2-c(0,t[-length(t2)]) <= 0}) > 0){
      stop("Analysis times must be ordered from smallest to largest.")
    }
   t3 <- t2
   t2 <- t2/max(t2)
   if ({{!is.null(pow)}&{!is.null(drft)}}|{{!is.null(pow)}&{!is.null(conf)}}|
       {{!is.null(conf)}&{!is.null(drft)}}){
      stop("Only one of power, drift, and confidence level can be given.")
    }
   drift1 <- NULL
   if (!is.null(pow)){
      if ({pow <= 0}|{pow > 1}){
         stop("Power must be in (0,1].")
       }
      type <- 1
      drift1 <- adrift(t2,za,zb,pow)
    }
   if (!is.null(drft)){
      type <- 2
      drift1 <- drft
    }
   if (!is.null(drift1)){
      gl <- glan(t2,za,zb,drift1)
      if (!is.null(drft)) pow <- gl$pr
      ans <- list(type=type,time=t,time2=t3,lower.bounds=za,upper.bounds=zb,power=pow,
                  drift=drift1,lower.probs=gl$qneg,upper.probs=gl$qpos,
                  exit.probs=gl$qneg+gl$qpos,cum.exit=cumsum(gl$qneg+gl$qpos))
    }
   if (!is.null(conf)){
      if (zval < 0){
         stop("Confidence interval is only for nonnegative final Z value.")
       }
      conf.limit <- ci(conf,zval,t2,za,zb)
      ans <- list(type=3,time=t,time2=t3,lower.bounds=za,upper.bounds=zb,
                  conf.level=conf,final.zvalue=zval,conf.interval=conf.limit)
    }
   class(ans) <- "drift"
   return(ans)
 }

