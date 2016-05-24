ucipBound <- function(RT, CR, stopping.rule=c("AND", "OR"), bound=c("upper", "lower"), Cspace =FALSE) {
  tvec <- c(0,sort(unique(unlist(RT))))
  returnlist <- vector("list", 0)
  if (agrep("upper", bound,  ignore.case=TRUE ) ) {
    if(agrep("or", stopping.rule, ignore.case=TRUE) ){ 

      if(Cspace) {
        numer <- 1-ecdf( RT[[2]][ CR[[2]]==1 ] )(tvec) +  1-ecdf( RT[[3]][ CR[[3]]==1 ] )(tvec) -1
        denom <- (1-ecdf( RT[[2]][ CR[[2]]==1 ] )(tvec) ) * (1-ecdf( RT[[3]][ CR[[3]]==1 ] )(tvec) )
        fn <- stepfun(tvec, c(log(numer) / log(denom), NA ))
        #  NEED TO FIGURE OUT PLOTTING LIMITS
        attr(fn, "call") <- "log(S1(t) + S2(t)-1)/log(S1(t)S2(t)"
      } else {
        fn <- stepfun(tvec, c(ecdf( RT[[2]][ CR[[2]]==1 ] )(tvec) +  ecdf( RT[[3]][ CR[[3]]==1 ] )(tvec), 2))
        attr(fn, "call") <- "F1(t) + F2(t)"
      }
      returnlist <- c(returnlist, upper.or = fn)
    }
    if(agrep("and", stopping.rule, ignore.case=TRUE) ) { 
      if(Cspace) {
      } else { 
        fn <- stepfun(tvec, c(pmin( ecdf( RT[[2]][ CR[[2]]==1 ] )(tvec),  ecdf( RT[[3]][ CR[[3]]==1 ] )(tvec) ),1) )
        attr(fn, "call") <- "MIN(F1(t), F2(t))"
      }
      returnlist <- c(returnlist, upper.and = fn)
    }
  }

  if (agrep("lower", bound, ignore.case=TRUE ) ) {
    if(agrep("or", stopping.rule, ignore.case=TRUE) ){ 
      if(Cspace) {
      } else { 
        fn <- stepfun(tvec, c(pmax( ecdf( RT[[2]][ CR[[2]]==1 ] )(tvec),  ecdf( RT[[3]][ CR[[3]]==1 ] )(tvec) ),1) )
        attr(fn, "call") <- "MAX(F1(t), F2(t))"
      }
      returnlist <- c(returnlist, lower.or = fn)
    }
    if(agrep("and", stopping.rule, ignore.case=TRUE)) { 
      if(Cspace) {
      } else { 
        fn <- stepfun(tvec, c(ecdf( RT[[2]][ CR[[2]]==1 ] )(tvec) +  ecdf( RT[[3]][ CR[[3]]==1 ] )(tvec) - 1, 1))
        attr(fn, "call") <- "F1(t) + F2(t)-1"
      }
      returnlist <- c(returnlist, lower.and = fn)
    }
  }

  return(returnlist)
}
