setAs("MCEstimate", "mle", def = function(from){
       crit.f0 <- from@criterion.fct
       start.f0 <- as.list(from@untransformed.estimate)
       if(!.isUnitMatrix(trafo(from)$mat)){
          
          ### we have to turn crit.f0 into a function in the
          ### transformed parameter; to this end, we specify
          ### an argument list argList.f0 with all transformed
          ### parameter coordinates

          ## store MLestimates
          est0 <- from@estimate
          est1 <- from@untransformed.estimate

          
          argList.f0 <- c(from@estimate, from@fixed)
          l <- length(argList.f0)
          crit.lst <- vector("list", l+1)
          crit.lst[1:l] <- argList.f0
          names(crit.lst) <- c(names(argList.f0),"")
          nx <- from@nuis.idx

          ft <-function(){          
                ## stack list of arguments into one vector th0
                ## (in transformation range) 
                mc <- as.list(match.call())[-1]
                th0 <- argList.f0
                th0[names(mc)] <- mc
                th0 <- c(unlist(th0))
       
                ## partition it into main, nuisance
                est.main <- est <- th0
                est.nuis <- NULL
                if(length(nx)){
                   est.main <- est[-nx]
                   est.nuis <- est[nx]
                }
                ## generate a valid ParamFamParameter object out of it
                param <- ParamFamParameter(main = est.main, nuisance = est.nuis,
                                           fixed = from@fixed)

                ## "invert" (locally!) the transformation,
                # i.e. th1 "=" trafo^-1(th0)                
                D1 <- (trafo(from)$fct)(th0)$mat
                th1 <- est1 + solve(D1, th0-est0)
                ## call critrium.fct with this transformed parameter
                do.call(from@criterion.fct,as.list(th1))
          }
          crit.lst[l+1] <- as.list(ft)[1]
          ## crit.f0 is now a function in the transformed parameter
          crit.f0 <- as.function(crit.lst)
          start.f0 <- as.list(from@estimate)
      }
      to <- new("mle")
      to@call <- substitute(mle(minuslogl = crit.f, start = startPar), 
                            list(crit.f = crit.f0,
                                 startPar = start.f0))
      to@coef <- from@estimate
      fe <- if(is.null(from@untransformed.estimate))
               from@estimate else from@untransformed.estimate
      to@fullcoef <- c(fe,from@fixed)
      to@vcov <- if(!is.null(from@asvar)) 
                 from@asvar/from@samplesize else matrix(NA,1,1)
      to@min <- from@criterion
      to@details <- as.list(c(from@Infos))
      to@method <- from@method
      to@minuslogl <- crit.f0
to})

setMethod("profile", "MCEstimate",
          function (fitted, which = 1:length(fitted@estimate), maxsteps = 100,
                    alpha = 0.01, zmax = sqrt(qchisq(1 - alpha, 1L)),
                    del = zmax/5, trace = FALSE, ...){
m.mle <- as(fitted,"mle")
profile(m.mle, which=which, maxsteps=maxsteps, alpha=alpha, zmax=zmax,
del=del, trace=trace, ...)
})
