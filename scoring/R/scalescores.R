scalescores <-
function(pars, fam="pow", ordered, nalts){
  ps <- c(0,1)

  ## Shortcut for scaling scores for >2 alts
  ## Min must be achieved when forecast 1 for lowest baseline
  ##     and lowest baseline occurs.
  ## Max must be achieved when forecast 0 for highest baseline
  ##     and 1 for lowest and highest baseline occurs.

  ## No ordered mods because this family can only be used
  ## for two-alternative forecasts
  if(fam=="beta"){
      xplier <- max(betafam(c(1,0), d=2, param=pars),
                    betafam(c(0,1), d=1, param=pars))

      if(xplier == Inf){
          warning("Scaling does not work because maximum possible score is Inf.")
          xplier <- 1
      }
      
  } else if(fam=="pow" | fam=="sph"){
      ## Check for baseline parameters
      if(length(pars) == 1){
        maxbase <- 1
        minbase <- 2
      } else {
          baselines <- pars[2:(nalts+1)]
          ## Which are largest and smallest?
          maxbase <- which(baselines==max(baselines))
          if(length(maxbase)>1){
              ## Choose the most extreme alternative
              maxbase <- maxbase[which.max(abs(maxbase - length(baselines)/2))[1]]
          }
          minbase <- which(baselines==min(baselines))
          if(length(minbase)>1){
              ## Choose the most extreme alternative
              extcat <- minbase[which.max(abs(minbase - length(baselines)/2))[1]]
              ## In case they are all equal, take a different category
              minbase <- ifelse(extcat == maxbase,
                                minbase[minbase!=maxbase][1],
                                extcat)
          }
      }

      if(ordered==FALSE){
          ## We can easily identify the min/max score
          fore <- out <- rep(0,nalts)
          fore[minbase] <- 1

          tmpsc <- calcscore(c(minbase,maxbase) ~ rbind(fore,fore), fam=fam, param=pars, ordered=ordered)
      
          xplier <- tmpsc #c(scmin, scmax)
      } else {
          ## Try a larger number of values (could probably derive
          ## the min/max to make this faster)
          fore <- diag(1, nalts)
          minbase <- rep(minbase, nrow(fore))
          maxbase <- rep(maxbase, nrow(fore))
          tmpsc <- calcscore(c(minbase, maxbase) ~ rbind(fore, fore), fam=fam, param=pars, ordered=ordered)

          xplier <- c(min(tmpsc), max(tmpsc))
      }
      ## Will get Inf for log scores and possibly others
      if(xplier[2]==Inf){
          warning("Scaling does not work because maximum possible score is Inf.")
          xplier <- c(0,1)
      }
  }

  xplier
}

if(FALSE){
    ## Proof that above yields min and max values for 3 alts:
    p <- seq(0,1,.01)
    y <- expand.grid(p,p,p)
    ysum <- apply(y,1,sum)
    y <- y[ysum==1,]
    out1 <- t(matrix(c(1,0,0),3,nrow(y)))
    out2 <- t(matrix(c(0,1,0),3,nrow(y)))
    out3 <- t(matrix(c(0,0,1),3,nrow(y)))

    sc1 <- calcscore(y ~ out1,fam="pow",param=c(3,.1,.5,.4))
    sc2 <- calcscore(y ~ out2,fam="pow",param=c(3,.1,.5,.4))
    sc3 <- calcscore(y ~ out3,fam="pow",param=c(3,.1,.5,.4))

    c(max(sc1),max(sc2),max(sc3))
    y[which(sc2==max(sc2)),]
    c(min(sc1),min(sc2),min(sc3))
    y[which(sc1==min(sc1)),]
}
