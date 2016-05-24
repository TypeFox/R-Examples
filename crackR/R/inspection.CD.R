inspection.CD <-
function(obj, inspection.type=1)
{
  ## function for conducting a future inspection with continuing damage involving two cracks
  ## cracks may be at the same hole or at two holes, according to parameter two.hole.bool
  ## only one repair type for this version, but use POD.threshold for reporting sizes
  ##   at detection and repair

  ## there's currently no bootstrapping with the CD code...do i need that?
  
  a.p <- obj$state$a.p
  a.s <- obj$state$a.s
  kc  <- obj$state$kc
  w   <- obj$state$w

  if(any(is.na(w)))
  {
    stop("at least one weight is NaN. this usually means all particles have reached the critical crack length.")
  }

  cg.cc.pc <- obj$parameters$cg.cc.pc
  cg.cc.sc <- obj$parameters$cg.cc.sc
  cg.cc.ph <- obj$parameters$cg.cc.ph
  cg.cc.sh <- obj$parameters$cg.cc.sh
  
  dta.pc <- obj$parameters$dta.pc
  dta.ph <- obj$parameters$dta.ph
  dta.sc <- obj$parameters$dta.sc
  dta.sh <- obj$parameters$dta.sh  

  Np <- obj$parameters$Np

  pod.fun <- obj$parameters$pod.func[[inspection.type]]

  ## there's only one type, use new parameter specific to one type model
  ## this is same as in RBDMS; values are the dividing points for repair types
  pod.threshold.1typ <- obj$parameters$pod.threshold.1typ

  ## if two holes, either crack may be repaired
  ## if one hole, both will be if either is found
  two.hole.bool <- obj$parameters$two.hole.bool

  ## remove particles with zero weight since we don't want to re-sample these
  ## this used to check the crack sizes against cg.cc.ph/sh, but that's not the issue
  living     <- ( w > 0 )
  living.num <- sum(living)
  
  a.p <- a.p[living]
  a.s <- a.s[living]
  kc  <- kc[living]
  w   <- w[living]
  
  ## probability of detection for each crack (independent)
  pod.p <- pod.fun(a.p)
  pod.s <- pod.fun(a.s)

  ## yields repair type for each particle (if each were found)
  rep.types.p <- findInterval(a.p, pod.threshold.1typ) + 1
  rep.types.s <- findInterval(a.s, pod.threshold.1typ) + 1

  ## currently treating the two hole and one hole cases separately
  if( two.hole.bool )
  {
      ## code for two holes, each with a crack

      ## new input parameter identifies which repair types fix both cracks
      two.hole.repair.both <- obj$parameters$two.hole.repair.both
      two.hole.repair.both <- c(FALSE, two.hole.repair.both) ## never repair both for smallest type repair

      ## there are four possible outcomes for repair: none, p.only, s.only, both
      ##   partition the weights of each particle according to the four outcomes

      pod.none <- ( 1 - pod.p ) * ( 1 - pod.s )

      pod.p.only <- pod.p * ( 1 - pod.s )
      pod.s.only <- pod.s * ( 1 - pod.p )

      ## note, if a.p is large enough to repair both cracks, p.only is zero (and vice versa)
      ##   first find which repair types will repair both at once (x1)
      ##   if any, loop over those repair types and get an index for which particles are repaired together
      ## rep.both.x1.index indicates which pairings will involve one repair on both (used later)
      rep.both.x1.test <- which( two.hole.repair.both )
      rep.both.x1.index <- rep(FALSE, living.num)
      if( length( rep.both.x1.test ) > 0 )
      {
          for(kkk in rep.both.x1.test)
              rep.both.x1.index[ ( rep.types.p == kkk ) | ( rep.types.s == kkk ) ] <- TRUE
          pod.p.only[ rep.both.x1.index ] <- 0
          pod.s.only[ rep.both.x1.index ] <- 0
      }

      ## pod.both must equal the remaining probability for each particle
      pod.both <- 1 - pod.none - pod.p.only - pod.s.only

      ## apportion the weights accordingly
      w.none   <- w * pod.none
      w.p.only <- w * pod.p.only
      w.s.only <- w * pod.s.only
      w.both   <- w * pod.both
      
      ## the above is correct as far as performing repairs is concerned
      ## however, when reporting PCD, part of pod.both should be double-counted
      pod.both.x2 <- pod.both
      pod.both.x2[ rep.both.x1.index ] <- 0
      w.both.x2     <- w * pod.both.x2
      w.both.x2.sum <- sum(w.both.x2)

      ## similar for the both repaired together case
      pod.both.x1 <- pod.both
      pod.both.x1[ !rep.both.x1.index ] <- 0
      w.both.x1     <- w * pod.both.x1
      w.both.x1.sum <- sum(w.both.x1)

      ## partition these by repair type for output
      n.rep.types <- length(pod.threshold.1typ) + 1
      pcd.by.type <- rep(0, n.rep.types)
      rep.types.larger <- apply(cbind(rep.types.p,rep.types.s), 1, max)
      for( iii in 1:n.rep.types )
      {
          pcd.by.type[iii] <- pcd.by.type[iii] + sum( w.p.only[ rep.types.p    == iii ] )
          pcd.by.type[iii] <- pcd.by.type[iii] + sum( w.s.only[ rep.types.s    == iii ] )

          ## both repaired together (repair type must be the larger since repairing both)
          pcd.by.type[iii] <- pcd.by.type[iii] + sum( w.both.x1[ rep.types.larger == iii ] )

          ## for the both repaired separately cases, there are two repairs
          pcd.by.type[iii] <- pcd.by.type[iii] + sum( w.both.x2[ rep.types.p == iii ] )
          pcd.by.type[iii] <- pcd.by.type[iii] + sum( w.both.x2[ rep.types.s == iii ] ) 
      }

      #############
      ## moving on to doing the repairs
      ## at this point we don't care about the distinction between x1 and x2

      ## number of unmodified particles would ideally be ~ Np*sum(w.none)
      ## number of p rep only particles would ideally be ~ Np*sum(w.p.only)
      ## number of s rep only particles would ideally be ~ Np*sum(w.s.only)
      ## number of both rep'd particles would ideally be ~ Np*sum(w.both)

      w.none.sum   <- sum( w.none   )
      w.p.only.sum <- sum( w.p.only )
      w.s.only.sum <- sum( w.s.only )
      w.both.sum   <- sum( w.both   )
      
      Np.none   <- round( Np * w.none.sum   )
      Np.p.only <- round( Np * w.p.only.sum )
      Np.s.only <- round( Np * w.s.only.sum )
      Np.both   <- round( Np * w.both.sum   )

      ## minimum number of particles in each group
      Np.min <- 10

      ## need to set Np partition according to these rules:
      ##   sum to Np
      ##   at least Np.min in each outcome
      ##   none, p.only, and s.only must have <= living.num in each
      ##   after sampling, at least 3 particles in each outcome need
      ##     to have positive weight (o.w. zero all weights)

      if( Np.none   > living.num ) Np.none   <- living.num
      if( Np.p.only > living.num ) Np.p.only <- living.num
      if( Np.s.only > living.num ) Np.s.only <- living.num

      if( Np.none   < Np.min ) Np.none   <- Np.min
      if( Np.p.only < Np.min ) Np.p.only <- Np.min
      if( Np.s.only < Np.min ) Np.s.only <- Np.min

      ## two possibilities:
      ##   1. sum leaves enough room for Np.both; i.e., sum < Np-Np.min
      ##      in this case, set Np.both to remainder and we're done
      ##   2. sum does not leave enough room for Np.both
      ##      in this case, find and reduce the largest bin to accommodate
      if( sum(Np.none,Np.p.only,Np.s.only) <= (Np-Np.min) )
          {
              Np.both <- Np - Np.none - Np.p.only - Np.s.only
          } else {
              temp <- which.max( c(Np.none, Np.p.only, Np.s.only) )
              if( temp == 1 ) Np.none   <- Np - Np.min - Np.p.only - Np.s.only
              if( temp == 2 ) Np.p.only <- Np - Np.min - Np.none   - Np.s.only
              if( temp == 3 ) Np.s.only <- Np - Np.min - Np.none   - Np.p.only
              Np.both <- Np - Np.none - Np.p.only - Np.s.only
          }

      ## sample for outcome none
      index.none <- sample(1:living.num, size=Np.none, replace=FALSE)
      a.p.none <- a.p[ index.none ]
      a.s.none <- a.s[ index.none ]
      kc.none  <- kc[ index.none ]
      w.none   <- w.none[ index.none ]
      w.none   <- w.none.sum * w.none / sum( w.none )

      ## sample for outcome p.only
      ## XXX if monte carlo does not match particle filter...w.p.only weighting could be the issue
      index.p.only <- sample(1:living.num, size=Np.p.only, replace=FALSE)
      a.p.p.only <- dta.pc$rfs.rsamp( Np.p.only ) ## XXX 1typ model needs repair eifs !!!
      a.s.p.only <- a.s[ index.p.only ]
      kc.p.only  <- kc[ index.p.only ]
      w.p.only   <- w.p.only[ index.p.only ] * dta.pc$rfs.dactual(a.p.p.only) / dta.pc$rfs.dsamp(a.p.p.only)
      w.p.only   <- w.p.only.sum * w.p.only / sum( w.p.only )

      ## sample for outcome s.only
      ## XXX if monte carlo does not match particle filter...w.s.only weighting could be the issue
      index.s.only <- sample(1:living.num, size=Np.s.only, replace=FALSE)
      a.p.s.only <- a.p[ index.s.only ]
      a.s.s.only <- dta.sc$rfs.rsamp( Np.s.only ) ## XXX 1typ model needs repair eifs !!!
      kc.s.only  <- kc[ index.s.only ]
      w.s.only   <- w.s.only[ index.s.only ] * dta.sc$rfs.dactual(a.s.s.only) / dta.sc$rfs.dsamp(a.s.s.only)
      w.s.only   <- w.s.only.sum * w.s.only / sum( w.s.only )

      ## sample for outcome both
      a.p.both <- dta.pc$rfs.rsamp( Np.both ) ## XXX 1typ model needs repair eifs !!!
      a.s.both <- dta.sc$rfs.rsamp( Np.both ) ## XXX 1typ model needs repair eifs !!!
      kc.both  <- dta.pc$kc.rsamp( Np.both )  ## only when both repaired is new Kc value (pri and sec same)
      w.both   <- ( dta.pc$rfs.dactual(a.p.both) * dta.sc$rfs.dactual(a.s.both) ) /
          ( dta.pc$rfs.dsamp(a.p.both) * dta.sc$rfs.dsamp(a.s.both) )
      w.both   <- w.both.sum * w.both / sum( w.both )

      ## verify that weights for several particles in each set are > 0
      ## if not, zero the outcome weights entirely
      if( sum( w.none   > 0 ) < 3 ) w.none   <- rep(0, Np.none)
      if( sum( w.p.only > 0 ) < 3 ) w.p.only <- rep(0, Np.p.only)
      if( sum( w.s.only > 0 ) < 3 ) w.s.only <- rep(0, Np.s.only)
      if( sum( w.both   > 0 ) < 3 ) w.both   <- rep(0, Np.both)

      ## collect the four outcomes into a single set
      obj$state$a.p <- c(a.p.none, a.p.p.only, a.p.s.only, a.p.both)
      obj$state$a.s <- c(a.s.none, a.s.p.only, a.s.s.only, a.s.both)
      obj$state$w   <- c(w.none, w.p.only, w.s.only, w.both)
      obj$state$kc  <- c(kc.none, kc.p.only, kc.s.only, kc.both)

      names(pcd.by.type) <- paste("pcd", 1:length(pcd.by.type), sep="")

      pcd.new <- data.frame(flight=max(obj$results$sfpof$flight),
                            insp.type=inspection.type,
                            t(pcd.by.type),
                            pcd=sum(pcd.by.type),
                            cd.none = w.none.sum,
                            cd.p.only = w.p.only.sum,
                            cd.s.only = w.s.only.sum,
                            cd.both.x1 = w.both.x1.sum,
                            cd.both.x2 = w.both.x2.sum)
      
      obj$results$pcd <- rbind(obj$results$pcd, pcd.new)

  } else {
      ## code for one hole with two cracks

      ## for one hole, we need the probability of finding either crack
      pod.either <- pod.p + pod.s - pod.p * pod.s

      pcd <- pod.either * w

      ## assume for cost purposes the more severe repair would occur
      rep.types <- apply( cbind(rep.types.p, rep.types.s), 1, max )
      
      ## partition PCD into repair types
      n.rep.types   <- length(pod.threshold.1typ) + 1
      pcd.by.type <- rep(0, n.rep.types)
      for( iii in 1:n.rep.types )
          pcd.by.type[iii] <- sum( pcd[ rep.types == iii ] )

      ## there are TWO possible outcomes for repair: none, both
      ## partition the weights of each particle according to the TWO outcomes
      w.none   <- w * ( 1 - pod.either )
      w.both   <- w * pod.either

      ## number of unmodified particles would ideally be ~ Np*sum(w.none)
      ## number of both rep'd particles would ideally be ~ Np*sum(w.both)

      w.none.sum   <- sum( w.none   )
      w.both.sum   <- sum( w.both   )
      
      Np.none   <- round( Np * w.none.sum   )
      Np.both   <- round( Np * w.both.sum   )

      ## minimum number of particles in each group
      Np.min <- 10

      ## need to set Np partition according to these rules:
      ##   sum to Np
      ##   at least Np.min in each outcome
      ##   none, p.only, and s.only must have <= living.num in each
      ##   after sampling, at least 3 particles in each outcome need
      ##     to have positive weight (o.w. zero all weights)

      if( Np.none   > living.num ) Np.none <- living.num
      if( Np.none   < Np.min )     Np.none <- Np.min

      ## two possibilities:
      ##   1. sum leaves enough room for Np.both; i.e., sum < Np-Np.min
      ##      in this case, set Np.both to remainder and we're done
      ##   2. sum does not leave enough room for Np.both
      ##      in this case, reduce none to accommodate
      if( Np.none <= (Np-Np.min) )
          {
              Np.both <- Np - Np.none
          } else {
              Np.none <- Np - Np.min
              Np.both <- Np - Np.none
          }

      ## sample for outcome none
      index.none <- sample(1:living.num, size=Np.none, replace=FALSE)
      a.p.none <- a.p[ index.none ]
      a.s.none <- a.s[ index.none ]
      kc.none  <- kc[ index.none ]
      w.none   <- w.none[ index.none ]
      w.none   <- w.none.sum * w.none / sum( w.none )

      ## sample for outcome both
      a.p.both <- dta.pc$rfs.rsamp( Np.both )
      a.s.both <- dta.sc$rfs.rsamp( Np.both )
      kc.both  <- dta.pc$kc.rsamp( Np.both )
      w.both   <- ( dta.pc$rfs.dactual(a.p.both) * dta.sc$rfs.dactual(a.s.both) ) /
          ( dta.pc$rfs.dsamp(a.p.both) * dta.sc$rfs.dsamp(a.s.both) )
      w.both   <- w.both.sum * w.both / sum( w.both )

      ## verify that weights for several particles in each set are > 0
      ## if not, zero the outcome weights entirely
      if( sum( w.none > 0 ) < 3 ) w.none   <- rep(0, Np.none)
      if( sum( w.both > 0 ) < 3 ) w.both   <- rep(0, Np.both)

      ## collect the two outcomes into a single set
      obj$state$a.p <- c(a.p.none, a.p.both)
      obj$state$a.s <- c(a.s.none, a.s.both)
      obj$state$w   <- c(w.none,   w.both)
      obj$state$kc  <- c(kc.none,  kc.both)

      names(pcd.by.type) <- paste("pcd", 1:length(pcd.by.type), sep="")

      pcd.new <- data.frame(flight=max(obj$results$sfpof$flight),
                            insp.type=inspection.type,
                            t(pcd.by.type),
                            pcd=sum(pcd.by.type))
      
      obj$results$pcd <- rbind(obj$results$pcd, pcd.new)
  }
  
  return(obj)
}
