`dumpData` <-
function(seg.ratio, model, stem="test", fix.one=TRUE,
                     data.file=paste(stem,"-data.R",sep="") )
{

  if (class(seg.ratio) != "segRatio") {
    stop("'seg.ratio' must be of class 'segRatio'")
  }

  ## set some values of T so that each chain always has at least one
  ## member the easiest is to force value(s) closest to expected
  ## logit(ratio) to be in the appropriate dosage class

  n.comp <- model$n.components
  n.indiv <- length(seg.ratio$n)
  E <- model$E.segRatio$ratio[1:n.comp]
  
  T <- rep(NA,n.indiv)

  if (fix.one) {
    for (i in 1:n.comp){ # find closest to expected segregation ratios
      T[ abs(seg.ratio$seg.ratio - E[i]) ==
        min(abs(seg.ratio$seg.ratio - E[i])) ] <- i
      if ((neqi <- sum(T==i, na.rm=TRUE))>1) { # randomly drop all except 1 
        T[ sample((1:length(T))[T==i & ! is.na(T)],neqi-1) ] <-  NA
      }
    }
  }

  N <- n.indiv
  r <- seg.ratio$r
  n <- seg.ratio$n

  ## NB:  control="S_compatible", added for JAGS/R > 2.50 compatability
  dump(c("N","r","n","T"),  control="S_compatible", file=data.file)

  ## superseded now that JAGS Version 1.0 required
  ##if (.Platform$OS.type == "windows"){# to be fixed after JAGS 0.90 superseded
  ##  fix.up <- readLines(data.file)
  ##  fixed <- gsub("`","\"",fix.up)  # replace ` with "
  ##  writeLines(fixed, data.file)
  ##}
  
}

