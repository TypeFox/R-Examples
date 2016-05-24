calcSol <- function(lat, BTd, sample='hour', BTi, EoT=TRUE, keep.night=TRUE, method='michalsky'){

  if (missing(BTi)){
    solD<-fSolD(lat, BTd=BTd, method=method)
    solI<-fSolI(solD, sample=sample, EoT=EoT, keep.night=keep.night, method=method)
    match <- attr(solI, 'match')
    sample <- attr(solI, 'sample')
  } else { ##utilizo BTi
    BTd=unique(truncDay(BTi))
    solD <- fSolD(lat, BTd=BTd, method=method)
    solI <- fSolI(solD, BTi=BTi, EoT=EoT, keep.night=keep.night, method=method)
    match <- attr(solI, 'match')
    sample <- attr(solI, 'sample')
  }
  attr(solD, 'lat') <- NULL
  attr(solI, 'lat') <- NULL
  attr(solI, 'match') <- NULL
  attr(solI, 'sample') <- NULL
  result <- new('Sol',
                lat=lat,
                solD=solD,
                solI=solI,
                match=match,
                sample=sample,
                method=method)
  return(result)
}
