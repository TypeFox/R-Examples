 
  mattools.roc  <-  function(inFossil, inModern, modTaxa=c(), fosTaxa=modTaxa, numAnalogs=NULL)
  {
    nColsMatrix = length(inFossil[, 1])
    LocMinRow = matrix(NA, nrow = numAnalogs, ncol = nColsMatrix)
    modMatrix = as.matrix(inModern[, modTaxa])
    fossilMatrix = sqrt(as.matrix(inFossil[, fosTaxa]))
    dimnames(modMatrix) = NULL
    dimnames(fossilMatrix) = NULL
    modMatrix = sqrt(t(modMatrix))
    for(i in 1:length(inFossil[, 1])) {
      currSpectrum = fossilMatrix[i,  ]
      sqdistVec = (currSpectrum - modMatrix)
      sqdistVec= sqdistVec*sqdistVec
      x = colSums(sqdistVec)
      y = rank(x,ties.method="first")
      #x = sort(x[which(y <= numAnalogs)])
      x = sort(x[y <= numAnalogs])
      LocMinRow[, i] = x
    }
    LocMinRow
  }