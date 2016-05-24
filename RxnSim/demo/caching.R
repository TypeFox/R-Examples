{
  rct1 <- 'CC(C)(O)C(O)C(=O)[O-]>>CC(C(C([O-])=O)=O)C.[OH2]'
  rct2 <- 'CCC(O)(C)C(C([O-])=O)O>>CCC(C)C(=O)C([O-])=O.[OH2]'

  message('Reaction definitions:')
  message(paste('Reaction 1:', rct1))
  message(paste('Reaction 2:', rct2))
  message(paste('Similarity between rct1 and rct2: ', rs.compute(rct1, rct2)))
}

readline("\nType  <Return>\t to start : ")
{
  message('Time to compute similarity without fingerprint caching (200 iterations).')
  message(signif(system.time(
    for (i in 1:200) {
      rs.compute(rct1, rct1, fpCached = FALSE)
    }
  )[3]), 1)
}

readline("\nType  <Return>\t to continue : ")
{
  message('Time to compute similarity with fingerprint caching (200 iterations).')
  rs.clearCache()
  message(signif(system.time(
    for (i in 1:200) {
      rs.compute(rct1, rct1, fpCached = TRUE)
    }
  )[3]), 1)
}

readline("\nType  <Return>\t to continue : ")
{
  message('Fingerprint generated and cached using circular fingerprint type.')
  rs.clearCache()
  sim <- ms.compute('CCCCO', 'CCOCCC', fp.type = 'circular', fpCached = T)
  message(paste('Similarity:', sim))
}

readline("\nType  <Return>\t to continue : ")
{
  message('Similarity value using default fingerprint type')
  sim <- ms.compute('CCCCO', 'CCOCCC')
  message(paste('Similarity:', sim))
}

readline("\nType  <Return>\t to continue : ")
{
  message('Similarity value obtained by using cached fingerprints')
  sim <- ms.compute('CCCCO', 'CCOCCC', fpCached = T)
  message(paste('Similarity:', sim))
}

readline("\nType  <Return>\t to continue : ")
{
  message('Similarity post clearing cache.')
  message('Fingerprints are re-generated and cached using default fingerprint type')
  rs.clearCache()
  sim <- ms.compute('CCCCO', 'CCOCCC', fpCached = T)
  message(paste('Similarity:', sim))
}