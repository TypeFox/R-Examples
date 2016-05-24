{
  rct1 <- 'CC(C)(O)C(O)C(=O)[O-]>>CC(C(C([O-])=O)=O)C.O'
  rct2 <- 'CCC(O)(C)C(C([O-])=O)O>>CCC(C)C(=O)C([O-])=O.[H]O[H]'
  
  message('Reaction definitions:')
  message(paste('Reaction 1:', rct1))
  message(paste('Reaction 2:', rct2))
}

readline("\nType  <Return>\t to start : ")
{
  message('Reaction similarity between rct1 and rct2 (standardize = TRUE): ')
  sim <- rs.compute(rct1, rct2, verbose = TRUE)
  message(sim)
}

readline("\nType  <Return>\t to continue : ")
{
  message('Reaction similarity between rct1 and rct2 (standardize = FALSE): ')
  sim <- rs.compute(rct1, rct2, standardize = FALSE, verbose = TRUE)
  message(sim)
}

readline("\nType  <Return>\t to continue : ")
{
  message('Compound similarity between: O and [H]O[H] (without removal of explicit [H])')
  sim <- ms.compute('O', '[H]O[H]', standardize = FALSE)
  message(sim)
}