{
  rct1 <- 'CC(C)(O)C(O)C(=O)[O-]>>CC(C(C([O-])=O)=O)C.O'
  rct2 <- 'CCC(O)(C)C(C([O-])=O)O>>CCC(C)C(=O)C([O-])=O'
  
  message('Reaction definitions:')
  message(paste('Reaction 1:', rct1))
  message(paste('Reaction 2:', rct2))
}

readline("\nType  <Return>\t to start : ")
{
  message('Reaction similarity between rct1 and rct2 (algo = msim): ')
  message('Adds 0 similarity score for unpaired O in Rct1.')
  sim <- rs.compute(rct1, rct2, verbose = TRUE)
  message(sim)
}

readline("\nType  <Return>\t to start : ")
{
  message('Reaction similarity between rct1 and rct2 using algo msim_max: ')
  sim <- rs.compute(rct1, rct2, verbose = TRUE, algo = 'msim_max')
  message(sim)
}