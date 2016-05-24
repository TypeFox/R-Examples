{
  rct1 <- 'CC(C)(O)C(O)C(=O)[O-]>>CC(C(C([O-])=O)=O)C.O'
  rct2 <- 'CCC(O)(C)C(C([O-])=O)O>>CCC(C)C(=O)C([O-])=O.O'
  
  message('Reaction definitions:')
  message(paste('Reaction 1:', rct1))
  message(paste('Reaction 2:', rct2))
}

readline("\nType  <Return>\t to start : ")
{
  message('Reaction similarity between rct1 and rct2 using \'msim\' algorithm: ')
  sim <- rs.compute(rct1, rct2, verbose = T)
  message(sim)
}

readline("\nType  <Return>\t to start : ")
{
  message('Reaction similarity between rct1 and rct2 using \'rsim\' algorithm: ')
  sim <- rs.compute(rct1, rct2, algo = 'rsim', verbose = T)
  message(sim)
}

readline("\nType  <Return>\t to start : ")
{
  message('Reaction similarity between rct1 and rct2 using \'rsim2\' algorithm: ')
  sim <- rs.compute(rct1, rct2, algo = 'rsim2')
  message(sim)
}