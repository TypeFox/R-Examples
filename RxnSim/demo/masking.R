{
  rct1 <- 'CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(OP(=O)(OCC1(C(OP([O-])(=O)[O-])C(O)C(O1)N3(C2(=C(C(N)=NC=N2)N=C3))))[O-])[O-].O>>CC([O-])=O.CC(C)(C(O)C(=O)NCCC(=O)NCCS)COP(=O)(OP(=O)(OCC1(OC(C(C1OP([O-])(=O)[O-])O)N3(C2(=C(C(N)=NC=N2)N=C3))))[O-])[O-]'
  rct2 <- 'CC(C)(C(O)C(=O)NCCC(=O)NCCSC(C1(=CC=C(O)C=C1))=O)COP(=O)(OP(=O)(OCC2(C(OP([O-])(=O)[O-])C(O)C(O2)N4(C3(=C(C(N)=NC=N3)N=C4))))[O-])[O-].O>>C(C1(C=CC(=CC=1)O))(=O)[O-].CC(C)(C(O)C(=O)NCCC(=O)NCCS)COP(=O)(OP(=O)(OCC1(OC(C(C1OP([O-])(=O)[O-])O)N3(C2(=C(C(N)=NC=N2)N=C3))))[O-])[O-]'
  
  message('Reaction definitions:')
  message(paste('Reaction 1:', rct1))
  message(paste('Reaction 2:', rct2))
}

readline("\nType  <Return>\t to start : ")
{
  message('Reaction similarity between rct1 and rct2: ')
  sim <- rs.compute(rct1, rct2, verbose = TRUE)
  message(sim)
}

readline("\nType  <Return>\t to start : ")
{
  r1 <- rs.mask('CC(C)(C(O)C(=O)NCCC(=O)NCCS)COP(=O)(OP(=O)(OCC1(OC(C(C1OP([O-])(=O)[O-])O)n3(c2(c(c(N)ncn2)nc3))))[O-])[O-]', '[Cs]', rct1)
  r2 <- rs.mask('CC(C)(C(O)C(=O)NCCC(=O)NCCS)COP(=O)(OP(=O)(OCC1(OC(C(C1OP([O-])(=O)[O-])O)n3(c2(c(c(N)ncn2)nc3))))[O-])[O-]', '[Cs]', rct2)
  
  message('Reaction SMILES post masking CoA:')
  message(paste('Reaction 1:', r1))
  message(paste('Reaction 2:', r2))
}

readline("\nType  <Return>\t to start : ")
{
  message('Reaction similarity between rct1 and rct2 (with CoA masked): ')
  sim <- rs.compute(r1, r2, verbose = TRUE)
  message(sim)
}

readline("\nType  <Return>\t to start : ")
{
  mol1 <- 'CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(OP(=O)(OCC1(C(OP([O-])(=O)[O-])C(O)C(O1)N3(C2(=C(C(N)=NC=N2)N=C3))))[O-])[O-]'
  mol2 <- 'CC(C)(C(O)C(=O)NCCC(=O)NCCSC(C1(=CC=C(O)C=C1))=O)COP(=O)(OP(=O)(OCC2(C(OP([O-])(=O)[O-])C(O)C(O2)N4(C3(=C(C(N)=NC=N3)N=C4))))[O-])[O-]'
  
  message('Molecule SMILES:')
  message(paste('Molecule 1:', mol1))
  message(paste('Molecule 2:', mol2))
}

readline("\nType  <Return>\t to start : ")
{
  message('Similarity between mol1 and mol2: ')
  sim <- ms.compute(mol1, mol2)
  message(sim)
}

readline("\nType  <Return>\t to start : ")
{
  m1 <- ms.mask('CC(C)(C(O)C(=O)NCCC(=O)NCCS)COP(=O)(OP(=O)(OCC1(OC(C(C1OP([O-])(=O)[O-])O)n3(c2(c(c(N)ncn2)nc3))))[O-])[O-]', '[Cs]', mol1)
  m2 <- ms.mask('CC(C)(C(O)C(=O)NCCC(=O)NCCS)COP(=O)(OP(=O)(OCC1(OC(C(C1OP([O-])(=O)[O-])O)n3(c2(c(c(N)ncn2)nc3))))[O-])[O-]', '[Cs]', mol2)
  
  message('Molecule SMILES post masking CoA:')
  message(paste('Molecule 1:', m1))
  message(paste('Molecule 2:', m2))
}

readline("\nType  <Return>\t to start : ")
{
  message('Similarity between mol1 and mol2 (with CoA masked): ')
  sim <- ms.compute(m1, m2)
  message(sim)
}