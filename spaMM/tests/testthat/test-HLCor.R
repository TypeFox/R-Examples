cat("\ntest HLCor:")
# HLCor

data(Loaloa)
hl <- HLCor(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
      +Matern(1|longitude+latitude),data=Loaloa,
      family=binomial(),ranPars=list(nu=0.5,rho=1/0.7)) ## takes ~ 6s

expect_equal(hl$APHLs$p_v,-645.7328,tolerance=1e-04)

data(scotlip)
hl <- HLCor(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(scotlip$expec)),
      ranPars=list(rho=0.174),
      adjMatrix=Nmatrix,family=poisson(),data=scotlip)

expect_equal(hl$APHLs$p_v,-161.5143,tolerance=1e-04)
