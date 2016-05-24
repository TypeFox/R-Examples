cat("\ntest Nugget")

data(Loaloa)
## use 1st 30 obs as this is slow:
loafit <- corrHLfit(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
          +Matern(1|longitude+latitude),
          data=Loaloa[1:30,],family=binomial(),
          init.corrHLfit=list(Nugget=0.1),ranFix=list(nu=0.5)) 

expect_equal(attr(loafit$corrPars,"type"),list(nu="fix",Nugget="outer",rho="outer"))
expect_equal(loafit$corrPars$Nugget,0.04578246,tolerance=1e-5) ## NB flat p_bv for Nugget
