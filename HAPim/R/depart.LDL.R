`depart.LDL` <-
function(moyenne.pere,perf,CD,PLA,desc.pere){

  nbre.pere  =    length(desc.pere[,1])
  Y.all	     =	  CD*perf
  X.all      =    CD*(2*PLA-1)
  theta      =    rep(NA,nbre.pere)
  sigma      =    rep(NA,nbre.pere)


  for(i in 1:nbre.pere){

     deb      = desc.pere[i,1]
     fin      = desc.pere[i,2]
     n.ind    = fin-deb+1 
     X        = cbind(CD[deb:fin],X.all[deb:fin])
     Y        = Y.all[deb:fin]

        av      =	aov(Y~0+X)
        sigma[i]=sum(av$residuals*av$residuals)/n.ind
        theta[i]=av$coefficients[2]

   }


  alpha.Q	=	0

  if(sum(!is.na(theta))>0){

    alpha.Q  = max(abs(theta[!is.na(theta)]))

  }

     s         = sum(sigma)/nbre.pere

     depart    = c(s,alpha.Q)


     depart



}

