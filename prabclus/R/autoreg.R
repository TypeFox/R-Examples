"autoreg" <-
function(x, probs, ejprob, plot=TRUE, nperp=4, species.fixed=TRUE, pdfnb=FALSE){
  tjumps <- matrix(nrow=nperp*length(probs),ncol=3)
  for (j in 1:length(probs)){
    cat("    Estimating disj. parameter: Simulations for p= ",probs[j],"\n")
    for (i in 1:nperp){
#      print(x$n.species)
#      print(x$specperreg)
#      print(sum(x$specperreg))
      test <- randpop.nb(x$nb,p.nb=probs[j], n.species=x$n.species,
                         vector.species=x$regperspec,
                         species.fixed=species.fixed,
                         pdf.regions=x$specperreg/sum(x$specperreg),
                         count=FALSE, pdfnb=pdfnb)
#      print("Test generated")
#      print(test)
#      print(dim(test))
#      print(x$regperspec)
      nst <- apply(test,2,sum)
#      print(nst)
#      print(apply(test,1,sum) - x$specperreg)
      tcn <- con.regmat(test,x$nb,count=FALSE)
#      print("tcn computed")
      ind <- (j-1)*nperp+i
      tjumps[ind,1] <- probs[j]
      tjumps[ind,2] <- sum(tcn-1)
#      print(tjumps[ind,2])
      tjumps[ind,3] <- tjumps[ind,2]/sum(nst-1)
    }
  }
  ejlm <- lm(tjumps[,3]~tjumps[,1])
  if (plot){
    plot(tjumps[,1],tjumps[,3],xlab="pdisj",
         ylab="qdisj")
#    print(tjumps[,1])
#    print(tjumps[,3])
    abline(ejlm$coef)
    abline(c(ejprob,0), lty=2)
  }
  pd <- (ejprob-ejlm$coef[1])/ejlm$coef[2]
  cat("  Estimated disjunction parameter =",pd,"\n")
  out <- list(pd=pd, coef=ejlm$coef)
  out
}
