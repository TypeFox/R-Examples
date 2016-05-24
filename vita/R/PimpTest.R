PimpTest = function(Pimp, ...)  UseMethod("PimpTest")

PimpTest.default = function(Pimp, para = FALSE, ...){
  ##################################################################
  # randomForest?
  if (!inherits(Pimp, "PIMP")) stop("Pimp is not of class PIMP")
  p = nrow(Pimp$PerVarImp)
  if(para){

        # mean and sample variance
        mean.PerVarImp = apply(Pimp$PerVarImp, 1, mean)
        sd.PerVarImp = apply(Pimp$PerVarImp, 1, sd)

        # Kolmogorov-Smirnov Test
        i=1:p
        test.norm = sapply(i, function(i){
          ks.test(unique(Pimp$PerVarImp[i,]), "pnorm",mean = mean.PerVarImp[i],sd = sd.PerVarImp[i])$p.value
        })

        # the p-value is the probability of observing the original VarImp or one more extreme,
        p.val = sapply(i, function(i){
          pnorm(Pimp$VarImp[i], mean.PerVarImp[i], sd.PerVarImp[i], lower.tail=F)
        })

  } else{
        i=1:p
        # The empirical cumulative null distribution functions Fn0^{*s} of VarImp^{*s}
        Fn0 = sapply(i, function(i){
          ecdf(Pimp$PerVarImp[i,])
        })

        # the p-value is the probability of observing the original VarImp or one more extreme,
        p.val = sapply(i, function(i){
          1-Fn0[[i]](Pimp$VarImp[i])
        })
  }
  #
  dimNames=dimnames(Pimp$VarImp)[[1]]
  out=list(VarImp = Pimp$VarImp,
           PerVarImp = Pimp$PerVarImp,
           para = para ,
           meanPerVarImp =  if(para){ matrix(mean.PerVarImp, ncol=1, dimnames=list(dimNames, "mean(PerVarImp)")) } else { NULL },
           sdPerVarImp =if(para){ matrix(sd.PerVarImp, ncol=1, dimnames=list(dimNames, "sd(PerVarImp)")) } else { NULL },
           p.ks.test =  if(para){ matrix(test.norm,ncol=1, dimnames=list(dimNames, "ks.test"))  } else { NULL},
           pvalue = matrix(p.val, ncol=1, dimnames=list(dimNames, "p-value")),
           type=Pimp$type,
           call.PIMP=Pimp$call,
           call = match.call())


  class(out) = "PimpTest"
  return(out)

}

print.PimpTest = function(x, ...){
  if (!inherits(x, "PimpTest")) stop(" is not of class PIMP.Test")
  cat("Call:\n \n")
  print(x$call)
  cat("type: ")
  print(x$type)
  cat("\noriginal VarImp:\n")
  print(x$VarImp)
  cat("\np-value:\n")
  print(x$pvalue)
}
