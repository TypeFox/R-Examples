#' NNTintervalsProspective
#'
#' Produce Bayesian and classical intervals for NNT from
#' observations in a prospective study.
#' Useful for "anticipated results" when designing a study,
#' The setting: patients will be tested immediately, and followed to determine the BestToTreat/BestToWait classification.
#' as well as analyzing study results.
#' There were (or will be) Npositives patients with a positive test, Nnegatives with a negative test.
#' The observed NNTs in each group were (or will be) NNTpos and NNTneg.
#'
#' @param Npositives Total number of observed positives.
#' @param Nnegatives Total number of observed negatives.
#' @param NtruePositives Observed or anticipated number of "BestToTreat" among the positives.
#' @param NtrueNegatives Observed or anticipated number of "BestToWait"  among the negatives.
#' @param prev = 0.15 Prevalence of "BestToTreat" characteristic.
#' @param alpha = 0.025 Significance level (one side).
#' @param prior Beta parameters for prior. Default is the Jeffreys prior = c(1/2,1/2). Jaynes prior = c(0,0) won't work when #fp=1.
#' @return The Bayesian predictive intervals for NNTpos and NNTneg. These are obtained from predictive intervals
#' for PPV and NPV, based on Jeffreys' beta(1/2,1/2) prior.

NNTintervalsProspective = function(
  Npositives,
  Nnegatives,
  NtruePositives,
  NtrueNegatives,
  prev = 0.15,
  alpha = 0.025,
  prior = c(1/2, 1/2)
){


  tp = NtruePositives
  fn = Nnegatives - NtrueNegatives
  tn = NtrueNegatives
  fp = Npositives - NtruePositives
  ppv = tp/(tp+fp)
  npv = tn/(tn+fn)
  NNTpos = 1/ppv
  NNTneg = 1/(1-npv)

  # predictive intervals
  ppvPI = qbeta(c(alpha, 1-alpha),
                tp+prior[1]-1, fp+prior[2]-1)
  #ifVerboseCat("  #CI for SN", NNTbiomarker::binom.confint(tp, Npositives))
  NNTposPI <- 1/ppvPI[2:1]

  npvPI = qbeta(c(alpha, 1-alpha),
                tn+prior[1]-1, fn+prior[2]-1)
  #ifVerboseCat("  #CI for SP", NNTbiomarker::binom.confint(tn, Nnegatives))
  NNTnegPI <- 1/(1-npvPI)
  result = cbind(NNTposPI, NNTnegPI, ppvPI, npvPI)
  result = rbind(c(NNTpos, NNTneg, ppv, npv), result)
  dimnames(result)[[2]] = c("NNTpos", "NNTneg", "PPV", "NPV")
  dimnames(result)[[1]] = c("estimates", "lower boundary", "upper boundary")
  result = result[ c(2,1,3), ]
  return(result)
}



#' NNTintervalsRetrospective
#'
#'Bayes predictive intervals for sensitivity, specificity, NNTpos and NNTneg
#'in a case-control retrospective study.
#'
#' @param Ncases Number of cases in the study
#' @param Ncontrols Number of controls in the study
#' @param NposCases Number of cases with positive test
#' @param NposControls Number of controls with positive test
#' @param prev Prevalence of the BestToTreat (versus BestToWait)
#' @param alpha Significance level for interval.
#' @param prior Beta parameters for prior. Default is the Jeffreys prior = c(1/2,1/2). Jaynes prior = c(0,0) won't work when #fp=1.
#' @return A list with 3 components containing intervals (predictive or otherwise), with names intervalsForSN, intervalsForSP, intervalsForNNT.
#' The intervals derive from assuming independent Jeffreys priors for SN and SP,
#'sampling from joint independent posteriors for SN and SP incorporating the
#'anticipated results, and applying NNT.from.sesp (Bayes theorem) to each
#'sampled pair to obtain a sample of NNTpos and NNTneg.

NNTintervalsRetrospective = function(
  Ncases = 10,
  Ncontrols = 30,
  NposCases = 6,
  NposControls = 2,
  prev = 0.15,
  alpha = 0.025,
  prior = c(1/2,1/2)
){
  # For a given sensitivity,
  spFunctionPos = function(sn, nntvalue, prev.=prev)
    c(sp=1 - sn*(nntvalue-1)*prev./(1-prev.))
  spFunctionNeg = function(sn, nntvalue, prev.=prev)
    c(sp=(1-sn)*(nntvalue-1)*prev./(1-prev.))

  # Validations.
  #   NNT.from.sesp(se=.3, sp=spFunctionPos(.3, 2), prev=prev)[1]
  #   NNT.from.sesp(se=.3, sp=spFunctionNeg(.3, 2), prev=prev)[2]
  #   # this confirms that spFunctionPos and Neg are OK.
  #ptemp <- (1:99)/100
  #plot(spFunctionPos(ptemp, 2), spFunctionNeg(ptemp, 2))
  #   print(  pv.from.sesp(se=.8, sp=spFunctionPos(.8, 2), prev=prev))
  #   print(  pv.from.NNT(NNTpos = NNTpos, NNTneg = NNTneg, prev=prev))
  ### note- specificity of "no test" is already 1 - prevalence!
  # print(sesp.from.NNT(NNTpos = NNTpos, NNTneg=NNTneg, prev=prev))
  # The same answer. Good.

  # Estimated sensitivity.
  snHat = NposCases/Ncases
  tp = NposCases
  fn = Ncases - NposCases

  # Estimated specificity.
  spHat = 1 - NposControls/Ncontrols
  tn = Ncontrols - NposControls
  fp = NposControls

  NNThat = NNT.from.sesp(se = snHat, sp=spHat,prev=prev)
  #
  # Same prior for sensitivity and specificity.
  Asn1=tp+prior[1]; Asn2=fn+prior[2];
  Asp1=fp+prior[2]; Asp2=tn+prior[2];

  #   #  Probability density
  #   kernelPos = function(sn, nntvalue.. ){
  #     sp = spFunctionPos(sn, nntvalue=nntvalue..)
  #     dbeta(sn, shape1=Asn1-1, shape2=Asn2-1) *
  #       pbeta(sp, shape1=Asp1-1, shape2=Asp2-1) *
  #       (prev/(1-prev))^2 *sn^2 *(1-sn)^2 *
  #       beta(tn,fp-1)/beta(tn,fp)
  #   }
  #
  #   kernelNegInner = function(sp)
  #     dbeta(sp, shape1=Asp1-1,
  #           shape2=Asp2-1) /
  #     abs(2*sp-1)
  #   #    beta(tn,fp-1)/beta(tn,fp)
  #
  #   kernelNeg = function(sn, nntvalue..){
  #     #spUpper = spFunctionNeg(sn, nntvalue=nntvalue..)
  #     upperLimit = spFunctionNeg(sn, nntvalue=nntvalue..)
  #     multiplier = (prev/(1-prev) * sn * (1-sn)) ^2 *
  #       dbeta(sn, shape1=Asn1-1, shape2=Asn2-1)
  #     integrals = sapply(sn, function(sn1)
  #       integrate(
  #         lower=0,
  #         upper=upperLimit,
  #         kernelNegInner
  #       )$value
  #     )
  #     return(multiplier * integrals)
  #   }
  #
  #   marginalNNTpos = function(nntvalue.)
  #     integrate(f=kernelPos,
  #               nntvalue..=nntvalue.,
  #               lower = 0.001, upper = 0.999
  #     )$value
  #
  #
  #   marginalNNTneg = function(nntvalue.)
  #     integrate(f=kernelNeg,
  #               nntvalue..=nntvalue.,
  #               lower = 0.001, upper = 0.999
  #     )$value
  #
  #
  #   plot(sn<-seq(.01,.99,.01), sapply(sn, kernelPos, nnt=2))
  #   lines(sn<-seq(.01,.99,.01), sapply(sn, kernelPos, nnt=3))
  #
  #   plot(nnt<-seq(1.1,10,.1),
  #        1-(sapply(nnt, marginalNNTpos ))
  #   )
  #   marginalNNTpos(1)   # 0.04280385
  #   plot(nnt<-seq(1,50,.1),
  #        sapply(nnt, marginalNNTneg ))
  #   marginalNNTneg(1000)  # 0.001662432
  #   marginalNNTneg(1)  # 0

  snRandom<-rbeta(100000, shape1=Asn1, shape2=Asn2)
  intervalForSN_sim = quantile(snRandom, probs = c(0.025, 0.975))
  intervalForSN_q = qbeta( p = c(0.025, 0.975), Asn1, Asn2)
  ciForSN = binom.confint(NposCases, Ncases)
  intervalsForSN = cbind(t(data.frame(intervalForSN_sim, intervalForSN_q, ciForSN)), snHat)

  spRandom<-1-rbeta(100000, shape1=Asp1, shape2=Asp2)
  intervalForSP_sim = quantile(spRandom, probs = c(0.025, 0.975))
  intervalForSP_q = qbeta( p = c(0.025, 0.975), Asp2, Asp1)
  ciForSP = binom.confint(Ncontrols - NposControls, Ncontrols)
  intervalsForSP = cbind(t(data.frame(intervalForSP_sim, intervalForSP_q, ciForSP)), spHat)

  # plot(snRandom, spRandom)
  NNTrandom = NNT.from.sesp(se=snRandom, sp=spRandom, prev=prev)
  NNTposRandom = NNTrandom[[1]]
  NNTnegRandom = NNTrandom[[2]]
  summary(NNTposRandom)
  summary(NNTnegRandom)
  summary(NNTnegRandom - NNTposRandom)  # A tiny number have NNTneg < NNTpos
  ###  All this looks pretty good.
  intervalForNNTpos = quantile(NNTposRandom, probs = c(alpha, 1-alpha))
  intervalForNNTneg = quantile(NNTnegRandom, probs = c(alpha, 1-alpha))
  ###  Reasonable.
  # cor(NNTpos, NNTneg) ## Negative
  intervalsForNNT = cbind(t(data.frame(intervalForNNTpos,
                                       intervalForNNTneg)), NNThat)
  result = cbind(
    intervalForNNTpos,
    intervalForNNTneg,
    intervalForSN_q,
    intervalForSP_q
  )
  result = rbind(c(NNThat, snHat, spHat), result)
  dimnames(result)[[2]] = c("NNTpos", "NNTneg", "sensitivity", "specificity")
  dimnames(result)[[1]] = c("estimates", "lower boundary", "upper boundary")
  result = result[c(2,1,3), ]
  return(result)
  # plot(NNTpos, NNTneg, log="y", pch=".", cex=0.1 , col=c("black", "green")[1+(NNTpos < NNTneg)])
  # plot(spRandom, snRandom, col=c("black", "green")[1+(NNTpos < NNTneg)])
  #   par(mfrow=c(1,2))
  #   plot(density(NNTpos), xlab=NA,main=NA)
  #   plot(density(NNTneg), xlab=NA,main=NA, xlim=c(1,70))
  #
  #   findCDFforNNTneg = function(plimit)
  #     uniroot(f = function(nnt) plimit -
  #               marginalNNTneg(nnt),
  #             interval=c(1.1,10000)
  #     )
  #   lbNNTneg = findCDFforNNTneg(alpha)$root    ## should be near 16
  #   ubNNTneg = findCDFforNNTneg(1-alpha)$root  ## should be near 88
  #   return(lbNNTneg = lbNNTneg, ubNNTneg = ubNNTneg)
}
