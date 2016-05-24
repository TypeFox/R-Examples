`simulateRegression.fnc` <-
function(
  beta = c(400, 2, 6, 4), 
  nitem = 20, nsubj = 10, 
  stdevItem = 40, stdevSubj = 80, stdevError = 50,
  nruns = 100, 
  learn = FALSE, learnRate = 10, ...) 
{

  stop("this function is not working due to changes in lme4.
        an update using lmerTest is in progress")

  require("lme4", quietly = TRUE)

  # some matrices for storing results
  pLmer = matrix(0, nruns, 4)         # p's for lmer correct
  pLmerS = matrix(0, nruns, 4)        # p's for lmer Subject
  #pLmerN = matrix(0, nruns, 4)       # p's for lmer nested Items
  pLmList = matrix(0, nruns, 4)       # p's for lmList
  pItemList = matrix(0, nruns, 4)     # p's for by-item regression

  MCMCpLmer = matrix(0, nruns, 4)     # MCMC p-values
  MCMCpLmerS = matrix(0, nruns, 4)    # MCMC p-values

  cLmer = matrix(0, nruns, 7)         # coef's for lmer correct
  cLmerS = matrix(0, nruns, 6)        # coef's for lmer Subject
  #cLmerN = matrix(0, nruns, 7)       # coef's for lmer nested Items
  cLmList = matrix(0, nruns, 4)       # coef's for lmList
  cItemList = matrix(0, nruns, 4)     # coef's for by-item regression

  for (run in 1:nruns) {

    data2 = make.reg.fnc(
      nsubj = nsubj, nitem = nitem,
      stdevItem = stdevItem, stdevSubj = stdevSubj, 
      stdevError = stdevError, 
      beta = beta, learn = learn, learnRate = learnRate) 

    # data2$Item2 = 
    #   factor(paste(data2$Subject, data2$Item, sep="."))
    # correct lmer analysis
    if (learn) {
      data2.lmer = lmer(RT~X+Y+Z+Trial+(1|Subject)+(1|Item), data=data2)
    } else {
      data2.lmer = lmer(RT~X+Y+Z+(1|Subject)+(1|Item), data=data2)
    }
    coef.ranefs = c(VarCorr(data2.lmer)[[1]]@factors$correlation@sd,
                    VarCorr(data2.lmer)[[2]]@factors$correlation@sd,
                    attr(summary(data2.lmer),"sigma"))
		names(coef.ranefs) = c("Item", "Subject", "Residual")
    x = pvals.fnc(data2.lmer, withMCMC = FALSE)  # does not return mcmc object
    coef.fixefs = x$fixed$Estimate[1:4]
    names(coef.fixefs) = rownames(x$fixed)
    cLmer[run,] = c(coef.ranefs, coef.fixefs)
    pLmer[run,] = x$fixed$"Pr(>|t|)"[1:4]
    MCMCpLmer[run,] = as.numeric(x$fixed$pMCMC)[1:4]
    if (run == 1) {
      colnames(pLmer) = names(coef.fixefs)
      colnames(MCMCpLmer) = names(coef.fixefs)
    }

    # only subject
    if (learn) {
      data2.lmerS = lmer(RT~X+Y+Z+Trial+(1|Subject), data=data2)
    } else {
      data2.lmerS = lmer(RT~X+Y+Z+(1|Subject), data=data2)
    }
    coef.ranefs = c(VarCorr(data2.lmerS)[[1]]@factors$correlation@sd,
                    attr(summary(data2.lmerS),"sigma"))
		names(coef.ranefs) = c("Subject", "Residual")
    x = pvals.fnc(data2.lmerS, withMCMC = FALSE)  # does not return mcmc object
    coef.fixefs = x$fixed$Estimate[1:4]
    names(coef.fixefs) = rownames(x$fixed)
    cLmerS[run,] = c(coef.ranefs, coef.fixefs)
    pLmerS[run, ] = x$fixed$"Pr(>|t|)"[1:4]
    MCMCpLmerS[run, ] = as.numeric(x$fixed$pMCMC)[1:4]
    if (run == 1) {
      colnames(pLmerS) = names(coef.fixefs)
      colnames(MCMCpLmerS) = names(coef.fixefs)
    }

    ## items nested
    #data2.lmerN = lmer(RT~X+Y+Z+(1|Subject)+(1|Item2), data=data2)
    #x = pvals.fnc(data2.lmerN,nsamp=0) 
    #coef.ranefs = as.numeric(as.character(x$random[,"Std.Dev"]))
    #coef.fixefs = x$summary$Estimate
    #cLmerN[run,] = c(coef.ranefs, coef.fixefs)
    #pLmerN[run, ] = x$summary$pvals
    #if (run == 1) {
    #  colnames(cLmerN) = c(as.character(x$random$Groups), 
    #                       rownames(x$summary))
    #  colnames(pLmerN) = rownames(x$summary)
    #}

    if (learn) {
      data2.lmList = lmList(RT~X+Y+Z+Trial|Subject, data=data2)
    } else {
      data2.lmList = lmList(RT~X+Y+Z|Subject,data=data2)
    }

    x = apply(coef(data2.lmList),2,t.test) 
    pLmList[run,] = c(x[[1]]$p.value, x[[2]]$p.value, 
      x[[3]]$p.value, x[[4]]$p.value)
    cLmList[run,]=as.numeric(apply(coef(data2.lmList),2,mean)[1:4])

    x = summary(item.fnc(data2))$coef
    pItemList[run,]=x[,4]
    cItemList[run,]=x[,1]
    cat(run, " ")
  }
  cat("\n")

  a05lmer = apply(pLmer < 0.05, 2, sum)/nruns
  a01lmer = apply(pLmer < 0.01, 2, sum)/nruns
  a05lmerS = apply(pLmerS < 0.05, 2, sum)/nruns
  a01lmerS = apply(pLmerS < 0.01, 2, sum)/nruns
  #a05lmerN = apply(pLmerN < 0.05, 2, sum)/nruns
  #a01lmerN = apply(pLmerN < 0.01, 2, sum)/nruns
  a05lmlist = apply(pLmList < 0.05, 2, sum)/nruns
  a01lmlist = apply(pLmList < 0.01, 2, sum)/nruns
  a05Itemlist = apply(pItemList < 0.05, 2, sum)/nruns
  a01Itemlist = apply(pItemList < 0.01, 2, sum)/nruns

  mcmclmer05  = apply(MCMCpLmer < 0.05, 2, sum)/nruns
  mcmclmer01  = apply(MCMCpLmer < 0.01, 2, sum)/nruns
  mcmclmerS05 = apply(MCMCpLmerS < 0.05, 2, sum)/nruns
  mcmclmerS01 = apply(MCMCpLmerS < 0.01, 2, sum)/nruns

  #alpha05 = rbind(a05lmer, a05lmerS, a05lmerN, a05lmlist, a05Itemlist)
  #alpha01 = rbind(a01lmer, a01lmerS, a01lmerN, a01lmlist, a01Itemlist)
  alpha05 = rbind(a05lmer, mcmclmer05, a05lmerS, mcmclmerS05, a05lmlist, a05Itemlist)
  alpha01 = rbind(a01lmer, mcmclmer01, a01lmerS, mcmclmerS01, a01lmlist, a01Itemlist)
  colnames(alpha05) = c("Intercept", "X", "Y", "Z")
  colnames(alpha01) = c("Intercept", "X", "Y", "Z")
  #rownames(alpha05) = c("lmer", "lmerS", "lmerN", "lmList", "item")
  #rownames(alpha01) = c("lmer", "lmerS", "lmerN", "lmList", "item")
  rownames(alpha05) = c("lmer", "lmer-mcmc", "lmerS", "lmerS-mcmc", "lmList", "item")
  rownames(alpha01) = c("lmer", "lmer-mcmc", "lmerS", "lmerS-mcmc", "lmList", "item")

  
	rn = rbind(
	 c(mean(as.numeric(cLmer[,1])),
	   mean(as.numeric(cLmer[,2])),
	   mean(as.numeric(cLmer[,3]))),
	 c(NA, mean(as.numeric(cLmerS[,1])),
     mean(as.numeric(cLmerS[,2]))))
  colnames(rn) = c("Item", "Subject", "Residual")
  rownames(rn) = c("lmer", "lmerS")
   
  return(list(alpha05 = alpha05, alpha01 = alpha01, ranef = rn))
}

