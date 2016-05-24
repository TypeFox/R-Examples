gofHybrid = function(copula, x, testset = c("gofPIOSRn", "gofKernel"), margins = "ranks", M = 1000, execute.times.comp = T, param = 0.5, param.est = T, df = 4, df.est = T, m = 1, MJ = 100, delta.J = 0.5, nodes.Integration = 12, m_b = 0.5, zeta.m = 0, b_Rn = 0.05){
  if(any(apply(as.matrix(testset), 2, function(x){is.element(x,c("gofPIOSRn", "gofPIOSTn", "gofKernel", "gofRosenblattSnB", 
                                                                 "gofRosenblattSnC", "gofADChisq", "gofADGamma", 
                                                                 "gofSn", "gofKendallCvM", "gofKendallKS",
                                                                 "gofWhite", "gofRn"))}) == F)==T){stop("At least one of the tests in 'testset' is not implemented")}
  if (dim(x)[2] > 2 & any(apply(as.matrix(testset), 2, function(x){is.element(x,c("gofPIOSRn", "gofPIOSTn", "gofKernel", "gofKendallCvM", "gofKendallKS", "gofWhite"))}) == T)){
    stop("At least one test in the testset can't handle dimensions greater 2. Please use the tests gofRosenblattSnB, 
          gofRosenblattSnC, gofADChisq, gofADGamma, gofSn or gofRn.")
  }
  
  if (any(x > 1) || any(x < 0)){
    if (margins == "ranks"){
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " of the observations.", sep = ""))
    } else {
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " distribution.", sep = ""))
    }
    
    x = .margins(x, margins)
  }
  
  if (execute.times.comp == T & M >= 100){
  times.comp = gofCheckTime(copula = copula, x=x, test = testset, M = M, print.res = F, param = param, param.est = param.est, df = df, df.est = df.est, m = m, MJ = MJ, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn)
  print(.get.time(times.comp))
  }
  
  pres1 = c()
  pres = c()
  tres1 = c()
  tres = c()
  if (is.element("gofPIOSRn", testset) == T){
    cop = try(gofPIOSRn(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est))
    if (class(cop) == "try-error"){pres = tres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
    }
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
  }
  if (is.element("gofPIOSTn", testset) == T){
    cop = try(gofPIOSTn(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, m = m))
    if (class(cop) == "try-error"){pres = tres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
    }
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
  }
  if (is.element("gofKernel", testset) == T){
    cop = try(gofKernel(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, MJ = MJ, delta.J = delta.J, nodes.Integration = nodes.Integration))
    if (class(cop) == "try-error"){pres = tres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
    }
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
  }
  if (is.element("gofRn", testset) == T){
    cop = try(gofRn(copula = copula, x = x, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn))
    if (class(cop) == "try-error"){pres = tres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
    }
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
  }
  if (is.element("gofKendallCvM", testset) == T){
    cop = try(gofKendallCvM(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est))
    if (class(cop) == "try-error"){pres = tres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
    }
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
  }
  if (is.element("gofKendallKS", testset) == T){
    cop = try(gofKendallKS(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est))
    if (class(cop) == "try-error"){pres = tres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
    }
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
  }
  if (is.element("gofWhite", testset) == T){
    cop = try(gofWhite(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est))
    if (class(cop) == "try-error"){pres = tres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
    }
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
  }
  if (is.element("gofRosenblattSnB", testset) == T){
    cop = try(gofRosenblattSnB(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est))
    if (class(cop) == "try-error"){pres = tres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
    }
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
  }
  if (is.element("gofRosenblattSnC", testset) == T){
    cop = try(gofRosenblattSnC(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est))
    if (class(cop) == "try-error"){pres = tres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
    }
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
  }
  if (is.element("gofADChisq", testset) == T){
    cop = try(gofADChisq(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est))
    if (class(cop) == "try-error"){pres = tres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
    }
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
  }
  if (is.element("gofADGamma", testset) == T){
    cop = try(gofADGamma(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est))
    if (class(cop) == "try-error"){pres = tres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
    }
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
  }
  if (is.element("gofSn", testset) == T){
    cop = try(gofSn(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est))
    if (class(cop) == "try-error"){pres = tres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
    }
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
  }
  
  if (length(testset) > 1) {
    which_comb =list()
    for(i in 1:(2^length(testset))){
      which_comb[[i]]=which(as.integer(intToBits(i)) == 1)
    }
    comb_exist = which_comb[which(unlist(lapply(which_comb, length)) > 1)]
    
    for (i in 1:length(comb_exist)){
      pres1 = c(pres1, min(length(pres1[comb_exist[[i]]]) * min(pres1[comb_exist[[i]]]), 1))
      tres1 = c(tres1, NaN)
    }
    hybrid_comb_names = paste("hybrid(", lapply(comb_exist, paste, collapse = ", "), ")", sep="")
    matrix_names = c(substring(testset, 4), hybrid_comb_names)
    structure(class = "gofCOP",
              list(method = sprintf("Hybrid Goodness-of-fit test with the testset = c(%s) and the %s copula.", 
                                    toString(testset), copula), 
                   erg.tests = matrix(c(pres1, tres1), ncol = 2, 
                                      dimnames = list(matrix_names, c("p.value", "test statistic")))))
  } else {
  matrix_names = substring(testset, 4)
  structure(class = "gofCOP",
            list(method = sprintf("Goodness-of-fit test with the test = c(%s) and the %s copula.", 
                                  toString(testset), copula), 
                 erg.tests = matrix(c(pres1, tres1), ncol = 2, 
                                    dimnames = list(matrix_names, c("p.value", "test statistic")))))
  }
}
