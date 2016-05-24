gof = function(x, priority = "tests", copula = NULL, tests = NULL, M = 50, MJ = 50, param = 0.5, param.est = T, df = 4, df.est = T, m = 1, delta.J = 0.5, nodes.Integration = 12, m_b = 0.5, zeta.m = 0, b_Rn = 0.05){
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if(any(lapply(tests, function(x){is.element(x,c("gofPIOSRn", "gofPIOSTn", "gofKernel", "gofRosenblattSnB", 
                                                                 "gofRosenblattSnC", "gofADChisq", "gofADGamma", 
                                                                 "gofSn", "gofKendallCvM", "gofKendallKS",
                                                                 "gofWhite", "gofRn"))}) == F)==T){stop("At least one of the tests in 'testset' is not implemented")}
  if(any(lapply(copula, function(x){is.element(x,c(# "independence", 
                                                   "gaussian", "t", "clayton", "frank", "gumbel"
                                                   #, "joe", "amh", "survival clayton", "survival gumbel", "survival joe", "90 clayton", "90 gumbel", "90 joe", "270 clayton", "270 gumbel", "270 joe", "bb1", "bb6", "bb7", "bb8", "survival bb1", "survival bb6", "survival bb7", "survival bb8", "90 bb1", "90 bb6", "90 bb7", "90 bb8", "270 bb1", "270 bb6", "270 bb7", "270 bb8"
  ))}) == F)==T){stop("At least one of the copulae is not implemented")}
  if (is.element(priority, c("tests", "copula")) == F){
    stop("Please insert a valid character string for the argument priority. It shall be either tests or copula.")
  }
  
  if (any(x > 1) || any(x < 0)){
      warning("The observations aren't in [0,1]. The margins will be estimated by the ranks of the observations.", sep = "")
    x = .margins(x, margins = "ranks")
  }
  
  if (!is.null(tests)){
    if (!is.null(copula)){
      if (M >= 100 | MJ >= 100){
      times.comp = lapply(copula, gofCheckTime, x=x, test = tests, M = M, MJ = MJ, print.res = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m)
      print(.get.time(Reduce("+", times.comp)))
      }
      lapply(copula, gofHybrid, x=x, testset = tests, M = M, MJ = MJ, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m)
    } else {
      copula_list = lapply(tests, gofWhichCopula)
      copula = Reduce(intersect, copula_list)
      if (M >= 100 | MJ >= 100){
      times.comp = lapply(copula, gofCheckTime, x=x, test = tests, M = M, MJ = MJ, print.res = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m)
      print(.get.time(Reduce("+", times.comp)))
      }
      lapply(copula, gofHybrid, x=x, testset = tests, M = M, MJ = MJ, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m)
    }
  } else {
    if (!is.null(copula)){
      tests_list = lapply(copula, gofWhich, d = dim(x)[2])
      tests = Reduce(intersect, tests_list)
      tests = tests[-which(tests == "gofHybrid")]
#      copula_list = lapply(tests, gofWhichCopula)
#      copula = Reduce(intersect, copula_list)
      if (M >= 100 | MJ >= 100){
      times.comp = lapply(copula, gofCheckTime, x=x, test = tests, M = M, MJ = MJ, print.res = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m)
      print(.get.time(Reduce("+", times.comp)))
      }
      lapply(copula, gofHybrid, x=x, testset = tests, M = M, MJ = MJ, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m)
    } else if (priority == "tests") {
      tests = gofWhich("gaussian", dim(x)[2])
      tests = tests[-which(tests == "gofHybrid")]
      copula_list = lapply(tests, gofWhichCopula)
      copula = Reduce(intersect, copula_list)
      if (M >= 100 | MJ >= 100){
      times.comp = lapply(copula, gofCheckTime, x=x, test = tests, M = M, MJ = MJ, print.res = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m)
      print(.get.time(Reduce("+", times.comp)))
      }
      lapply(copula, gofHybrid, x=x, testset = tests, M = M, MJ = MJ, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m)
    } else if (priority == "copula"){
      copula = c("gaussian", "t", "frank", "gumbel", "clayton")
      tests_list = lapply(copula, gofWhich, d = dim(x)[2])
      tests = Reduce(intersect, tests_list)
      tests = tests[-which(tests == "gofHybrid")]
      if (M >= 100 | MJ >= 100){
      times.comp = lapply(copula, gofCheckTime, x=x, test = tests, M = M, MJ = MJ, print.res = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m)
      print(.get.time(Reduce("+", times.comp)))
      }
      lapply(copula, gofHybrid, x=x, testset = tests, M = M, MJ = MJ, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m)
    }
  }
}