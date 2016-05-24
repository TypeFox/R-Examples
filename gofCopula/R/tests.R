gofPIOSRn = function(copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", execute.times.comp = T){
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (dim(x)[2] != 2){stop("x must be of dimension 2")}
  if (is.element(copula, c("gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofPIOSRn.")}
  dims = dim(x)[2]
  n = dim(x)[1]
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofPIOSRn", M = M, print.res = F)
    print(.get.time(times.comp))
  }
  
  if (copula == "gaussian"){copula = "normal"}
  if (df.est == T){df.fixed = F} else if (df.est == F){df.fixed = T}
  
  if (any(x > 1) || any(x < 0)){
    if (margins == "ranks"){
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " of the observations.", sep = ""))
    } else {
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " distribution.", sep = ""))
    }
    
    x = .margins(x, margins)
  }
  
  if ("normal" == copula || "t" == copula){
    if (param.est == T){
      param = try(fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed), data = x, method = "mpl")@estimate, silent = T)
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed), data = x, method = "itau")@estimate}
    }
    if (copula == "t" & df.fixed == F){
      copula = ellipCopula(copula, param = param[1], dim = dim(x)[2], df = param[2], df.fixed = T)
    } else {
      copula = ellipCopula(copula, param = param, dim = dim(x)[2], df = df, df.fixed = df.fixed)
    }
  } else if("clayton" == copula || "frank" == copula || "gumbel" == copula || 
              "amh" == copula || "joe" == copula){
    if (param.est == T){
      param = try(fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "mpl")@estimate, silent = T)
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "itau")@estimate}
    }
    copula = archmCopula(copula, param = param, dim = dim(x)[2])
  }
  
  bs.ac.c = c()
  ac.c  = .Rn(x = x, copula = copula, dims = dims)
  for(i in 1:M) {
    bs.ac.c[i] = .Rn(rCopula(n, copula), copula = copula, dims = dims)
  }
  test = mean(abs(bs.ac.c) > abs(ac.c))
 
structure(class = "gofCOP", 
          list(method = sprintf("Parametric bootstrap goodness-of-fit test (PIOS)"), 
               erg.tests = matrix(c(test, ac.c), ncol = 2, 
                                  dimnames = list("PIOSRn", c("p.value", "test statistic")))))
}

#################################################################################
gofPIOSTn = function(copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", m = 1, execute.times.comp = T){
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (dim(x)[2] != 2){stop("x must be of dimension 2")}
  if (is.element(copula, c("gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofPIOSTn.")}
  dims = dim(x)[2]
  n = dim(x)[1]
  B = n / m
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofPIOSTn", M = M, print.res = F)
    print(.get.time(times.comp))
  }
  
  if (copula == "gaussian"){copula = "normal"}
  if (df.est == T){df.fixed = F} else if (df.est == F){df.fixed = T}
  
  if (any(x > 1) || any(x < 0)){
    if (margins == "ranks"){
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " of the observations.", sep = ""))
    } else {
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " distribution.", sep = ""))
    }
    
    x = .margins(x, margins)
  }
  
  if ("normal" == copula || "t" == copula){
    if (param.est == T){
      param = try(fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed), data = x, method = "mpl")@estimate, silent = T)
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed), data = x, method = "itau")@estimate}
    }
    if (copula == "t" & df.fixed == F){
      copula = ellipCopula(copula, param = param[1], dim = dim(x)[2], df = param[2], df.fixed = T)
    } else {
      copula = ellipCopula(copula, param = param, dim = dim(x)[2], df = df, df.fixed = df.fixed)
    }
  } else if("clayton" == copula || "frank" == copula || "gumbel" == copula || 
              "amh" == copula || "joe" == copula){
    if (param.est == T){
      param = try(fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "mpl")@estimate, silent = T)
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "itau")@estimate}
    }
    copula = archmCopula(copula, param = param, dim = dim(x)[2])
  }
  
  bs.ac.c = c()
  ac.c  = .Tn(x = x, copula = copula, B = B, m = m, dims = dims, param.est = param.est)
  for(i in 1:M) {
    bs.ac.c[i] = .Tn(rCopula(n, copula), copula = copula, B = B, m = m, dims = dims, param.est = param.est)
  }
  test = mean(abs(bs.ac.c) > abs(ac.c))
  
  structure(class = "gofCOP", 
            list(method = sprintf("Parametric bootstrap goodness-of-fit test (approximate PIOS)"), 
                 erg.tests = matrix(c(test, ac.c), ncol = 2, 
                                    dimnames = list("PIOSTn", c("p.value", "test statistic")))))
}

#################################################################################
gofKernel = function(copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", MJ = 100, delta.J = 0.5, nodes.Integration = 12, execute.times.comp = T){
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (dim(x)[2] != 2){stop("x must be of dimension 2")}
  if (is.element(copula, c("gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofKernel.")}
  dims = 2
  n = dim(x)[1]
  
  if (execute.times.comp == T & M >= 100 | execute.times.comp == T & MJ >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofKernel", M = M, MJ = MJ, print.res = F)
    print(.get.time(times.comp))
  }

  if (copula == "gaussian"){copula = "normal"}
  if (df.est == T){df.fixed = F} else if (df.est == F){df.fixed = T}
  
  if (any(x > 1) || any(x < 0)){
    if (margins == "ranks"){
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " of the observations.", sep = ""))
    } else {
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " distribution.", sep = ""))
    }
    
    x = .margins(x, margins)
  }
  
  if ("normal" == copula || "t" == copula){
    if (param.est == T){
      param = try(fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed), data = x, method = "mpl")@estimate, silent = T)
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed), data = x, method = "itau")@estimate}
    }
    if (copula == "t" & df.fixed == F){
      copula = ellipCopula(copula, param = param[1], dim = dim(x)[2], df = param[2], df.fixed = T)
    } else {
      copula = ellipCopula(copula, param = param, dim = dim(x)[2], df = df, df.fixed = df.fixed)
    }
  } else if("clayton" == copula || "frank" == copula || "gumbel" == copula || 
            "amh" == copula || "joe" == copula){
    if (param.est == T){
      param = try(fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "mpl")@estimate, silent = T)
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "itau")@estimate}
    }
    copula = archmCopula(copula, param = param, dim = dim(x)[2])
  }
  
  bs.ac.c = c()
  ac.c  = .Kernel(x = x, copula = copula, dims = dims, n = n, nodes.Integration = nodes.Integration, MJ = MJ, delta.J = delta.J)
  for(i in 1:M) {
    bs.ac.c[i] = .Kernel(rCopula(n, copula), copula = copula, dims = dims, n = n, nodes.Integration = nodes.Integration, MJ = MJ, delta.J = delta.J)
  }
  test = mean(abs(bs.ac.c) > abs(ac.c))

  structure(class = "gofCOP", 
            list(method = sprintf("Semiparametric bootstrap goodness-of-fit test with Scaillet test"),
                 erg.tests = matrix(c(test, ac.c), ncol = 2, 
                                    dimnames = list("Kernel", c("p.value", "test statistic")))))
}
