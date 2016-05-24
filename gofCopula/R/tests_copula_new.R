########################################################## SnB
gofRosenblattSnB = function (copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", execute.times.comp = T) {
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (is.element(copula, c("gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofRosenblattSnB.")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofRosenblattSnB", M = M, print.res = F)
    print(.get.time(times.comp))
  }
#  .gofRosenblattSnB(copula = copula, x = x, M = M, param = param, param.est = param.est, df = df, df.est = df.est, margins = margins, ...)
  res = try(.gofRosenblattSnB(copula = copula, x = x, M = M, param = param, param.est = param.est, df = df, df.est = df.est, margins = margins), silent = T)
  if (res[1] == "Error in claytonCopula(param, dim = dim, ...) : \n  param can be negative only for dim = 2\n"){
    stop("The dependence parameter is negative for the dataset. For the clayton copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")
  } else if (res[1] == "Error in frankCopula(param, dim = dim, ...) : \n  param can be negative only for dim = 2\n") {
    stop("The dependence parameter is negative for the dataset. For the frank copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")
  } else {
    res
  }
}
.gofRosenblattSnB = function (copula, x, M, param, param.est, df, df.est, margins) {
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
  res = try(.gofCopulapb(copula = copula, x = x, N = M, method = "SnB", estim.method = "mpl"), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); .gofCopulapb(copula = copula, x = x, N = M, method = "SnB", estim.method = "itau")} else {res}
}

########################################################## SnC
gofRosenblattSnC = function (copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", execute.times.comp = T) {
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (is.element(copula, c("gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofRosenblattSnC.")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofRosenblattSnC", M = M, print.res = F)
    print(.get.time(times.comp))
  }
  res = try(.gofRosenblattSnC(copula = copula, x = x, M = M, param = param, param.est = param.est, df = df, df.est = df.est, margins = margins), silent = T)
  if (res[1] == "Error in claytonCopula(param, dim = dim, ...) : \n  param can be negative only for dim = 2\n"){
    stop("The dependence parameter is negative for the dataset. For the clayton copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")
  } else if (res[1] == "Error in frankCopula(param, dim = dim, ...) : \n  param can be negative only for dim = 2\n") {
    stop("The dependence parameter is negative for the dataset. For the frank copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")
  } else {
    res
  }
}
  .gofRosenblattSnC = function (copula, x, M, param, param.est, df, df.est, margins) {
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
  res = try(.gofCopulapb(copula = copula, x = x, N = M, method = "SnC", estim.method = "mpl"), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); .gofCopulapb(copula = copula, x = x, N = M, method = "SnC", estim.method = "itau")} else {res}
  
}

########################################################## AnChisq
gofADChisq = function (copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", execute.times.comp = T) {
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (is.element(copula, c("gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofADChisq.")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofADChisq", M = M, print.res = F)
    print(.get.time(times.comp))
  }
  res = try(.gofADChisq(copula = copula, x = x, M = M, param = param, param.est = param.est, df = df, df.est = df.est, margins = margins), silent = T)
  if (res[1] == "Error in claytonCopula(param, dim = dim, ...) : \n  param can be negative only for dim = 2\n"){
    stop("The dependence parameter is negative for the dataset. For the clayton copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")
  } else if (res[1] == "Error in frankCopula(param, dim = dim, ...) : \n  param can be negative only for dim = 2\n") {
    stop("The dependence parameter is negative for the dataset. For the frank copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")
  } else {
    res
  }
}
.gofADChisq = function (copula, x, M, param, param.est, df, df.est, margins) {
  if (df.est == T){df.fixed = F} else if (df.est == F){df.fixed = T}
  if (copula == "gaussian"){copula = "normal"}
  
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
  res = try(.gofCopulapb(copula = copula, x = x, N = M, method = "AnChisq", estim.method = "mpl"), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); .gofCopulapb(copula = copula, x = x, N = M, method = "AnChisq", estim.method = "itau")} else {res}
}

########################################################## AnGamma
gofADGamma = function (copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", execute.times.comp = T) {
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (is.element(copula, c("gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofADGamma.")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofADGamma", M = M, print.res = F)
    print(.get.time(times.comp))
  }
  res = try(.gofADGamma(copula = copula, x = x, M = M, param = param, param.est = param.est, df = df, df.est = df.est, margins = margins), silent = T)
  if (res[1] == "Error in claytonCopula(param, dim = dim, ...) : \n  param can be negative only for dim = 2\n"){
    stop("The dependence parameter is negative for the dataset. For the clayton copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")
  } else if (res[1] == "Error in frankCopula(param, dim = dim, ...) : \n  param can be negative only for dim = 2\n") {
    stop("The dependence parameter is negative for the dataset. For the frank copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")
  } else {
    res
  }
}
.gofADGamma = function (copula, x, M, param, param.est, df, df.est, margins) {
  if (df.est == T){df.fixed = F} else if (df.est == F){df.fixed = T}
  if (copula == "gaussian"){copula = "normal"}
  
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
  res = try(.gofCopulapb(copula = copula, x = x, N = M, method = "AnGamma", estim.method = "mpl"), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); .gofCopulapb(copula = copula, x = x, N = M, method = "AnGamma", estim.method = "itau")} else {res}
}

########################################################## Rn
gofRn = function (copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, m_b = 0.5, zeta.m = 0, b_Rn = 0.05, execute.times.comp = T) {
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (is.element(copula, c("gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofRn.")}
  if (dim(x)[2] != 2){stop("x must be of dimension 2")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofRn", M = M, print.res = F)
    print(.get.time(times.comp))
  }
  res = try(.gofRn(copula = copula, x = x, M = M, param = param, param.est = param.est, df = df, df.est = df.est, m = m_b, zeta.m = zeta.m, b_Rn = b_Rn), silent = T)
  if (res[1] == "Error in claytonCopula(param, dim = dim, ...) : \n  param can be negative only for dim = 2\n"){
    stop("The dependence parameter is negative for the dataset. For the clayton copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")
  } else if (res[1] == "Error in frankCopula(param, dim = dim, ...) : \n  param can be negative only for dim = 2\n") {
    stop("The dependence parameter is negative for the dataset. For the frank copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")
  } else {
    res
  }
}
.gofRn = function (copula, x, M, param, param.est, df, df.est, m, zeta.m, b_Rn) {  
  if (df.est == T){df.fixed = F} else if (df.est == F){df.fixed = T}
  if (copula == "gaussian"){copula = "normal"}
  
  if ("normal" == copula || "t" == copula){
    if (param.est == T){
      param = fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed), data = pobs(x), method = "itau")@estimate
    }
    if (copula == "t" & df.fixed == F){
      copula = ellipCopula(copula, param = param[1], dim = dim(x)[2], df = param[2], df.fixed = T)
    } else {
      copula = ellipCopula(copula, param = param, dim = dim(x)[2], df = df, df.fixed = df.fixed)
    }
  } else if("clayton" == copula || "frank" == copula || "gumbel" == copula || 
              "amh" == copula || "joe" == copula){
    if (param.est == T){
      param = fitCopula(archmCopula(copula, dim = dim(x)[2]), data = pobs(x), method = "itau")@estimate
    }
    copula = archmCopula(copula, param = param, dim = dim(x)[2])
  }
  erg = gofCopula(copula = copula, x = x, N = M, method = "Rn", estim.method = "itau", simulation = "mult", m = m, zeta.m = zeta.m, b = b_Rn)
  structure(class = "gofCOP", 
            list(method = "Multiplier bootstrap goodness-of-fit test with Rn test",
                 erg.tests = matrix(c(erg$p.value, erg$statistic), ncol = 2, 
                                    dimnames = list("Rn", c("p.value", "test statistic")))))
}

########################################################## Sn
gofSn = function (copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", execute.times.comp = T) {
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (is.element(copula, c("gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofSn.")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofSn", M = M, print.res = F)
    print(.get.time(times.comp))
  }
  res = try(.gofSn(copula = copula, x = x, M = M, param = param, param.est = param.est, df = df, df.est = df.est, margins = margins), silent = T)
  if (res[1] == "Error in claytonCopula(param, dim = dim, ...) : \n  param can be negative only for dim = 2\n"){
    stop("The dependence parameter is negative for the dataset. For the clayton copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")
  } else if (res[1] == "Error in frankCopula(param, dim = dim, ...) : \n  param can be negative only for dim = 2\n") {
    stop("The dependence parameter is negative for the dataset. For the frank copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")
  } else {
    res
  }
}
.gofSn = function (copula, x, M, param, param.est, df, df.est, margins) {
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
  res = try(.gofCopulapb(copula = copula, x = x, N = M, method = "Sn", estim.method = "mpl"), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); .gofCopulapb(copula = copula, x = x, N = M, method = "Sn", estim.method = "itau")} else {res}
}

