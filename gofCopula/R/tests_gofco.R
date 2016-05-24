gofco = function(copulaobject, x, testset = c("gofPIOSRn", "gofKernel"), margins = "ranks", M = 1000, 
                 execute.times.comp = T, m = 1, MJ = 100, delta.J = 0.5, nodes.Integration = 12, 
                 m_b = 0.5, zeta.m = 0, b_Rn = 0.05) {
  if (is.matrix(x) == F){stop("x must be a matrix")}
  switch(class(copulaobject),
         normalCopula = {copula = "gaussian"; param = copulaobject@parameters; param.est = if (is.na(copulaobject@parameters)){T} else {F}; df = 4; df.est = T},
         tCopula = {copula = "t"; param = copulaobject@parameters[1]; param.est = if (is.na(copulaobject@parameters[1])){T} else {F}; df = copulaobject@parameters[2]; df.est = if (copulaobject@df.fixed == F) {T} else if (copulaobject@df.fixed == T) {F}},
         claytonCopula = {copula = "clayton"; param = copulaobject@parameters; param.est = if (is.na(copulaobject@parameters)){T} else {F}; df = 4; df.est = T},
         frankCopula = {copula = "frank"; param = copulaobject@parameters; param.est = if (is.na(copulaobject@parameters)){T} else {F}; df = 4; df.est = T},
         gumbelCopula = {copula = "gumbel"; param = copulaobject@parameters; param.est = if (is.na(copulaobject@parameters)){T} else {F}; df = 4; df.est = T},
         stop("The class of the object is not supported."))
  gofHybrid(copula = copula, x = x, testset = testset, margins = margins, M = M, execute.times.comp = execute.times.comp, param = param, param.est = param.est, df = df, df.est = df.est, m = m, MJ = MJ, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn)
}