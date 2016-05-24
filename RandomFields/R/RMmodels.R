
## This file has been created automatically by 'rfGenerateModels'.


RMtrend <- function(mean) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('mean') && !is.null(subst <- substitute(mean))) 
	par.model[['mean']] <- CheckArg(mean, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RMtrend', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMtrend <- new('RMmodelgenerator',
	.Data = RMtrend,
	type = c('trend'),
	isotropy = c('parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = FALSE,
	maxdim = Inf,
	vdim = -1
	)



RMplus <- function(C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(C0)) submodels[['C0']] <- C0
  if (hasArg(C1)) submodels[['C1']] <- C1
  if (hasArg(C2)) submodels[['C2']] <- C2
  if (hasArg(C3)) submodels[['C3']] <- C3
  if (hasArg(C4)) submodels[['C4']] <- C4
  if (hasArg(C5)) submodels[['C5']] <- C5
  if (hasArg(C6)) submodels[['C6']] <- C6
  if (hasArg(C7)) submodels[['C7']] <- C7
  if (hasArg(C8)) submodels[['C8']] <- C8
  if (hasArg(C9)) submodels[['C9']] <- C9
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMplus', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMplus <- new('RMmodelgenerator',
	.Data = RMplus,
	type = c('undefined'),
	isotropy = c('parameter dependent'),
	domain = c('single variable', 'kernel'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



RMmult <- function(C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(C0)) submodels[['C0']] <- C0
  if (hasArg(C1)) submodels[['C1']] <- C1
  if (hasArg(C2)) submodels[['C2']] <- C2
  if (hasArg(C3)) submodels[['C3']] <- C3
  if (hasArg(C4)) submodels[['C4']] <- C4
  if (hasArg(C5)) submodels[['C5']] <- C5
  if (hasArg(C6)) submodels[['C6']] <- C6
  if (hasArg(C7)) submodels[['C7']] <- C7
  if (hasArg(C8)) submodels[['C8']] <- C8
  if (hasArg(C9)) submodels[['C9']] <- C9
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMmult', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmult <- new('RMmodelgenerator',
	.Data = RMmult,
	type = c('undefined'),
	isotropy = c('parameter dependent'),
	domain = c('single variable', 'kernel'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



RMS  <- function(phi, var, scale, Aniso, proj, anisoT) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.model[['var']] <- CheckArg(var, subst, TRUE)
  if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.model[['scale']] <- CheckArg(scale, subst, TRUE)
  if (hasArg('anisoT') && !is.null(subst <- substitute(anisoT))) 
	par.model[['anisoT']] <- CheckArg(anisoT, subst, TRUE)
  if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.model[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
  if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.model[['proj']] <- CheckProj(proj, subst)
  
  model <- new('RMmodel', call = cl, name = 'RMS', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMS <- new('RMmodelgenerator',
	.Data = RMS,
	type = c('undefined', 'undefined'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable', 'kernel'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



RMave <- function(phi, A, z, spacetime, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('A') && !is.null(subst <- substitute(A))) 
	par.model[['A']] <- CheckArg(A, subst, TRUE)
  if (hasArg('z') && !is.null(subst <- substitute(z))) 
	par.model[['z']] <- CheckArg(z, subst, TRUE)
  if (hasArg('spacetime') && !is.null(subst <- substitute(spacetime))) 
	par.model[['spacetime']] <- CheckArg(spacetime, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMave', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMave <- new('RMmodelgenerator',
	.Data = RMave,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 10,
	vdim = 1
	)



RMbcw <- function(alpha, beta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
  if (hasArg('beta') && !is.null(subst <- substitute(beta))) 
	par.model[['beta']] <- CheckArg(beta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMbcw', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbcw <- new('RMmodelgenerator',
	.Data = RMbcw,
	type = c('variogram', 'positive definite', 'tail correlation', 'positive definite'),
	isotropy = c('isotropic', 'isotropic', 'isotropic', 'spherical isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'normal mixture',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMbessel <- function(nu, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('nu') && !is.null(subst <- substitute(nu))) 
	par.model[['nu']] <- CheckArg(nu, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMbessel', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbessel <- new('RMmodelgenerator',
	.Data = RMbessel,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMbigneiting <- function(kappa, mu, s, sred12, gamma, cdiag, rhored, c, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('kappa') && !is.null(subst <- substitute(kappa))) 
	par.model[['kappa']] <- CheckArg(kappa, subst, TRUE)
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, TRUE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, TRUE)
  if (hasArg('sred12') && !is.null(subst <- substitute(sred12))) 
	par.model[['sred12']] <- CheckArg(sred12, subst, TRUE)
  if (hasArg('gamma') && !is.null(subst <- substitute(gamma))) 
	par.model[['gamma']] <- CheckArg(gamma, subst, TRUE)
  if (hasArg('cdiag') && !is.null(subst <- substitute(cdiag))) 
	par.model[['cdiag']] <- CheckArg(cdiag, subst, TRUE)
  if (hasArg('rhored') && !is.null(subst <- substitute(rhored))) 
	par.model[['rhored']] <- CheckArg(rhored, subst, TRUE)
  if (hasArg('c') && !is.null(subst <- substitute(c))) 
	par.model[['c']] <- CheckArg(c, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMbigneiting', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbigneiting <- new('RMmodelgenerator',
	.Data = RMbigneiting,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -1,
	vdim = 2
	)



RMbernoulli <- function(phi, threshold, correlation, centred, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('threshold') && !is.null(subst <- substitute(threshold))) 
	par.model[['threshold']] <- CheckArg(threshold, subst, TRUE)
  if (hasArg('correlation') && !is.null(subst <- substitute(correlation))) 
	par.model[['correlation']] <- CheckArg(correlation, subst, TRUE)
  if (hasArg('centred') && !is.null(subst <- substitute(centred))) 
	par.model[['centred']] <- CheckArg(centred, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMbernoulli', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbernoulli <- new('RMmodelgenerator',
	.Data = RMbernoulli,
	type = c('tail correlation'),
	isotropy = c('parameter dependent'),
	domain = c('single variable', 'kernel'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMbiwm <- function(nudiag, nured12, nu, s, cdiag, rhored, c, notinvnu, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('nudiag') && !is.null(subst <- substitute(nudiag))) 
	par.model[['nudiag']] <- CheckArg(nudiag, subst, TRUE)
  if (hasArg('nured12') && !is.null(subst <- substitute(nured12))) 
	par.model[['nured12']] <- CheckArg(nured12, subst, TRUE)
  if (hasArg('nu') && !is.null(subst <- substitute(nu))) 
	par.model[['nu']] <- CheckArg(nu, subst, TRUE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, TRUE)
  if (hasArg('cdiag') && !is.null(subst <- substitute(cdiag))) 
	par.model[['cdiag']] <- CheckArg(cdiag, subst, TRUE)
  if (hasArg('rhored') && !is.null(subst <- substitute(rhored))) 
	par.model[['rhored']] <- CheckArg(rhored, subst, TRUE)
  if (hasArg('c') && !is.null(subst <- substitute(c))) 
	par.model[['c']] <- CheckArg(c, subst, TRUE)
  if (hasArg('notinvnu') && !is.null(subst <- substitute(notinvnu))) 
	par.model[['notinvnu']] <- CheckArg(notinvnu, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMbiwm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbiwm <- new('RMmodelgenerator',
	.Data = RMbiwm,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 2
	)



RMbrownresnick <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMbrownresnick', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbrownresnick <- new('RMmodelgenerator',
	.Data = RMbrownresnick,
	type = c('tail correlation'),
	isotropy = c('parameter dependent'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMbr2bg <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMbr2bg', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbr2bg <- new('RMmodelgenerator',
	.Data = RMbr2bg,
	type = c('positive definite'),
	isotropy = c('parameter dependent'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMbr2eg <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMbr2eg', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbr2eg <- new('RMmodelgenerator',
	.Data = RMbr2eg,
	type = c('positive definite'),
	isotropy = c('parameter dependent'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMcauchy <- function(gamma, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('gamma') && !is.null(subst <- substitute(gamma))) 
	par.model[['gamma']] <- CheckArg(gamma, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMcauchy', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMcauchy <- new('RMmodelgenerator',
	.Data = RMcauchy,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'normal mixture',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMcircular <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMcircular', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMcircular <- new('RMmodelgenerator',
	.Data = RMcircular,
	type = c('tail correlation'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'Gneiting-Schaback class',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 2,
	vdim = 1
	)



RMconstant <- function(M, var) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('M') && !is.null(subst <- substitute(M))) 
	par.model[['M']] <- CheckArg(M, subst, FALSE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
  model <- new('RMmodel', call = cl, name = 'RMconstant', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMconstant <- new('RMmodelgenerator',
	.Data = RMconstant,
	type = c('negative definite'),
	isotropy = c('parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'submodel dependent monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -1
	)



iRMcovariate <- function(norm, c, x, raw, addNA, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(norm)) submodels[['norm']] <- norm
  
  if (hasArg('c') && !is.null(subst <- substitute(c))) 
	par.model[['c']] <- CheckArg(c, subst, TRUE)
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckArg(x, subst, TRUE)
  if (hasArg('raw') && !is.null(subst <- substitute(raw))) 
	par.model[['raw']] <- CheckArg(raw, subst, TRUE)
  if (hasArg('addNA') && !is.null(subst <- substitute(addNA))) 
	par.model[['addNA']] <- CheckArg(addNA, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMcovariate', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

iRMcovariate <- new('RMmodelgenerator',
	.Data = iRMcovariate,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = FALSE,
	maxdim = Inf,
	vdim = -1
	)



iRMfixcov <- function(norm, M, x, raw, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(norm)) submodels[['norm']] <- norm
  
  if (hasArg('M') && !is.null(subst <- substitute(M))) 
	par.model[['M']] <- CheckArg(M, subst, TRUE)
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckArg(x, subst, TRUE)
  if (hasArg('raw') && !is.null(subst <- substitute(raw))) 
	par.model[['raw']] <- CheckArg(raw, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMfixcov', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

iRMfixcov <- new('RMmodelgenerator',
	.Data = iRMfixcov,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = c('kernel'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = FALSE,
	maxdim = Inf,
	vdim = -1
	)



RMcoxisham <- function(phi, mu, D, beta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, TRUE)
  if (hasArg('D') && !is.null(subst <- substitute(D))) 
	par.model[['D']] <- CheckArg(D, subst, TRUE)
  if (hasArg('beta') && !is.null(subst <- substitute(beta))) 
	par.model[['beta']] <- CheckArg(beta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMcoxisham', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMcoxisham <- new('RMmodelgenerator',
	.Data = RMcoxisham,
	type = c('positive definite'),
	isotropy = c('zero-space-isotropic'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 1
	)



RMcubic <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMcubic', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMcubic <- new('RMmodelgenerator',
	.Data = RMcubic,
	type = c('tail correlation'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 1
	)



RMcurlfree <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMcurlfree', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMcurlfree <- new('RMmodelgenerator',
	.Data = RMcurlfree,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RMcutoff <- function(phi, diameter, a, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('diameter') && !is.null(subst <- substitute(diameter))) 
	par.model[['diameter']] <- CheckArg(diameter, subst, TRUE)
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMcutoff', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMcutoff <- new('RMmodelgenerator',
	.Data = RMcutoff,
	type = c('positive definite', 'positive definite'),
	isotropy = c('isotropic', 'spherical isotropic'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 13,
	vdim = 1
	)



RMdagum <- function(beta, gamma, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('beta') && !is.null(subst <- substitute(beta))) 
	par.model[['beta']] <- CheckArg(beta, subst, TRUE)
  if (hasArg('gamma') && !is.null(subst <- substitute(gamma))) 
	par.model[['gamma']] <- CheckArg(gamma, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMdagum', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMdagum <- new('RMmodelgenerator',
	.Data = RMdagum,
	type = c('positive definite', 'tail correlation', 'positive definite'),
	isotropy = c('isotropic', 'isotropic', 'spherical isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'parameter dependent monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMdampedcos <- function(lambda, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('lambda') && !is.null(subst <- substitute(lambda))) 
	par.model[['lambda']] <- CheckArg(lambda, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMdampedcos', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMdampedcos <- new('RMmodelgenerator',
	.Data = RMdampedcos,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -1,
	vdim = 1
	)



RMdewijsian <- function(alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMdewijsian', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMdewijsian <- new('RMmodelgenerator',
	.Data = RMdewijsian,
	type = c('variogram'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMdivfree <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMdivfree', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMdivfree <- new('RMmodelgenerator',
	.Data = RMdivfree,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RMepscauchy <- function(alpha, beta, eps, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
  if (hasArg('beta') && !is.null(subst <- substitute(beta))) 
	par.model[['beta']] <- CheckArg(beta, subst, TRUE)
  if (hasArg('eps') && !is.null(subst <- substitute(eps))) 
	par.model[['eps']] <- CheckArg(eps, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMepscauchy', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMepscauchy <- new('RMmodelgenerator',
	.Data = RMepscauchy,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'normal mixture',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMexp <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMexp', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMexp <- new('RMmodelgenerator',
	.Data = RMexp,
	type = c('tail correlation', 'positive definite'),
	isotropy = c('isotropic', 'spherical isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'completely monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMexponential <- function(phi, n, standardised, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('n') && !is.null(subst <- substitute(n))) 
	par.model[['n']] <- CheckArg(n, subst, TRUE)
  if (hasArg('standardised') && !is.null(subst <- substitute(standardised))) 
	par.model[['standardised']] <- CheckArg(standardised, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMexponential', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMexponential <- new('RMmodelgenerator',
	.Data = RMexponential,
	type = c('positive definite'),
	isotropy = c('parameter dependent'),
	domain = c('single variable', 'kernel'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMschlather <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMschlather', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMschlather <- new('RMmodelgenerator',
	.Data = RMschlather,
	type = c('tail correlation'),
	isotropy = c('parameter dependent'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMfractdiff <- function(a, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMfractdiff', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMfractdiff <- new('RMmodelgenerator',
	.Data = RMfractdiff,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 1,
	vdim = 1
	)



RMflatpower <- function(alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMflatpower', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMflatpower <- new('RMmodelgenerator',
	.Data = RMflatpower,
	type = c('variogram'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'Bernstein',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMfbm <- function(alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMfbm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMfbm <- new('RMmodelgenerator',
	.Data = RMfbm,
	type = c('variogram'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'Bernstein',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMfractgauss <- function(alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMfractgauss', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMfractgauss <- new('RMmodelgenerator',
	.Data = RMfractgauss,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 1,
	vdim = 1
	)



RMgauss <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMgauss', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMgauss <- new('RMmodelgenerator',
	.Data = RMgauss,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'normal mixture',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMgenfbm <- function(alpha, beta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
  if (hasArg('beta') && !is.null(subst <- substitute(beta))) 
	par.model[['beta']] <- CheckArg(beta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMgenfbm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMgenfbm <- new('RMmodelgenerator',
	.Data = RMgenfbm,
	type = c('variogram'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMgencauchy <- function(alpha, beta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
  if (hasArg('beta') && !is.null(subst <- substitute(beta))) 
	par.model[['beta']] <- CheckArg(beta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMgencauchy', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMgencauchy <- new('RMmodelgenerator',
	.Data = RMgencauchy,
	type = c('positive definite', 'tail correlation', 'positive definite'),
	isotropy = c('isotropic', 'isotropic', 'spherical isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'parameter dependent monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMgengneiting <- function(kappa, mu, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('kappa') && !is.null(subst <- substitute(kappa))) 
	par.model[['kappa']] <- CheckArg(kappa, subst, TRUE)
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMgengneiting', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMgengneiting <- new('RMmodelgenerator',
	.Data = RMgengneiting,
	type = c('positive definite', 'positive definite', 'positive definite', 'positive definite', 'positive definite', 'positive definite'),
	isotropy = c('isotropic', 'spherical isotropic', 'isotropic', 'spherical isotropic', 'isotropic', 'spherical isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMgneiting <- function(orig, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('orig') && !is.null(subst <- substitute(orig))) 
	par.model[['orig']] <- CheckArg(orig, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMgneiting', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMgneiting <- new('RMmodelgenerator',
	.Data = RMgneiting,
	type = c('positive definite', 'positive definite'),
	isotropy = c('isotropic', 'spherical isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -1,
	vdim = 1
	)



RMgennsst <- function(phi, psi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(psi)) submodels[['psi']] <- psi
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMgennsst', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMgennsst <- new('RMmodelgenerator',
	.Data = RMgennsst,
	type = c('positive definite', 'positive definite'),
	isotropy = c('parameter dependent', 'symmetric'),
	domain = c('kernel'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMhyperbolic <- function(nu, lambda, delta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('nu') && !is.null(subst <- substitute(nu))) 
	par.model[['nu']] <- CheckArg(nu, subst, TRUE)
  if (hasArg('lambda') && !is.null(subst <- substitute(lambda))) 
	par.model[['lambda']] <- CheckArg(lambda, subst, TRUE)
  if (hasArg('delta') && !is.null(subst <- substitute(delta))) 
	par.model[['delta']] <- CheckArg(delta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMhyperbolic', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMhyperbolic <- new('RMmodelgenerator',
	.Data = RMhyperbolic,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'normal mixture',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMiaco <- function(nu, lambda, delta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('nu') && !is.null(subst <- substitute(nu))) 
	par.model[['nu']] <- CheckArg(nu, subst, TRUE)
  if (hasArg('lambda') && !is.null(subst <- substitute(lambda))) 
	par.model[['lambda']] <- CheckArg(lambda, subst, TRUE)
  if (hasArg('delta') && !is.null(subst <- substitute(delta))) 
	par.model[['delta']] <- CheckArg(delta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMiaco', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMiaco <- new('RMmodelgenerator',
	.Data = RMiaco,
	type = c('positive definite'),
	isotropy = c('space-isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMid <- function(phi, vdim, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('vdim') && !is.null(subst <- substitute(vdim))) 
	par.model[['vdim']] <- CheckArg(vdim, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMid', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMid <- new('RMmodelgenerator',
	.Data = RMid,
	type = c('undefined'),
	isotropy = c('parameter dependent'),
	domain = c('single variable', 'kernel'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



RMkolmogorov <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMkolmogorov', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMkolmogorov <- new('RMmodelgenerator',
	.Data = RMkolmogorov,
	type = c('variogram'),
	isotropy = c('vector-isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 3
	)



RMlgd <- function(alpha, beta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
  if (hasArg('beta') && !is.null(subst <- substitute(beta))) 
	par.model[['beta']] <- CheckArg(beta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMlgd', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMlgd <- new('RMmodelgenerator',
	.Data = RMlgd,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -1,
	vdim = 1
	)



RMmastein <- function(phi, nu, delta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('nu') && !is.null(subst <- substitute(nu))) 
	par.model[['nu']] <- CheckArg(nu, subst, TRUE)
  if (hasArg('delta') && !is.null(subst <- substitute(delta))) 
	par.model[['delta']] <- CheckArg(delta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMmastein', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmastein <- new('RMmodelgenerator',
	.Data = RMmastein,
	type = c('positive definite'),
	isotropy = c('space-isotropic'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMma <- function(phi, alpha, theta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
  if (hasArg('theta') && !is.null(subst <- substitute(theta))) 
	par.model[['theta']] <- CheckArg(theta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMma', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMma <- new('RMmodelgenerator',
	.Data = RMma,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMintexp <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMintexp', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMintexp <- new('RMmodelgenerator',
	.Data = RMintexp,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMmatrix <- function(phi, M, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('M') && !is.null(subst <- substitute(M))) 
	par.model[['M']] <- CheckArg(M, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMmatrix', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmatrix <- new('RMmodelgenerator',
	.Data = RMmatrix,
	type = c('positive definite'),
	isotropy = c('parameter dependent'),
	domain = c('single variable', 'kernel'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RMmatern <- function(nu, notinvnu, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('nu') && !is.null(subst <- substitute(nu))) 
	par.model[['nu']] <- CheckArg(nu, subst, TRUE)
  if (hasArg('notinvnu') && !is.null(subst <- substitute(notinvnu))) 
	par.model[['notinvnu']] <- CheckArg(notinvnu, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMmatern', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmatern <- new('RMmodelgenerator',
	.Data = RMmatern,
	type = c('positive definite', 'tail correlation', 'positive definite'),
	isotropy = c('parameter dependent', 'isotropic', 'spherical isotropic'),
	domain = c('single variable', 'kernel'),
	operator = FALSE,
	monotone = 'submodel dependent monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMmqam <- function(phi, C1, C2, C3, C4, C5, C6, C7, C8, C9, theta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(C1)) submodels[['C1']] <- C1
  if (hasArg(C2)) submodels[['C2']] <- C2
  if (hasArg(C3)) submodels[['C3']] <- C3
  if (hasArg(C4)) submodels[['C4']] <- C4
  if (hasArg(C5)) submodels[['C5']] <- C5
  if (hasArg(C6)) submodels[['C6']] <- C6
  if (hasArg(C7)) submodels[['C7']] <- C7
  if (hasArg(C8)) submodels[['C8']] <- C8
  if (hasArg(C9)) submodels[['C9']] <- C9
  
  if (hasArg('theta') && !is.null(subst <- substitute(theta))) 
	par.model[['theta']] <- CheckArg(theta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMmqam', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmqam <- new('RMmodelgenerator',
	.Data = RMmqam,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RMnatsc <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMnatsc', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMnatsc <- new('RMmodelgenerator',
	.Data = RMnatsc,
	type = c('tail correlation'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMnsst <- function(phi, psi, delta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(psi)) submodels[['psi']] <- psi
  
  if (hasArg('delta') && !is.null(subst <- substitute(delta))) 
	par.model[['delta']] <- CheckArg(delta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMnsst', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMnsst <- new('RMmodelgenerator',
	.Data = RMnsst,
	type = c('positive definite'),
	isotropy = c('space-isotropic'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMnugget <- function(tol, vdim, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('tol') && !is.null(subst <- substitute(tol))) 
	par.model[['tol']] <- CheckArg(tol, subst, TRUE)
  if (hasArg('vdim') && !is.null(subst <- substitute(vdim))) 
	par.model[['vdim']] <- CheckArg(vdim, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMnugget', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMnugget <- new('RMmodelgenerator',
	.Data = RMnugget,
	type = c('tail correlation', 'tail correlation', 'tail correlation'),
	isotropy = c('isotropic', 'earth isotropic', 'spherical isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -2
	)



RMparswm <- function(nudiag, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('nudiag') && !is.null(subst <- substitute(nudiag))) 
	par.model[['nudiag']] <- CheckArg(nudiag, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMparswm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMparswm <- new('RMmodelgenerator',
	.Data = RMparswm,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -1
	)



RMpenta <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMpenta', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMpenta <- new('RMmodelgenerator',
	.Data = RMpenta,
	type = c('positive definite', 'positive definite'),
	isotropy = c('isotropic', 'spherical isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 1
	)



RMaskey <- function(alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMaskey', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMaskey <- new('RMmodelgenerator',
	.Data = RMaskey,
	type = c('positive definite', 'positive definite', 'tail correlation'),
	isotropy = c('isotropic', 'spherical isotropic', 'isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMpower <- function(phi, alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMpower', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMpower <- new('RMmodelgenerator',
	.Data = RMpower,
	type = c('positive definite', 'positive definite'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable', 'kernel'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMprod <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMprod', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMprod <- new('RMmodelgenerator',
	.Data = RMprod,
	type = c('positive definite', 'positive definite', 'positive definite'),
	isotropy = c('symmetric', 'spherical symmetric', 'earth symmetric'),
	domain = c('kernel'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -3
	)



RMqam <- function(phi, C1, C2, C3, C4, C5, C6, C7, C8, C9, theta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(C1)) submodels[['C1']] <- C1
  if (hasArg(C2)) submodels[['C2']] <- C2
  if (hasArg(C3)) submodels[['C3']] <- C3
  if (hasArg(C4)) submodels[['C4']] <- C4
  if (hasArg(C5)) submodels[['C5']] <- C5
  if (hasArg(C6)) submodels[['C6']] <- C6
  if (hasArg(C7)) submodels[['C7']] <- C7
  if (hasArg(C8)) submodels[['C8']] <- C8
  if (hasArg(C9)) submodels[['C9']] <- C9
  
  if (hasArg('theta') && !is.null(subst <- substitute(theta))) 
	par.model[['theta']] <- CheckArg(theta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMqam', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMqam <- new('RMmodelgenerator',
	.Data = RMqam,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMqexp <- function(alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMqexp', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMqexp <- new('RMmodelgenerator',
	.Data = RMqexp,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMschur <- function(phi, M, diag, rhored, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('M') && !is.null(subst <- substitute(M))) 
	par.model[['M']] <- CheckArg(M, subst, TRUE)
  if (hasArg('diag') && !is.null(subst <- substitute(diag))) 
	par.model[['diag']] <- CheckArg(diag, subst, TRUE)
  if (hasArg('rhored') && !is.null(subst <- substitute(rhored))) 
	par.model[['rhored']] <- CheckArg(rhored, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMschur', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMschur <- new('RMmodelgenerator',
	.Data = RMschur,
	type = c('positive definite'),
	isotropy = c('parameter dependent'),
	domain = c('single variable', 'kernel'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



RMdelay <- function(phi, s, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMdelay', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMdelay <- new('RMmodelgenerator',
	.Data = RMdelay,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RMspheric <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMspheric', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMspheric <- new('RMmodelgenerator',
	.Data = RMspheric,
	type = c('tail correlation', 'positive definite'),
	isotropy = c('isotropic', 'spherical isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'Gneiting-Schaback class',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 1
	)



RMstable <- function(alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMstable', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMstable <- new('RMmodelgenerator',
	.Data = RMstable,
	type = c('positive definite', 'tail correlation', 'positive definite'),
	isotropy = c('isotropic', 'isotropic', 'spherical isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'parameter dependent monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMintrinsic <- function(phi, diameter, rawR, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('diameter') && !is.null(subst <- substitute(diameter))) 
	par.model[['diameter']] <- CheckArg(diameter, subst, TRUE)
  if (hasArg('rawR') && !is.null(subst <- substitute(rawR))) 
	par.model[['rawR']] <- CheckArg(rawR, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMintrinsic', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMintrinsic <- new('RMmodelgenerator',
	.Data = RMintrinsic,
	type = c('positive definite', 'positive definite'),
	isotropy = c('isotropic', 'spherical isotropic'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 13,
	vdim = 1
	)



RMstein <- function(nu, z, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('nu') && !is.null(subst <- substitute(nu))) 
	par.model[['nu']] <- CheckArg(nu, subst, TRUE)
  if (hasArg('z') && !is.null(subst <- substitute(z))) 
	par.model[['z']] <- CheckArg(z, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMstein', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMstein <- new('RMmodelgenerator',
	.Data = RMstein,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMstp <- function(xi, phi, S, z, M, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(xi)) submodels[['xi']] <- xi
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('S') && !is.null(subst <- substitute(S))) 
	par.model[['S']] <- CheckArg(S, subst, TRUE)
  if (hasArg('z') && !is.null(subst <- substitute(z))) 
	par.model[['z']] <- CheckArg(z, subst, TRUE)
  if (hasArg('M') && !is.null(subst <- substitute(M))) 
	par.model[['M']] <- CheckArg(M, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMstp', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMstp <- new('RMmodelgenerator',
	.Data = RMstp,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = c('kernel'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 10,
	vdim = 1
	)



RMtbm <- function(phi, fulldim, reduceddim, layers, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('fulldim') && !is.null(subst <- substitute(fulldim))) 
	par.model[['fulldim']] <- CheckArg(fulldim, subst, TRUE)
  if (hasArg('reduceddim') && !is.null(subst <- substitute(reduceddim))) 
	par.model[['reduceddim']] <- CheckArg(reduceddim, subst, TRUE)
  if (hasArg('layers') && !is.null(subst <- substitute(layers))) 
	par.model[['layers']] <- CheckArg(layers, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMtbm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMtbm <- new('RMmodelgenerator',
	.Data = RMtbm,
	type = c('positive definite'),
	isotropy = c('parameter dependent'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -1,
	vdim = -3
	)



RMsum <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMsum', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMsum <- new('RMmodelgenerator',
	.Data = RMsum,
	type = c('negative definite', 'negative definite'),
	isotropy = c('symmetric', 'earth symmetric'),
	domain = c('kernel'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -3
	)



RMvector <- function(phi, a, Dspace, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  if (hasArg('Dspace') && !is.null(subst <- substitute(Dspace))) 
	par.model[['Dspace']] <- CheckArg(Dspace, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMvector', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMvector <- new('RMmodelgenerator',
	.Data = RMvector,
	type = c('positive definite', 'positive definite'),
	isotropy = c('symmetric', 'symmetric'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RMwave <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMwave', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMwave <- new('RMmodelgenerator',
	.Data = RMwave,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 1
	)



RMwhittle <- function(nu, notinvnu, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('nu') && !is.null(subst <- substitute(nu))) 
	par.model[['nu']] <- CheckArg(nu, subst, TRUE)
  if (hasArg('notinvnu') && !is.null(subst <- substitute(notinvnu))) 
	par.model[['notinvnu']] <- CheckArg(notinvnu, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMwhittle', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMwhittle <- new('RMmodelgenerator',
	.Data = RMwhittle,
	type = c('positive definite', 'tail correlation', 'positive definite'),
	isotropy = c('parameter dependent', 'isotropic', 'spherical isotropic'),
	domain = c('single variable', 'kernel'),
	operator = FALSE,
	monotone = 'normal mixture',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMsinepower <- function(alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMsinepower', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMsinepower <- new('RMmodelgenerator',
	.Data = RMsinepower,
	type = c('positive definite'),
	isotropy = c('spherical isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 2,
	vdim = 1
	)



RMmultiquad <- function(delta, tau, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('delta') && !is.null(subst <- substitute(delta))) 
	par.model[['delta']] <- CheckArg(delta, subst, TRUE)
  if (hasArg('tau') && !is.null(subst <- substitute(tau))) 
	par.model[['tau']] <- CheckArg(tau, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMmultiquad', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmultiquad <- new('RMmodelgenerator',
	.Data = RMmultiquad,
	type = c('positive definite'),
	isotropy = c('spherical isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 2,
	vdim = 1
	)



RMangle <- function(angle, lat.angle, ratio, diag) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('angle') && !is.null(subst <- substitute(angle))) 
	par.model[['angle']] <- CheckArg(angle, subst, TRUE)
  if (hasArg('lat.angle') && !is.null(subst <- substitute(lat.angle))) 
	par.model[['lat.angle']] <- CheckArg(lat.angle, subst, TRUE)
  if (hasArg('ratio') && !is.null(subst <- substitute(ratio))) 
	par.model[['ratio']] <- CheckArg(ratio, subst, TRUE)
  if (hasArg('diag') && !is.null(subst <- substitute(diag))) 
	par.model[['diag']] <- CheckArg(diag, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'RMangle', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMangle <- new('RMmodelgenerator',
	.Data = RMangle,
	type = c('shape function'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -1
	)



RMball <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckProj(proj, subst)
  model <- new('RMmodel', call = cl, name = 'RMball', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMball <- new('RMmodelgenerator',
	.Data = RMball,
	type = c('shape function'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMeaxxa <- function(E, A) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('E') && !is.null(subst <- substitute(E))) 
	par.model[['E']] <- CheckArg(E, subst, TRUE)
  if (hasArg('A') && !is.null(subst <- substitute(A))) 
	par.model[['A']] <- CheckArg(A, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'RMeaxxa', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMeaxxa <- new('RMmodelgenerator',
	.Data = RMeaxxa,
	type = c('shape function'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 10,
	vdim = -1
	)



RMetaxxa <- function(E, A, alpha) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('E') && !is.null(subst <- substitute(E))) 
	par.model[['E']] <- CheckArg(E, subst, TRUE)
  if (hasArg('A') && !is.null(subst <- substitute(A))) 
	par.model[['A']] <- CheckArg(A, subst, TRUE)
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'RMetaxxa', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMetaxxa <- new('RMmodelgenerator',
	.Data = RMetaxxa,
	type = c('shape function'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 10,
	vdim = 3
	)



RMtrafo <- function(phi, new) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('new') && !is.null(subst <- substitute(new))) 
	par.model[['new']] <- CheckChar(new, subst, ISONAMES, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'RMtrafo', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMtrafo <- new('RMmodelgenerator',
	.Data = RMtrafo,
	type = c('undefined'),
	isotropy = c('parameter dependent'),
	domain = c('single variable', 'kernel'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -1
	)



RMpolygon <- function(lambda) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('lambda') && !is.null(subst <- substitute(lambda))) 
	par.model[['lambda']] <- CheckArg(lambda, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'RMpolygon', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMpolygon <- new('RMmodelgenerator',
	.Data = RMpolygon,
	type = c('shape function'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 2,
	vdim = 1
	)



RMrational <- function(A, a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('A') && !is.null(subst <- substitute(A))) 
	par.model[['A']] <- CheckArg(A, subst, TRUE)
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'RMrational', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMrational <- new('RMmodelgenerator',
	.Data = RMrational,
	type = c('shape function'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMrotat <- function(speed, phi) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('speed') && !is.null(subst <- substitute(speed))) 
	par.model[['speed']] <- CheckArg(speed, subst, TRUE)
  if (hasArg('phi') && !is.null(subst <- substitute(phi))) 
	par.model[['phi']] <- CheckArg(phi, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'RMrotat', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMrotat <- new('RMmodelgenerator',
	.Data = RMrotat,
	type = c('shape function'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 1
	)



RMrotation <- function(phi) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('phi') && !is.null(subst <- substitute(phi))) 
	par.model[['phi']] <- CheckArg(phi, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'RMrotation', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMrotation <- new('RMmodelgenerator',
	.Data = RMrotation,
	type = c('shape function'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = -1
	)



RMsign <- function(phi, p) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('p') && !is.null(subst <- substitute(p))) 
	par.model[['p']] <- CheckArg(p, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'RMsign', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMsign <- new('RMmodelgenerator',
	.Data = RMsign,
	type = c('shape function'),
	isotropy = c('parameter dependent'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMm2r <- function(phi) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  

  model <- new('RMmodel', call = cl, name = 'RMm2r', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMm2r <- new('RMmodelgenerator',
	.Data = RMm2r,
	type = c('shape function'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 1
	)



RMm3b <- function(phi) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  

  model <- new('RMmodel', call = cl, name = 'RMm3b', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMm3b <- new('RMmodelgenerator',
	.Data = RMm3b,
	type = c('shape function'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 1
	)



RMmps <- function(phi) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  

  model <- new('RMmodel', call = cl, name = 'RMmps', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmps <- new('RMmodelgenerator',
	.Data = RMmps,
	type = c('shape function'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 2,
	vdim = 1
	)



RMtruncsupport <- function(phi, radius) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('radius') && !is.null(subst <- substitute(radius))) 
	par.model[['radius']] <- CheckArg(radius, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'RMtruncsupport', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMtruncsupport <- new('RMmodelgenerator',
	.Data = RMtruncsupport,
	type = c('shape function'),
	isotropy = c('parameter dependent'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RRdeterm <- function(mean) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('mean') && !is.null(subst <- substitute(mean))) 
	par.model[['mean']] <- CheckArg(mean, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RRdeterm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRdeterm <- new('RMmodelgenerator',
	.Data = RRdeterm,
	type = c('distribution family'),
	isotropy = c('cartesian system'),
	domain = c('single variable', 'kernel'),
	operator = FALSE,
	monotone = 'mismatch in monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -3
	)



RRgauss <- function(mu, sd, log) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, FALSE)
  if (hasArg('sd') && !is.null(subst <- substitute(sd))) 
	par.model[['sd']] <- CheckArg(sd, subst, FALSE)
  if (hasArg('log') && !is.null(subst <- substitute(log))) 
	par.model[['log']] <- CheckArg(log, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RRgauss', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRgauss <- new('RMmodelgenerator',
	.Data = RRgauss,
	type = c('distribution family'),
	isotropy = c('cartesian system'),
	domain = c('single variable', 'kernel'),
	operator = FALSE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -1
	)



RRloc <- function(phi, mu, scale, pow) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, FALSE)
  if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.model[['scale']] <- CheckArg(scale, subst, FALSE)
  if (hasArg('pow') && !is.null(subst <- substitute(pow))) 
	par.model[['pow']] <- CheckArg(pow, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RRloc', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRloc <- new('RMmodelgenerator',
	.Data = RRloc,
	type = c('distribution family'),
	isotropy = c('cartesian system'),
	domain = c('single variable', 'kernel'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



RRmcmc <- function(phi, mcmc_n, sigma, normed, maxdensity, rand.loc, gibbs) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('mcmc_n') && !is.null(subst <- substitute(mcmc_n))) 
	par.model[['mcmc_n']] <- CheckArg(mcmc_n, subst, FALSE)
  if (hasArg('sigma') && !is.null(subst <- substitute(sigma))) 
	par.model[['sigma']] <- CheckArg(sigma, subst, FALSE)
  if (hasArg('normed') && !is.null(subst <- substitute(normed))) 
	par.model[['normed']] <- CheckArg(normed, subst, FALSE)
  if (hasArg('maxdensity') && !is.null(subst <- substitute(maxdensity))) 
	par.model[['maxdensity']] <- CheckArg(maxdensity, subst, FALSE)
  if (hasArg('rand.loc') && !is.null(subst <- substitute(rand.loc))) 
	par.model[['rand.loc']] <- CheckArg(rand.loc, subst, FALSE)
  if (hasArg('gibbs') && !is.null(subst <- substitute(gibbs))) 
	par.model[['gibbs']] <- CheckArg(gibbs, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RRmcmc', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRmcmc <- new('RMmodelgenerator',
	.Data = RRmcmc,
	type = c('distribution family'),
	isotropy = c('cartesian system'),
	domain = c('single variable', 'kernel'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -1
	)



RRrectangular <- function(phi, safety, minsteplen, maxsteps, parts, maxit, innermin, outermax, mcmc_n, normed, approx, onesided) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('safety') && !is.null(subst <- substitute(safety))) 
	par.model[['safety']] <- CheckArg(safety, subst, FALSE)
  if (hasArg('minsteplen') && !is.null(subst <- substitute(minsteplen))) 
	par.model[['minsteplen']] <- CheckArg(minsteplen, subst, FALSE)
  if (hasArg('maxsteps') && !is.null(subst <- substitute(maxsteps))) 
	par.model[['maxsteps']] <- CheckArg(maxsteps, subst, FALSE)
  if (hasArg('parts') && !is.null(subst <- substitute(parts))) 
	par.model[['parts']] <- CheckArg(parts, subst, FALSE)
  if (hasArg('maxit') && !is.null(subst <- substitute(maxit))) 
	par.model[['maxit']] <- CheckArg(maxit, subst, FALSE)
  if (hasArg('innermin') && !is.null(subst <- substitute(innermin))) 
	par.model[['innermin']] <- CheckArg(innermin, subst, FALSE)
  if (hasArg('outermax') && !is.null(subst <- substitute(outermax))) 
	par.model[['outermax']] <- CheckArg(outermax, subst, FALSE)
  if (hasArg('mcmc_n') && !is.null(subst <- substitute(mcmc_n))) 
	par.model[['mcmc_n']] <- CheckArg(mcmc_n, subst, FALSE)
  if (hasArg('normed') && !is.null(subst <- substitute(normed))) 
	par.model[['normed']] <- CheckArg(normed, subst, FALSE)
  if (hasArg('approx') && !is.null(subst <- substitute(approx))) 
	par.model[['approx']] <- CheckArg(approx, subst, FALSE)
  if (hasArg('onesided') && !is.null(subst <- substitute(onesided))) 
	par.model[['onesided']] <- CheckArg(onesided, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RRrectangular', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRrectangular <- new('RMmodelgenerator',
	.Data = RRrectangular,
	type = c('distribution family'),
	isotropy = c('cartesian system'),
	domain = c('single variable', 'kernel'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -1
	)



RRspheric <- function(spacedim, balldim, R) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('spacedim') && !is.null(subst <- substitute(spacedim))) 
	par.model[['spacedim']] <- CheckArg(spacedim, subst, FALSE)
  if (hasArg('balldim') && !is.null(subst <- substitute(balldim))) 
	par.model[['balldim']] <- CheckArg(balldim, subst, FALSE)
  if (hasArg('R') && !is.null(subst <- substitute(R))) 
	par.model[['R']] <- CheckArg(R, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RRspheric', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRspheric <- new('RMmodelgenerator',
	.Data = RRspheric,
	type = c('distribution family'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'mismatch in monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 1,
	vdim = 1
	)



RRunif <- function(min, max, normed) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('min') && !is.null(subst <- substitute(min))) 
	par.model[['min']] <- CheckArg(min, subst, FALSE)
  if (hasArg('max') && !is.null(subst <- substitute(max))) 
	par.model[['max']] <- CheckArg(max, subst, FALSE)
  if (hasArg('normed') && !is.null(subst <- substitute(normed))) 
	par.model[['normed']] <- CheckArg(normed, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RRunif', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRunif <- new('RMmodelgenerator',
	.Data = RRunif,
	type = c('distribution family'),
	isotropy = c('cartesian system'),
	domain = c('single variable', 'kernel'),
	operator = FALSE,
	monotone = 'mismatch in monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -1
	)



RMmppplus <- function(C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, p) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(C0)) submodels[['C0']] <- C0
  if (hasArg(C1)) submodels[['C1']] <- C1
  if (hasArg(C2)) submodels[['C2']] <- C2
  if (hasArg(C3)) submodels[['C3']] <- C3
  if (hasArg(C4)) submodels[['C4']] <- C4
  if (hasArg(C5)) submodels[['C5']] <- C5
  if (hasArg(C6)) submodels[['C6']] <- C6
  if (hasArg(C7)) submodels[['C7']] <- C7
  if (hasArg(C8)) submodels[['C8']] <- C8
  if (hasArg(C9)) submodels[['C9']] <- C9
  
  if (hasArg('p') && !is.null(subst <- substitute(p))) 
	par.model[['p']] <- CheckArg(p, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RMmppplus', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmppplus <- new('RMmodelgenerator',
	.Data = RMmppplus,
	type = c('point-shape function'),
	isotropy = c('parameter dependent'),
	domain = c('single variable', 'kernel'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



RPaverage <- function(phi, shape, boxcox, intensity, method) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(shape)) submodels[['shape']] <- shape
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('intensity') && !is.null(subst <- substitute(intensity))) 
	par.model[['intensity']] <- CheckArg(intensity, subst, FALSE)
  if (hasArg('method') && !is.null(subst <- substitute(method))) 
	par.model[['method']] <- CheckArg(method, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPaverage', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPaverage <- new('RMmodelgenerator',
	.Data = RPaverage,
	type = c('method for Gauss process', 'method for Gauss process'),
	isotropy = c('non-dimension-reducing', 'non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = 1
	)



RPcirculant <- function(phi, boxcox, force, mmin, strategy, maxGB, maxmem, tolIm, tolRe, trials, useprimes, dependent, approx_step, approx_maxgrid) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('force') && !is.null(subst <- substitute(force))) 
	par.model[['force']] <- CheckArg(force, subst, FALSE)
  if (hasArg('mmin') && !is.null(subst <- substitute(mmin))) 
	par.model[['mmin']] <- CheckArg(mmin, subst, FALSE)
  if (hasArg('strategy') && !is.null(subst <- substitute(strategy))) 
	par.model[['strategy']] <- CheckArg(strategy, subst, FALSE)
  if (hasArg('maxGB') && !is.null(subst <- substitute(maxGB))) 
	par.model[['maxGB']] <- CheckArg(maxGB, subst, FALSE)
  if (hasArg('maxmem') && !is.null(subst <- substitute(maxmem))) 
	par.model[['maxmem']] <- CheckArg(maxmem, subst, FALSE)
  if (hasArg('tolIm') && !is.null(subst <- substitute(tolIm))) 
	par.model[['tolIm']] <- CheckArg(tolIm, subst, FALSE)
  if (hasArg('tolRe') && !is.null(subst <- substitute(tolRe))) 
	par.model[['tolRe']] <- CheckArg(tolRe, subst, FALSE)
  if (hasArg('trials') && !is.null(subst <- substitute(trials))) 
	par.model[['trials']] <- CheckArg(trials, subst, FALSE)
  if (hasArg('useprimes') && !is.null(subst <- substitute(useprimes))) 
	par.model[['useprimes']] <- CheckArg(useprimes, subst, FALSE)
  if (hasArg('dependent') && !is.null(subst <- substitute(dependent))) 
	par.model[['dependent']] <- CheckArg(dependent, subst, FALSE)
  if (hasArg('approx_step') && !is.null(subst <- substitute(approx_step))) 
	par.model[['approx_step']] <- CheckArg(approx_step, subst, FALSE)
  if (hasArg('approx_maxgrid') && !is.null(subst <- substitute(approx_maxgrid))) 
	par.model[['approx_maxgrid']] <- CheckArg(approx_maxgrid, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPcirculant', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPcirculant <- new('RMmodelgenerator',
	.Data = RPcirculant,
	type = c('method for Gauss process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 13,
	vdim = -3
	)



RPcutoff <- function(phi, boxcox, force, mmin, strategy, maxGB, maxmem, tolIm, tolRe, trials, useprimes, dependent, approx_step, approx_maxgrid, diameter, a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('force') && !is.null(subst <- substitute(force))) 
	par.model[['force']] <- CheckArg(force, subst, FALSE)
  if (hasArg('mmin') && !is.null(subst <- substitute(mmin))) 
	par.model[['mmin']] <- CheckArg(mmin, subst, FALSE)
  if (hasArg('strategy') && !is.null(subst <- substitute(strategy))) 
	par.model[['strategy']] <- CheckArg(strategy, subst, FALSE)
  if (hasArg('maxGB') && !is.null(subst <- substitute(maxGB))) 
	par.model[['maxGB']] <- CheckArg(maxGB, subst, FALSE)
  if (hasArg('maxmem') && !is.null(subst <- substitute(maxmem))) 
	par.model[['maxmem']] <- CheckArg(maxmem, subst, FALSE)
  if (hasArg('tolIm') && !is.null(subst <- substitute(tolIm))) 
	par.model[['tolIm']] <- CheckArg(tolIm, subst, FALSE)
  if (hasArg('tolRe') && !is.null(subst <- substitute(tolRe))) 
	par.model[['tolRe']] <- CheckArg(tolRe, subst, FALSE)
  if (hasArg('trials') && !is.null(subst <- substitute(trials))) 
	par.model[['trials']] <- CheckArg(trials, subst, FALSE)
  if (hasArg('useprimes') && !is.null(subst <- substitute(useprimes))) 
	par.model[['useprimes']] <- CheckArg(useprimes, subst, FALSE)
  if (hasArg('dependent') && !is.null(subst <- substitute(dependent))) 
	par.model[['dependent']] <- CheckArg(dependent, subst, FALSE)
  if (hasArg('approx_step') && !is.null(subst <- substitute(approx_step))) 
	par.model[['approx_step']] <- CheckArg(approx_step, subst, FALSE)
  if (hasArg('approx_maxgrid') && !is.null(subst <- substitute(approx_maxgrid))) 
	par.model[['approx_maxgrid']] <- CheckArg(approx_maxgrid, subst, FALSE)
  if (hasArg('diameter') && !is.null(subst <- substitute(diameter))) 
	par.model[['diameter']] <- CheckArg(diameter, subst, FALSE)
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPcutoff', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPcutoff <- new('RMmodelgenerator',
	.Data = RPcutoff,
	type = c('method for Gauss process', 'method for Gauss process'),
	isotropy = c('non-dimension-reducing', 'non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 13,
	vdim = 1
	)



RPintrinsic <- function(phi, boxcox, force, mmin, strategy, maxGB, maxmem, tolIm, tolRe, trials, useprimes, dependent, approx_step, approx_maxgrid, diameter, rawR) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('force') && !is.null(subst <- substitute(force))) 
	par.model[['force']] <- CheckArg(force, subst, FALSE)
  if (hasArg('mmin') && !is.null(subst <- substitute(mmin))) 
	par.model[['mmin']] <- CheckArg(mmin, subst, FALSE)
  if (hasArg('strategy') && !is.null(subst <- substitute(strategy))) 
	par.model[['strategy']] <- CheckArg(strategy, subst, FALSE)
  if (hasArg('maxGB') && !is.null(subst <- substitute(maxGB))) 
	par.model[['maxGB']] <- CheckArg(maxGB, subst, FALSE)
  if (hasArg('maxmem') && !is.null(subst <- substitute(maxmem))) 
	par.model[['maxmem']] <- CheckArg(maxmem, subst, FALSE)
  if (hasArg('tolIm') && !is.null(subst <- substitute(tolIm))) 
	par.model[['tolIm']] <- CheckArg(tolIm, subst, FALSE)
  if (hasArg('tolRe') && !is.null(subst <- substitute(tolRe))) 
	par.model[['tolRe']] <- CheckArg(tolRe, subst, FALSE)
  if (hasArg('trials') && !is.null(subst <- substitute(trials))) 
	par.model[['trials']] <- CheckArg(trials, subst, FALSE)
  if (hasArg('useprimes') && !is.null(subst <- substitute(useprimes))) 
	par.model[['useprimes']] <- CheckArg(useprimes, subst, FALSE)
  if (hasArg('dependent') && !is.null(subst <- substitute(dependent))) 
	par.model[['dependent']] <- CheckArg(dependent, subst, FALSE)
  if (hasArg('approx_step') && !is.null(subst <- substitute(approx_step))) 
	par.model[['approx_step']] <- CheckArg(approx_step, subst, FALSE)
  if (hasArg('approx_maxgrid') && !is.null(subst <- substitute(approx_maxgrid))) 
	par.model[['approx_maxgrid']] <- CheckArg(approx_maxgrid, subst, FALSE)
  if (hasArg('diameter') && !is.null(subst <- substitute(diameter))) 
	par.model[['diameter']] <- CheckArg(diameter, subst, FALSE)
  if (hasArg('rawR') && !is.null(subst <- substitute(rawR))) 
	par.model[['rawR']] <- CheckArg(rawR, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPintrinsic', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPintrinsic <- new('RMmodelgenerator',
	.Data = RPintrinsic,
	type = c('method for Gauss process', 'method for Gauss process'),
	isotropy = c('non-dimension-reducing', 'non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 13,
	vdim = 1
	)



RPdirect <- function(phi, boxcox, max_variab) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('max_variab') && !is.null(subst <- substitute(max_variab))) 
	par.model[['max_variab']] <- CheckArg(max_variab, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPdirect', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPdirect <- new('RMmodelgenerator',
	.Data = RPdirect,
	type = c('method for Gauss process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -3
	)



RPhyperplane <- function(phi, boxcox, superpos, maxlines, mar_distr, mar_param, additive) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('superpos') && !is.null(subst <- substitute(superpos))) 
	par.model[['superpos']] <- CheckArg(superpos, subst, FALSE)
  if (hasArg('maxlines') && !is.null(subst <- substitute(maxlines))) 
	par.model[['maxlines']] <- CheckArg(maxlines, subst, FALSE)
  if (hasArg('mar_distr') && !is.null(subst <- substitute(mar_distr))) 
	par.model[['mar_distr']] <- CheckArg(mar_distr, subst, FALSE)
  if (hasArg('mar_param') && !is.null(subst <- substitute(mar_param))) 
	par.model[['mar_param']] <- CheckArg(mar_param, subst, FALSE)
  if (hasArg('additive') && !is.null(subst <- substitute(additive))) 
	par.model[['additive']] <- CheckArg(additive, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPhyperplane', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPhyperplane <- new('RMmodelgenerator',
	.Data = RPhyperplane,
	type = c('method for Gauss process', 'method for Gauss process'),
	isotropy = c('non-dimension-reducing', 'non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 2,
	vdim = 1
	)



RPnugget <- function(phi, boxcox, tol, vdim) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('tol') && !is.null(subst <- substitute(tol))) 
	par.model[['tol']] <- CheckArg(tol, subst, FALSE)
  if (hasArg('vdim') && !is.null(subst <- substitute(vdim))) 
	par.model[['vdim']] <- CheckArg(vdim, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPnugget', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPnugget <- new('RMmodelgenerator',
	.Data = RPnugget,
	type = c('method for Gauss process', 'method for Gauss process'),
	isotropy = c('non-dimension-reducing', 'non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -2
	)



RPcoins <- function(phi, shape, boxcox, intensity, method) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(shape)) submodels[['shape']] <- shape
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('intensity') && !is.null(subst <- substitute(intensity))) 
	par.model[['intensity']] <- CheckArg(intensity, subst, FALSE)
  if (hasArg('method') && !is.null(subst <- substitute(method))) 
	par.model[['method']] <- CheckArg(method, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPcoins', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPcoins <- new('RMmodelgenerator',
	.Data = RPcoins,
	type = c('method for Gauss process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = 1
	)



RPsequential <- function(phi, boxcox, max_variables, back_steps, initial) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('max_variables') && !is.null(subst <- substitute(max_variables))) 
	par.model[['max_variables']] <- CheckArg(max_variables, subst, FALSE)
  if (hasArg('back_steps') && !is.null(subst <- substitute(back_steps))) 
	par.model[['back_steps']] <- CheckArg(back_steps, subst, FALSE)
  if (hasArg('initial') && !is.null(subst <- substitute(initial))) 
	par.model[['initial']] <- CheckArg(initial, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPsequential', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPsequential <- new('RMmodelgenerator',
	.Data = RPsequential,
	type = c('method for Gauss process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RPspectral <- function(phi, boxcox, sp_lines, sp_grid, prop_factor, sigma) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('sp_lines') && !is.null(subst <- substitute(sp_lines))) 
	par.model[['sp_lines']] <- CheckArg(sp_lines, subst, FALSE)
  if (hasArg('sp_grid') && !is.null(subst <- substitute(sp_grid))) 
	par.model[['sp_grid']] <- CheckArg(sp_grid, subst, FALSE)
  if (hasArg('prop_factor') && !is.null(subst <- substitute(prop_factor))) 
	par.model[['prop_factor']] <- CheckArg(prop_factor, subst, FALSE)
  if (hasArg('sigma') && !is.null(subst <- substitute(sigma))) 
	par.model[['sigma']] <- CheckArg(sigma, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPspectral', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPspectral <- new('RMmodelgenerator',
	.Data = RPspectral,
	type = c('method for Gauss process', 'method for Gauss process'),
	isotropy = c('non-dimension-reducing', 'non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = 1
	)



RPspecific <- function(phi, boxcox) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPspecific', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPspecific <- new('RMmodelgenerator',
	.Data = RPspecific,
	type = c('method for Gauss process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = -3
	)



RPtbm <- function(phi, boxcox, fulldim, reduceddim, layers, lines, linessimufactor, linesimustep, center, points) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('fulldim') && !is.null(subst <- substitute(fulldim))) 
	par.model[['fulldim']] <- CheckArg(fulldim, subst, FALSE)
  if (hasArg('reduceddim') && !is.null(subst <- substitute(reduceddim))) 
	par.model[['reduceddim']] <- CheckArg(reduceddim, subst, FALSE)
  if (hasArg('layers') && !is.null(subst <- substitute(layers))) 
	par.model[['layers']] <- CheckArg(layers, subst, FALSE)
  if (hasArg('lines') && !is.null(subst <- substitute(lines))) 
	par.model[['lines']] <- CheckArg(lines, subst, FALSE)
  if (hasArg('linessimufactor') && !is.null(subst <- substitute(linessimufactor))) 
	par.model[['linessimufactor']] <- CheckArg(linessimufactor, subst, FALSE)
  if (hasArg('linesimustep') && !is.null(subst <- substitute(linesimustep))) 
	par.model[['linesimustep']] <- CheckArg(linesimustep, subst, FALSE)
  if (hasArg('center') && !is.null(subst <- substitute(center))) 
	par.model[['center']] <- CheckArg(center, subst, FALSE)
  if (hasArg('points') && !is.null(subst <- substitute(points))) 
	par.model[['points']] <- CheckArg(points, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPtbm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPtbm <- new('RMmodelgenerator',
	.Data = RPtbm,
	type = c('method for Gauss process', 'method for Gauss process'),
	isotropy = c('non-dimension-reducing', 'non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RPtrend <- function(phi, boxcox) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPtrend', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPtrend <- new('RMmodelgenerator',
	.Data = RPtrend,
	type = c('method for Gauss process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -3
	)



RPbrorig <- function(phi, tcf, xi, mu, s) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(tcf)) submodels[['tcf']] <- tcf
  
  if (hasArg('xi') && !is.null(subst <- substitute(xi))) 
	par.model[['xi']] <- CheckArg(xi, subst, FALSE)
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, FALSE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPbrorig', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPbrorig <- new('RMmodelgenerator',
	.Data = RPbrorig,
	type = c('method for Brown-Resnick process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = -3
	)



RPbrmixed <- function(phi, tcf, xi, mu, s, meshsize, vertnumber, optim_mixed, optim_mixed_tol, optim_mixed_maxpo, lambda, areamat, variobound) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(tcf)) submodels[['tcf']] <- tcf
  
  if (hasArg('xi') && !is.null(subst <- substitute(xi))) 
	par.model[['xi']] <- CheckArg(xi, subst, FALSE)
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, FALSE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, FALSE)
  if (hasArg('meshsize') && !is.null(subst <- substitute(meshsize))) 
	par.model[['meshsize']] <- CheckArg(meshsize, subst, FALSE)
  if (hasArg('vertnumber') && !is.null(subst <- substitute(vertnumber))) 
	par.model[['vertnumber']] <- CheckArg(vertnumber, subst, FALSE)
  if (hasArg('optim_mixed') && !is.null(subst <- substitute(optim_mixed))) 
	par.model[['optim_mixed']] <- CheckArg(optim_mixed, subst, FALSE)
  if (hasArg('optim_mixed_tol') && !is.null(subst <- substitute(optim_mixed_tol))) 
	par.model[['optim_mixed_tol']] <- CheckArg(optim_mixed_tol, subst, FALSE)
  if (hasArg('optim_mixed_maxpo') && !is.null(subst <- substitute(optim_mixed_maxpo))) 
	par.model[['optim_mixed_maxpo']] <- CheckArg(optim_mixed_maxpo, subst, FALSE)
  if (hasArg('lambda') && !is.null(subst <- substitute(lambda))) 
	par.model[['lambda']] <- CheckArg(lambda, subst, FALSE)
  if (hasArg('areamat') && !is.null(subst <- substitute(areamat))) 
	par.model[['areamat']] <- CheckArg(areamat, subst, FALSE)
  if (hasArg('variobound') && !is.null(subst <- substitute(variobound))) 
	par.model[['variobound']] <- CheckArg(variobound, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPbrmixed', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPbrmixed <- new('RMmodelgenerator',
	.Data = RPbrmixed,
	type = c('method for Brown-Resnick process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = -3
	)



RPbrshifted <- function(phi, tcf, xi, mu, s) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(tcf)) submodels[['tcf']] <- tcf
  
  if (hasArg('xi') && !is.null(subst <- substitute(xi))) 
	par.model[['xi']] <- CheckArg(xi, subst, FALSE)
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, FALSE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPbrshifted', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPbrshifted <- new('RMmodelgenerator',
	.Data = RPbrshifted,
	type = c('method for Brown-Resnick process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = -3
	)



RPbernoulli <- function(phi, stationary_only, threshold) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('stationary_only') && !is.null(subst <- substitute(stationary_only))) 
	par.model[['stationary_only']] <- CheckArg(stationary_only, subst, FALSE)
  if (hasArg('threshold') && !is.null(subst <- substitute(threshold))) 
	par.model[['threshold']] <- CheckArg(threshold, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPbernoulli', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPbernoulli <- new('RMmodelgenerator',
	.Data = RPbernoulli,
	type = c('process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 10,
	vdim = -3
	)



RPbrownresnick <- function(phi, tcf, xi, mu, s) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(tcf)) submodels[['tcf']] <- tcf
  
  if (hasArg('xi') && !is.null(subst <- substitute(xi))) 
	par.model[['xi']] <- CheckArg(xi, subst, FALSE)
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, FALSE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPbrownresnick', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPbrownresnick <- new('RMmodelgenerator',
	.Data = RPbrownresnick,
	type = c('process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = -3
	)



RPgauss <- function(phi, boxcox, stationary_only) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('stationary_only') && !is.null(subst <- substitute(stationary_only))) 
	par.model[['stationary_only']] <- CheckArg(stationary_only, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPgauss', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPgauss <- new('RMmodelgenerator',
	.Data = RPgauss,
	type = c('process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 10,
	vdim = -3
	)



RPpoisson <- function(phi, intensity) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('intensity') && !is.null(subst <- substitute(intensity))) 
	par.model[['intensity']] <- CheckArg(intensity, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPpoisson', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPpoisson <- new('RMmodelgenerator',
	.Data = RPpoisson,
	type = c('process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = -3
	)



RPschlather <- function(phi, tcf, xi, mu, s) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(tcf)) submodels[['tcf']] <- tcf
  
  if (hasArg('xi') && !is.null(subst <- substitute(xi))) 
	par.model[['xi']] <- CheckArg(xi, subst, FALSE)
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, FALSE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPschlather', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPschlather <- new('RMmodelgenerator',
	.Data = RPschlather,
	type = c('process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 10,
	vdim = 1
	)



RPopitz <- function(phi, xi, mu, s, alpha) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('xi') && !is.null(subst <- substitute(xi))) 
	par.model[['xi']] <- CheckArg(xi, subst, FALSE)
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, FALSE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, FALSE)
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPopitz', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPopitz <- new('RMmodelgenerator',
	.Data = RPopitz,
	type = c('process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 10,
	vdim = 1
	)



RPsmith <- function(shape, tcf, xi, mu, s) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(shape)) submodels[['shape']] <- shape
  if (hasArg(tcf)) submodels[['tcf']] <- tcf
  
  if (hasArg('xi') && !is.null(subst <- substitute(xi))) 
	par.model[['xi']] <- CheckArg(xi, subst, FALSE)
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, FALSE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPsmith', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPsmith <- new('RMmodelgenerator',
	.Data = RPsmith,
	type = c('process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = 1
	)



RPchi2 <- function(phi, boxcox, f) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('f') && !is.null(subst <- substitute(f))) 
	par.model[['f']] <- CheckArg(f, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPchi2', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPchi2 <- new('RMmodelgenerator',
	.Data = RPchi2,
	type = c('process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 10,
	vdim = -3
	)



RPt <- function(phi, boxcox, nu) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('nu') && !is.null(subst <- substitute(nu))) 
	par.model[['nu']] <- CheckArg(nu, subst, FALSE)
  
  model <- new('RMmodel', call = cl, name = 'RPt', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPt <- new('RMmodelgenerator',
	.Data = RPt,
	type = c('process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 10,
	vdim = -3
	)



R.minus <- function(a, b, factor) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  if (hasArg('b') && !is.null(subst <- substitute(b))) 
	par.model[['b']] <- CheckArg(b, subst, TRUE)
  if (hasArg('factor') && !is.null(subst <- substitute(factor))) 
	par.model[['factor']] <- CheckArg(factor, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.minus', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.minus <- new('RMmodelgenerator',
	.Data = R.minus,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.plus <- function(a, b, factor) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  if (hasArg('b') && !is.null(subst <- substitute(b))) 
	par.model[['b']] <- CheckArg(b, subst, TRUE)
  if (hasArg('factor') && !is.null(subst <- substitute(factor))) 
	par.model[['factor']] <- CheckArg(factor, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.plus', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.plus <- new('RMmodelgenerator',
	.Data = R.plus,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 1,
	vdim = 1
	)



R.div <- function(a, b, factor) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  if (hasArg('b') && !is.null(subst <- substitute(b))) 
	par.model[['b']] <- CheckArg(b, subst, TRUE)
  if (hasArg('factor') && !is.null(subst <- substitute(factor))) 
	par.model[['factor']] <- CheckArg(factor, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.div', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.div <- new('RMmodelgenerator',
	.Data = R.div,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 1,
	vdim = 1
	)



R.mult <- function(a, b, factor) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  if (hasArg('b') && !is.null(subst <- substitute(b))) 
	par.model[['b']] <- CheckArg(b, subst, TRUE)
  if (hasArg('factor') && !is.null(subst <- substitute(factor))) 
	par.model[['factor']] <- CheckArg(factor, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.mult', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.mult <- new('RMmodelgenerator',
	.Data = R.mult,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 1,
	vdim = 1
	)



R.const <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.const', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.const <- new('RMmodelgenerator',
	.Data = R.const,
	type = c('shape function', 'trend', 'tail correlation'),
	isotropy = c('parameter dependent', 'parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.p <- function(proj, new, factor) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.model[['proj']] <- CheckArg(proj, subst, TRUE)
  if (!(hasArg('new') && !is.null(subst <- substitute(new)))) new <- UNREDUCED
	par.model[['new']] <- CheckChar(new, subst, ISONAMES, TRUE)
  if (hasArg('factor') && !is.null(subst <- substitute(factor))) 
	par.model[['factor']] <- CheckArg(factor, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.p', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.p <- new('RMmodelgenerator',
	.Data = R.p,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



R.c <- function(a, b, c, d, e, f, g, h, i, j, factor) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  if (hasArg('b') && !is.null(subst <- substitute(b))) 
	par.model[['b']] <- CheckArg(b, subst, TRUE)
  if (hasArg('c') && !is.null(subst <- substitute(c))) 
	par.model[['c']] <- CheckArg(c, subst, TRUE)
  if (hasArg('d') && !is.null(subst <- substitute(d))) 
	par.model[['d']] <- CheckArg(d, subst, TRUE)
  if (hasArg('e') && !is.null(subst <- substitute(e))) 
	par.model[['e']] <- CheckArg(e, subst, TRUE)
  if (hasArg('f') && !is.null(subst <- substitute(f))) 
	par.model[['f']] <- CheckArg(f, subst, TRUE)
  if (hasArg('g') && !is.null(subst <- substitute(g))) 
	par.model[['g']] <- CheckArg(g, subst, TRUE)
  if (hasArg('h') && !is.null(subst <- substitute(h))) 
	par.model[['h']] <- CheckArg(h, subst, TRUE)
  if (hasArg('i') && !is.null(subst <- substitute(i))) 
	par.model[['i']] <- CheckArg(i, subst, TRUE)
  if (hasArg('j') && !is.null(subst <- substitute(j))) 
	par.model[['j']] <- CheckArg(j, subst, TRUE)
  if (hasArg('factor') && !is.null(subst <- substitute(factor))) 
	par.model[['factor']] <- CheckArg(factor, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.c', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.c <- new('RMmodelgenerator',
	.Data = R.c,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 1,
	vdim = 1
	)



R.is <- function(a, is, b) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  if (hasArg('is') && !is.null(subst <- substitute(is))) 
	par.model[['is']] <- CheckChar(is, subst, EQNAMES, TRUE)
  if (hasArg('b') && !is.null(subst <- substitute(b))) 
	par.model[['b']] <- CheckArg(b, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.is', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.is <- new('RMmodelgenerator',
	.Data = R.is,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 1,
	vdim = 1
	)



R.acos <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.acos', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.acos <- new('RMmodelgenerator',
	.Data = R.acos,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.asin <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.asin', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.asin <- new('RMmodelgenerator',
	.Data = R.asin,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.atan <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.atan', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.atan <- new('RMmodelgenerator',
	.Data = R.atan,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.atan2 <- function(a, b) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  if (hasArg('b') && !is.null(subst <- substitute(b))) 
	par.model[['b']] <- CheckArg(b, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.atan2', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.atan2 <- new('RMmodelgenerator',
	.Data = R.atan2,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.cos <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.cos', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.cos <- new('RMmodelgenerator',
	.Data = R.cos,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.sin <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.sin', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.sin <- new('RMmodelgenerator',
	.Data = R.sin,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.tan <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.tan', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.tan <- new('RMmodelgenerator',
	.Data = R.tan,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.acosh <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.acosh', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.acosh <- new('RMmodelgenerator',
	.Data = R.acosh,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.asinh <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.asinh', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.asinh <- new('RMmodelgenerator',
	.Data = R.asinh,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.atanh <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.atanh', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.atanh <- new('RMmodelgenerator',
	.Data = R.atanh,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.cosh <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.cosh', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.cosh <- new('RMmodelgenerator',
	.Data = R.cosh,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.sinh <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.sinh', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.sinh <- new('RMmodelgenerator',
	.Data = R.sinh,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.tanh <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.tanh', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.tanh <- new('RMmodelgenerator',
	.Data = R.tanh,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.exp <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.exp', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.exp <- new('RMmodelgenerator',
	.Data = R.exp,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.log <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.log', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.log <- new('RMmodelgenerator',
	.Data = R.log,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.expm1 <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.expm1', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.expm1 <- new('RMmodelgenerator',
	.Data = R.expm1,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.log1p <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.log1p', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.log1p <- new('RMmodelgenerator',
	.Data = R.log1p,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.logb <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.logb', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.logb <- new('RMmodelgenerator',
	.Data = R.logb,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.exp2 <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.exp2', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.exp2 <- new('RMmodelgenerator',
	.Data = R.exp2,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.log2 <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.log2', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.log2 <- new('RMmodelgenerator',
	.Data = R.log2,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.pow <- function(a, b) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  if (hasArg('b') && !is.null(subst <- substitute(b))) 
	par.model[['b']] <- CheckArg(b, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.pow', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.pow <- new('RMmodelgenerator',
	.Data = R.pow,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.sqrt <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.sqrt', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.sqrt <- new('RMmodelgenerator',
	.Data = R.sqrt,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.hypot <- function(a, b) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  if (hasArg('b') && !is.null(subst <- substitute(b))) 
	par.model[['b']] <- CheckArg(b, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.hypot', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.hypot <- new('RMmodelgenerator',
	.Data = R.hypot,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.cbrt <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.cbrt', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.cbrt <- new('RMmodelgenerator',
	.Data = R.cbrt,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.ceil <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.ceil', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.ceil <- new('RMmodelgenerator',
	.Data = R.ceil,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.fabs <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.fabs', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.fabs <- new('RMmodelgenerator',
	.Data = R.fabs,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.floor <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.floor', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.floor <- new('RMmodelgenerator',
	.Data = R.floor,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.fmod <- function(a, b) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  if (hasArg('b') && !is.null(subst <- substitute(b))) 
	par.model[['b']] <- CheckArg(b, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.fmod', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.fmod <- new('RMmodelgenerator',
	.Data = R.fmod,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.nearbyint <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.nearbyint', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.nearbyint <- new('RMmodelgenerator',
	.Data = R.nearbyint,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.round <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.round', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.round <- new('RMmodelgenerator',
	.Data = R.round,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.trunc <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.trunc', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.trunc <- new('RMmodelgenerator',
	.Data = R.trunc,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.lrint <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.lrint', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.lrint <- new('RMmodelgenerator',
	.Data = R.lrint,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.llrint <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.llrint', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.llrint <- new('RMmodelgenerator',
	.Data = R.llrint,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.lround <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.lround', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.lround <- new('RMmodelgenerator',
	.Data = R.lround,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.llround <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.llround', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.llround <- new('RMmodelgenerator',
	.Data = R.llround,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.copysign <- function(a, b) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  if (hasArg('b') && !is.null(subst <- substitute(b))) 
	par.model[['b']] <- CheckArg(b, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.copysign', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.copysign <- new('RMmodelgenerator',
	.Data = R.copysign,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.erf <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.erf', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.erf <- new('RMmodelgenerator',
	.Data = R.erf,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.erfc <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.erfc', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.erfc <- new('RMmodelgenerator',
	.Data = R.erfc,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.tgamma <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.tgamma', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.tgamma <- new('RMmodelgenerator',
	.Data = R.tgamma,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.lgamma <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.lgamma', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.lgamma <- new('RMmodelgenerator',
	.Data = R.lgamma,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.rint <- function(a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.rint', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.rint <- new('RMmodelgenerator',
	.Data = R.rint,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.nextafter <- function(a, b) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  if (hasArg('b') && !is.null(subst <- substitute(b))) 
	par.model[['b']] <- CheckArg(b, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.nextafter', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.nextafter <- new('RMmodelgenerator',
	.Data = R.nextafter,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.nexttoward <- function(a, b) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  if (hasArg('b') && !is.null(subst <- substitute(b))) 
	par.model[['b']] <- CheckArg(b, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.nexttoward', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.nexttoward <- new('RMmodelgenerator',
	.Data = R.nexttoward,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.remainder <- function(a, b) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  if (hasArg('b') && !is.null(subst <- substitute(b))) 
	par.model[['b']] <- CheckArg(b, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.remainder', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.remainder <- new('RMmodelgenerator',
	.Data = R.remainder,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.fdim <- function(a, b) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  if (hasArg('b') && !is.null(subst <- substitute(b))) 
	par.model[['b']] <- CheckArg(b, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.fdim', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.fdim <- new('RMmodelgenerator',
	.Data = R.fdim,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.fmax <- function(a, b) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  if (hasArg('b') && !is.null(subst <- substitute(b))) 
	par.model[['b']] <- CheckArg(b, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.fmax', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.fmax <- new('RMmodelgenerator',
	.Data = R.fmax,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.fmin <- function(a, b) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  if (hasArg('b') && !is.null(subst <- substitute(b))) 
	par.model[['b']] <- CheckArg(b, subst, TRUE)
  
  model <- new('RMmodel', call = cl, name = 'R.fmin', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.fmin <- new('RMmodelgenerator',
	.Data = R.fmin,
	type = c('shape function', 'trend'),
	isotropy = c('parameter dependent', 'parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



