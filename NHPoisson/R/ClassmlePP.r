
setClass('mlePP', representation(call='language', coef='numeric', 
   fullcoef='numeric', vcov='matrix',  min='numeric',
   details='list',  minuslogl='function', nobs='integer', 
   method='character', detailsb='list',npar='integer',inddat='numeric',
	lambdafit='numeric',  LIlambda='numeric', UIlambda='numeric',convergence='integer',
	posE='numeric', covariates='matrix', fixed='list', tit='character',tind='logical', t='numeric'), contains='mle')
