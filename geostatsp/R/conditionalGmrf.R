
varOneCell = function(Dcell, DcellEnd, 
		Qchol, cholQptau, param) {
	Ncell = ncol(Qchol)
	Dcell = Dcell:DcellEnd
	# Dcell can be a vector of cells
	thisD = sparseMatrix(Dcell,1:length(Dcell),x=1,
			dims=c(Ncell,length(Dcell)))
	solveQ = solve(Qchol, thisD)
	solveQdiag = diag(solveQ[Dcell,])
	#		solveQp = solve(cholIcQ, solveQ,system='P')
	#		thisDp = solve(cholIcQ, thisD ,system='P')
	solvecqp = 	solve(cholQptau , thisD)
	
	result=param['variance_optimal']*(
				solveQdiag - 
				apply(solvecqp * solveQ, 2, sum)
				) / param['nugget']
	result
}


conditionalGmrf = function(param,
                           Yvec,Xmat, NN,
                           template=NULL, 
                           mc.cores=1, 
						   cellsPerLoop=10,
						   ...) {
  if(is.vector(param))
	  param = as.matrix(param)
					   
  rownames(param) = gsub("sigmasq","variance", rownames(param))
  rownames(param) = gsub("tausq","nugget", rownames(param))
  

  fixed= Xmat %*% param[colnames(Xmat),]
 
  residsOrig = Yvec -	fixed
		
 Q = list(maternGmrfPrec(NN,
                     param=c(param[c('oneminusar','shape'),1], variance=1),
                     ...)
	)
 Qchol = list(Cholesky(Q[[1]], LDL=FALSE))
  for(Dsim in seq(from=2, by=1, len=ncol(param)-1)) {
	 Q[[Dsim]] = maternGmrfPrec(NN,
			 param=c(param[c('oneminusar','shape'),Dsim], variance=1),
			 ...)
	 Qchol[[Dsim]] = update(Qchol[[1]], Q[[Dsim]])
 }	
	
  # Q/sigmasq = var(U)^(-1)
  

  # Qptau = tausq (1/tausq + Q)
  
  
  #	LofQ = expand(Qchol)$L
  #	pRev = as(ncol(LofQ):1, "pMatrix")		
  #	lQLL =  as( t(LofQ %*% pRev) %*% pRev,'dtCMatrix')
  #	QLL = forceSymmetric(tcrossprod(lQLL))
  
  # QLL = P' L  L' P 
  # QLLorig = Prev P' L   L' P Prev'
  # tausq L^(-1) IcQ L^(-1 ') = var(Y)
  
  #	cholIcQ = Cholesky(QLL,LDL=FALSE,perm=TRUE,
  #			Imult=param['variance']/param['nugget'])
  # need to multiply by param['nugget']
  #	ptwice =   as(expand(cholQLL)$P,'sparseMatrix') %*% 
  #			t(as(pRev,'sparseMatrix')) %*% 
  #			as(expand(Qchol)$P,'sparseMatrix')
  #	ptwice2 = as(ptwice, 'pMatrix')
  #	cholIcQ@perm = as.integer(ptwice2@perm-1)
  
  

  EUY = VUY= NULL
  cholQptau = list()
  
  SstartCell = seq(from=1, to=nrow(Xmat), by=cellsPerLoop)

  SstartCell = c(SstartCell, nrow(Xmat)+1)

  SendCell = SstartCell[-1]-1
  SstartCell = SstartCell[-length(SstartCell)]
  
  for(Dsim in 1:ncol(param)) {
	  cholQptau[[Dsim]] = update(Qchol[[Dsim]],
			  parent=Q[[Dsim]],
			  mult=1/param['nugget',Dsim]
	  )
	  
	  EUY = cbind(EUY, 
			  as.vector(
					  solve(cholQptau[[Dsim]], 
							  residsOrig[,Dsim],
							  system='A'
					  )
			  )/param['nugget',Dsim]
	  )
	  moreArgs = list(
			  Qchol=Qchol[[Dsim]], 
			  cholQptau=cholQptau[[Dsim]], 
			  param=param[,Dsim]
	  )
	  if(mc.cores==1) {
		 thediag = mapply(varOneCell, Dcell=SstartCell,
				 DcellEnd=SendCell,
				 MoreArgs = moreArgs,
				 SIMPLIFY=FALSE
		 )
	  } else {
		  thediag = parallel::mcmapply(
				  varOneCell, Dcell=SstartCell,
				  DcellEnd=SendCell,
				  MoreArgs = moreArgs,
				  mc.cores=mc.cores,SIMPLIFY=FALSE)
	  }
	  VUY = cbind(VUY,
			  unlist(thediag)
	  )
  }
  
   
  
  result=abind::abind(random=EUY,krigeSd= sqrt(VUY),
               fixed=fixed,predict=fixed+EUY,
			   resids=residsOrig,
			   along=3)
	if(is.null(dimnames(result)[[2]])) {
		dimnames(result)[[2]] = 
				paste("sim", 1:(dim(result)[[2]]),sep="")
	}
	   
  if(!is.null(template)){
    resRast = raster::brick(raster(template), nl=prod(dim(result)[-1]))
    names(resRast) = as.vector(do.call(
					function(a,b) outer(a,b,paste,sep="_"),
					dimnames(result)[-1] ))
	if(dim(result)[[2]]==1) {
		names(resRast) = dimnames(result)[[3]]
	}
	values(resRast) = as.vector(result)
    result = resRast		
  }
  result
}

