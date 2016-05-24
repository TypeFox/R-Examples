loglikGmrfOneRange = function(
    oneminusar,
    shape, NN,
    adjustEdges, obsCov, xisqTausq,
    reml, jacobian, control){
  
  Q =  geostatsp::maternGmrfPrec(NN,
      param=c(
          shape=as.vector(shape),
          oneminusar=as.vector(oneminusar),
          conditionalVariance=1),
      adjustEdges=adjustEdges)
  
  result = .Call('gmrfLik',
      Q, 
      obsCov, 
      as.double(xisqTausq), 
      reml,
      as.double(jacobian),
      as.double(
          control$logXisqTausq[
              c('lower','upper','tol')]
      )
  )
  
  attributes(result)$param = attributes(Q)$param
  
  result
}  
  
loglikGmrf = function(
    Yvec, Xmat, NN, 
    oneminusar = NULL,
    propNugget=NULL,
    boxcox=1,
    fixBoxcox=TRUE,
    seqBoxcox=outer(c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5), c(-1,1)),
    shape=1,
    reml=TRUE,
    adjustEdges=FALSE,
    jacobian = 0,
    control=list(
        oneminusar = c(lower=1e-3, upper=0.8, tol= .Machine$double.eps^0.25),
        propNugget = c(lower=1e-3, upper=1e3, tol= .Machine$double.eps^0.25)
    ),
    mc.cores=1
) {

  Yvec = as.matrix(Yvec)
	if(is.null(colnames(Yvec)))
		colnames(Yvec) = paste("y", 1:ncol(Yvec), sep='')
  Ny = ncol(Yvec)
	
	if(Ny > 1 & !(length(oneminusar) & length(propNugget)))
		stop("specify oneminusar and propNugget if there is more than one dataset in Y")
	
  Nobs = nrow(Yvec)
  Lcol = c('m2logLml', 'm2logLreml')[reml+1]
  
  # boxcox
  if(!fixBoxcox){
    seqBoxcox = sort(unique(c(0, seqBoxcox)))
    if(Ny>1) warning("need fixBoxcox=TRUE if Y has more than one column")
    bcStuff = forBoxCox(y=Yvec,x=Xmat,seqBoxCox=seqBoxcox)
    Ny = ncol(bcStuff$y)
    Yvec = bcStuff$y
    jacobian = bcStuff$jacobian
    boxcox=bcStuff$boxcox
  } else {
    if(length(jacobian)==1)
      jacobian = rep(jacobian, Ny)
    if(length(boxcox)==1){
      boxcox = rep(boxcox, Ny)
    } 
  }
	
	profBC = NULL # store profile likelihood of box-cox parameter
  
  if(length(jacobian)!= Ny | length(boxcox)!= Ny)
    warning('length of jacobian or boxcox must be same as columns of Y')
  
  obsCov = as.matrix(cbind(Yvec, Xmat))

  if(length(propNugget))
    propNugget = sort(unique(c(0, propNugget)))
  xisqTausq = 1/propNugget 
  
  
  Nxysq = ncol(obsCov)^2
  mlColNames = c(
      'det','detReml','m2logLml', 'm2logLreml',
      'profiledVarianceHatMl', 'profiledVarianceHatReml', 'xisqTausq','junk')
  Ncols = length(mlColNames)*Ny+Nxysq
  
  if(!length(control$xisqTausq)){
    control$xisqTausq = c(
        lower=as.vector(1/control$propNugget['upper']),
        upper=as.vector(1/control$propNugget['lower']),
        control$propNugget['tol']
        )
  }
  
  if(!length(control$logXisqTausq)){
    control$logXisqTausq = c(
        log(control$xisqTausq[c('lower','upper')]),
        control$xisqTausq['tol']
    )
  }


  if(length(oneminusar)){
    # not optimizing
    
  fromC = parallel::mcmapply(
      loglikGmrfOneRange,
      oneminusar=oneminusar,
      MoreArgs=list(
          shape=shape, NN=NN,
      adjustEdges=adjustEdges, 
      obsCov=obsCov, 
      xisqTausq=xisqTausq,
      reml=reml, 
      jacobian=jacobian, control=control),
  mc.cores=mc.cores, SIMPLIFY=FALSE
  )

  parInfo = lapply(fromC, function(qq) attributes(qq)$param)
  fromC = simplify2array(fromC)
  
  NxysqTausq = nrow(fromC)/(8*Ny+ncol(obsCov)^2)
  logLstart = NxysqTausq*Nxysq
  Lseq = 1:logLstart
  
  ml = fromC[-Lseq,]
  ml = array(
      as.vector(ml), 
      dim=c(Ny, NxysqTausq, length(mlColNames), ncol(fromC)),
      dimnames=list(
          colnames(obsCov)[1:Ny], 
          NULL, mlColNames, oneminusar
      )
  )
	if(length(xisqTausq))
	  dimnames(ml)[[2]] = paste('propNugget=',as.character(1/xisqTausq),sep='')

  dimnames(ml)[[3]] = gsub("^junk$", "boxcox",dimnames(ml)[[3]])
  
  ml[,,'boxcox',] = boxcox
	
	
  ssq = fromC[Lseq,]
  ssq=array(ssq, 
      dim=c(ncol(obsCov),ncol(obsCov),NxysqTausq, ncol(fromC)),
      dimnames=list(colnames(obsCov), colnames(obsCov), 
          dimnames(ml)[[2]], oneminusar))
  
  
  covDim = seq(Ny+1, ncol(obsCov))

  if( !length(propNugget)){
    # retain only propNugget with best likelihood
    bestN = apply(ml[,,Lcol,,drop=FALSE],c(1,3,4),which.min)
		
    newml = array(NA, 
        c(dim(ml)[1],1,dim(ml)[-(1:2)]),
        dimnames=c(dimnames(ml)[1],
            list('propNugget'),dimnames(ml)[-(1:2)])
    )
		
    newssq =array(NA, 
        c(dim(ssq)[1:2],1,dim(ssq)[4]),
        dimnames=list(
          	dimnames(ssq)[[1]],dimnames(ssq)[[2]],
            NULL, dimnames(ssq)[[4]])
    )
    
    for(Dbc in 1:dim(newml)[1]){
      for(Dar in 1:dim(newml)[4]){
        newml[Dbc,1,,Dar] = ml[
            ,
            bestN[Dbc,1,Dar],
            ,Dar] 
        newssq[,,1,Dar]=ssq[
          	,,
          	bestN[Dbc,1,Dar],
          	Dar]    
      }
    }
    ml = newml
    ssq=newssq
  } # end no propNugget
	
	
  # retain only the best box-cox parameter for each pair
# of oneminusar and propNugget
  if(length(boxcox)>1 & !fixBoxcox){
    
		
		profBC = cbind(
				boxcox=boxcox,
#				xisqTausq = ml[1,profBCind[,1],'xisqTausq',1],
#				oneminusar = oneminusar[profBCind[,2]], 
				logLstuff = -apply(ml[,,Lcol,,drop=FALSE], 1, min)/2
		)
		colnames(profBC) = gsub("logLstuff", 
				gsub("^m2", "", Lcol), 
				colnames(profBC))

		
      # only keep best box-cox
  	profBCind = arrayInd(
			apply(ml[,,Lcol,,drop=FALSE], 1, which.min),
			dim(ml)[c(2,4)])
	
		bestBcList = apply(ml[,,Lcol,,drop=FALSE], c(2,4), which.min)
		newml = ml[bestBcList[1,1],,,,drop=FALSE]
    newssq = ssq[c(bestBcList[1,1], covDim),
				c(bestBcList[1,1], covDim),,,drop=FALSE]
		
		
		if(length(unique(as.vector(bestBcList)))>=1){
			# no one box-cox parameter maximizes every likelihood
			# retain only the best box-cox for each par of parameters
  		bestBc = array(as.double(bestBcList), dim(bestBcList))
			dimnames(newml)[[1]] = dimnames(newssq)[[1]][1] = 
					dimnames(newssq)[[2]][1] = 'y'
			
      for(Dar in 1:dim(newml)[4]){
        for(Dnugget in 1:dim(newml)[2]){
          thisBc = bestBc[Dnugget,Dar]
          if(!is.na(thisBc)){
            newml[1,Dnugget,,Dar] =
              ml[thisBc,Dnugget,,Dar]
            newssq[,,Dnugget,Dar] = 
                ssq[c(thisBc,covDim),c(thisBc,covDim),Dnugget,Dar]
            }
        }
      }
      ml = newml
      ssq=newssq
			} # end retain best box cox
		
    } # end have boxcox
    
    
    theInf = ml[,,'xisqTausq',,drop=FALSE] != Inf
    
    forml = abind::abind(
        propNugget = 1/ml[,,'xisqTausq',,drop=FALSE],
        tausqHatMl=ml[,,'profiledVarianceHatMl',,drop=FALSE]*theInf,
        tausqHatReml=ml[,,'profiledVarianceHatReml',,drop=FALSE]*theInf,
        xisqHatMl = theInf*ml[,,'profiledVarianceHatMl',,drop=FALSE]*ml[,,'xisqTausq',,drop=FALSE],
        xisqHatReml = theInf*ml[,,'profiledVarianceHatReml',,drop=FALSE]*ml[,,'xisqTausq',,drop=FALSE],
        oneminusar = aperm(
            array(oneminusar, dim(ml)[c(4,2,1)]),3:1
    		),
				boxcox = ml[,,'boxcox',,drop=FALSE],
        along=3
    )
    dimnames(forml)[[3]] = c(
        'propNugget','tausqHatMl', 'tausqHatReml', 
        'xisqHatMl','xisqHatReml', 'oneminusar',
				'boxcox')
    
    Ny = dim(ml)[1]
    betaHat = ssq[-(1:Ny),1:Ny,,,drop=FALSE]
    dimnames(betaHat)[[1]] = paste(
        dimnames(betaHat)[[1]], 'BetaHat',sep=''
    )
    betaHat = aperm(betaHat, c(2,3,1,4))

    
    seBetaHat = apply(ssq[-(1:Ny),-(1:Ny),,,drop=FALSE],c(3,4),diag)
    seBetaHat = array(seBetaHat, c(dim(ml)[1],dim(seBetaHat)),
        dimnames=c(dimnames(ml)[1], dimnames(seBetaHat)))
    
    
    seBetaHat = aperm(seBetaHat, c(1,3,2,4))
    
    seBetaHatMl = seBetaHatReml = seBetaHat
    
    for(D in dimnames(seBetaHat)[[3]]){
      seBetaHatMl[,,D,] = seBetaHatMl[,,D,,drop=FALSE] * 
          ml[,,'profiledVarianceHatMl',,drop=FALSE]
      seBetaHatReml[,,D,] = seBetaHatReml[,,D,,drop=FALSE] * 
          ml[,,'profiledVarianceHatReml',,drop=FALSE]
    }

    dimnames(seBetaHatMl)[[3]] = paste(dimnames(seBetaHatMl)[[3]],'SeMl',sep='')
    dimnames(seBetaHatReml)[[3]] = paste(dimnames(seBetaHatReml)[[3]],'SeReml',sep='')
    
    logLml = -ml[,,c('m2logLreml','m2logLml'),,drop=FALSE]/2
    dimnames(logLml)[[3]] = gsub("^m2","", dimnames(logLml)[[3]])

    
    parMat = lapply(parInfo, function(qq){
          c(
          qq$theo[c('range','shape')],
					sigmasqXisq = as.numeric(qq$theo['variance']), 
          optimalShape=as.numeric(qq$optimalShape['shape']),
          optimalVarianceAndShape=qq$optimal[c('shape','variance')]
      )
        })
    parMat = simplify2array(parMat)
    parArray = aperm(
        array(parMat,c(dim(parMat),dim(ml)[1:2])),
        c(length(dim(parMat))+c(1,2), 1:length(dim(parMat))))
    dimnames(parArray)[[3]] = rownames(parMat)

		sigmasqHat = forml[,,c('xisqHatMl', 'xisqHatReml'),,drop=FALSE] *
				parArray[,,rep('sigmasqXisq',2),,drop=FALSE]
		dimnames(sigmasqHat)[[3]] = gsub("^xi", "sigma", 
				dimnames(sigmasqHat)[[3]])
		
    res = abind::abind(logLml, forml, sigmasqHat,
				betaHat, seBetaHatMl, seBetaHatReml, parArray,
        along=3)
		

		theMle = apply(res[,,gsub("^m2","",Lcol),,drop=FALSE], 
				1, which.max)
		mleIndex = arrayInd(theMle, dim(res)[-c(1,3)])
		
		mle = NULL
		for(D in 1:nrow(mleIndex))
	    mle = cbind(mle, res[D,mleIndex[D,1],,mleIndex[D,2]])
    
		colsToKeep = grep(paste(c('Reml','Ml')[1+reml],'$',sep=''),
				rownames(mle), invert=TRUE)
		mle = mle[colsToKeep,,drop=FALSE]
		mle = mle[grep("^logL", rownames(mle), invert=TRUE),,drop=FALSE]
		rownames(mle) = gsub("(Beta)?(Hat)?(Reml|Ml)?$", "", rownames(mle))

    
    res = list(
        mle=mle, 
        mlArray = res
    )
    
  } else {
    #oneminusar is NULL
	########################
    # optimizing
	######################################################
	resEnv = new.env()
    assign('allPar', list(), envir=resEnv)
    assign('allFromC',  NULL, envir=resEnv)
    assign('oneminusar',  NULL, envir=resEnv)
    
    oneL = function(x){
      
      Q =  geostatsp::maternGmrfPrec(NN,
          param=c(
              shape=as.vector(shape),
              oneminusar=as.vector(x),
              conditionalVariance=1),
          adjustEdges=adjustEdges)
      
      fromC = .Call('gmrfLik',
          Q, 
          obsCov, 
          as.double(xisqTausq), 
          reml,
          as.double(jacobian),
          as.double(
              control$logXisqTausq[
                  c('lower','upper','tol')]
          )
      )
			
      
      assign('allFromC',
          cbind(get('allFromC', resEnv), fromC),
          resEnv)
      assign('oneminusar',
          c(get('oneminusar', resEnv), x),
          resEnv)
      assign('allPar', 
          c(get('allPar', resEnv), 
              list(attributes(Q)$param)),
          resEnv
      )
      
      NxysqTausq = length(fromC)/Ncols
      
      logLstart = NxysqTausq*(
            Nxysq+Ny*(2 + reml)
            )
      # get the likelihood column
      ml = fromC[seq(
              logLstart, by=1, len=NxysqTausq
          )]
      # 
      ml = min(ml, na.rm=TRUE)
      
      ml
    }
    
    optimize(oneL, 
        lower=control$oneminusar['lower'], 
        upper=control$oneminusar['upper'], 
        tol=control$oneminusar['tol']
    )

  oneminusarRes = get('oneminusar', resEnv)
  fromC = as.matrix(get('allFromC', resEnv))
  NxysqTausq = nrow(fromC)/(8*Ny+ncol(obsCov)^2)
  logLstart = NxysqTausq*Nxysq
  Lseq = 1:logLstart
  
  ssq = fromC[Lseq,]
  ssq=array(ssq, 
      dim=c(ncol(obsCov),ncol(obsCov),NxysqTausq, ncol(fromC)),
      dimnames=list(colnames(obsCov), colnames(obsCov), NULL, oneminusarRes))

  ml = fromC[-Lseq,,drop=FALSE]
  ml = array(
      as.vector(ml), 
      dim=c(Ny, NxysqTausq, length(mlColNames), ncol(fromC)),
      dimnames=list(
          colnames(Yvec), 
          NULL, mlColNames, oneminusarRes
      )
  )

  minMl = min(ml[,,Lcol,],na.rm=TRUE)
  mleIndex = arrayInd(
          which((ml[,,Lcol,,drop=FALSE]-minMl) < 10),
          dim(ml)[-3])
  
  mlMat = NULL
  for(D in 1:nrow(mleIndex)){
    mlMat = rbind(mlMat,
        ml[mleIndex[1],mleIndex[2],,mleIndex[3]])
  }
  
  mlMat = cbind(
      mlMat, 
      boxcox=boxcox[mleIndex[1]],
      oneminusar = oneminusarRes[mleIndex[,3]],
      propNugget = 1/mlMat[,'xisqTausq']
  )


  mle = mlMat[which.min(mlMat[,Lcol]),]
  
  
  Qinfo = get('allPar', resEnv)
  thepar = Qinfo[[
      which.min(abs(oneminusarRes-mle['oneminusar']))
  ]]
  

  parMat = lapply(Qinfo, function(qq){
        c(
          	qq$theo[c('range','shape')],
						sigmasqXisq = as.numeric(qq$theo['variance']), 
          	optimalShape=as.numeric(qq$optimalShape['shape']),
          	optimalVarianceAndShape=qq$optimal[c('shape','variance')]
      	)
      })
	
	parMat =
        c(
            thepar$theo[c('range','shape')], 
						sigmasqXisq = as.numeric(thepar$theo['variance']), 
          	optimalShape=as.numeric(thepar$optimalShape['shape']),
          	optimalVarianceAndShape=thepar$optimal[c('shape','variance')]
		)
  
  varMat =   parMat[grep('variance', names(parMat))]
  varMle = varMat * mle[c('varMl','varReml')[1+reml]]
  nuggetMle = c(nugget=as.numeric(
          mle['propNugget'] * mle[c('varMl','varReml')[1+reml]]
      ))
  
  mle = c(mle, 
      varMle, nuggetMle,
      parMat[grep('variance', names(parMat), invert=TRUE)])       
  
  
  
  mleIndex = drop(arrayInd(
      which.min(ml[,,Lcol,,drop=FALSE]-minMl),
      dim(ml)[-3]))
  
  toKeep = c(mleIndex[1], seq(Ny+1, ncol(obsCov)))
  mleSsq = ssq[toKeep[-1],toKeep,
      mleIndex[2],
      mleIndex[3]
  ]
  colnames(mleSsq)[1] = 'mle'
  
  betaHat = mleSsq[,'mle']
  seBetaHat = diag(mleSsq[,-1])
  seBetaHat = seBetaHat * mle[c('varMl','varReml')[1+reml]]
  names(seBetaHat) = paste('stdErr',names(seBetaHat), sep='.')
  
  res = list(
      mle=c(mle[names(mle)!='junk'], 
      betaHat, 
      seBetaHat),
      mlMat = mlMat,
      Qinfo = thepar)
  
}

	res$profileBoxCox = profBC

  res
  
}



summaryGmrfFit= function(x) {
  UseMethod("summaryGmrfFit")  
}

summaryGmrfFit.array = function(x) {
  
  
  x2 = aperm(x,seq(length(dim(x)),1))
  x3 = array(x2, c(prod(dim(x2)[1:2]),dim(x2)[-(1:2)]),
		  dimnames=c(
				  list(as.vector(outer(dimnames(x2)[[1]], dimnames(x2)[[2]], paste))),
				  dimnames(x2)[-(1:2)]))
  if(length(dim(x3))==2){
	  res = summaryGmrfFit.matrix(x3,npar=2)
  } else if(length(dim(x3))==3){
	res=NULL
	for(D  in seq(dim(x3)[2],1)){
		res = abind::abind(summaryGmrfFit.matrix(x3[,D,]), 
				res, along=4)
	}
	dimnames(res)[[4]]=dimnames(x3)[[2]]
  }
	return(res)
}


summaryGmrfFit.matrix = function(x,npar=1) {
  
  if(any(rownames(x)=='logL.ml')){
    x = t(x)
  }
  result = list()
  someL = c('ml','reml')
  for(D1 in someL) {
    D=paste('logL.',D1,sep='')
    MLErow = x[which.max(x[ ,D]),]  
    withinCI = which(-2*(x[,D] - MLErow[D]) < qchisq(0.95, 2))
    
    xsub = x[,grep("^se\\.",colnames(x),invert=TRUE),drop=FALSE]
    parCI= t(apply(
      xsub[withinCI, ,drop=FALSE],
      2,range))
    colnames(parCI) = c('q0.025', 'q0.975')
    
    parMat = cbind(mle=MLErow[colnames(xsub)],parCI,se=NA) 
    notD = someL[someL != D1]
    notD = grep(
      paste("\\.",notD, "$",sep=''),
      rownames(parMat), 
      invert=TRUE,value=TRUE)
    
    parMat = parMat[notD,]
    rownames(parMat) = gsub(paste("\\.", D1, "$",sep=""),
                            "",rownames(parMat))
    
    se = grep("\\.se$", rownames(parMat), value=TRUE)
	betahats = gsub("se$", "betaHat", se)
	
	parMat[betahats,"se"] = parMat[se,"mle"]
	

	parMat[betahats,"q0.025"] = as.numeric(parMat[betahats,"mle"]-2*parMat[betahats,"se"])
	parMat[betahats,"q0.975"] = as.numeric(parMat[betahats,"mle"]+2*parMat[betahats,"se"])
	
 
    parKeep = parMat[!rownames(parMat) %in% se,]
	
    dele = grep("^shape|^sigmasq|^xisq|^tausq",rownames(parKeep))
    parKeep[dele,grep("^q0\\.", colnames(parKeep))] = NA
    result[[D1]] = parKeep
    
  }
  result = abind::abind(result[[1]], result[[2]], along=3)
  dimnames(result)[[3]] = someL
  return(result)
}


