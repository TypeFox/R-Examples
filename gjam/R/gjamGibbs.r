
gjamGibbs <-
function(formula, xdata, y, modelList){
  
  # xdata        - n by Q design data.frame
  # y            - n by S response, continuous and discrete depending on TYPE
  # ng           - no. iterations
  # maxBreaks    - if too many classes, aggregated to this value
  # holdoutN     - number to predict out-of-sample
  # holdoutIndex - which to predict out-of-sample
  # censor       - each list within list has a typeName followed by 
  #                $columns (vector of indices) and $partition (3 by K matrix)
  
  holdoutN      <-  0
  holdoutIndex  <- numeric(0)
  thetaPrior    <- NULL
  betaPrior     <- NULL
  pgPrior       <- NULL
  traitMat      <- traitTypes <- traitPart <- NULL
  modelSummary  <- NULL
  censor        <- NULL
  effort        <- NULL
  breakList     <- NULL
  standardX     <- NULL
  traitList     <- NULL
  pg            <- NULL
  ng            <- 2000
  burnin        <- 500
  flist         <- NULL
  ZEROINFL      <- F
  
  for(k in 1:length(modelList))assign( names(modelList)[k], modelList[[k]] )
  if(!is.null(traitList))for(k in 1:length(traitList))assign( names(traitList)[k], traitList[[k]] )
  
  S <- ncol(y)
  n <- nrow(y)
  
  if(is.data.frame(y))y <- as.matrix(y)
  
  if(length(typeNames) == 1)typeNames <- rep(typeNames,S)
  if(length(typeNames) != S) stop('typeNames must be one value or no. columns in y')
  if(!is.numeric(y))         stop('y must be numeric matrix')
  if(burnin >= ng)           stop( 'burnin must be > no. MCMC steps, ng' )
  if('censor' %in% names(modelList)){
    for(k in 1:length(censor)){
      if(!names(censor)[[k]] %in% c('CA','DA'))stop('censor name must be CA or DA')
      if( nrow(censor[[k]]$partition) != 3 )stop('censor matrix must have 3 rows for value, lo, hi')
      rownames(censor[[k]]$partition) <- c('value','lo','hi')
    }
  }
      
  tmp      <- .gjamGetTypes(typeNames)
  typeCols <- tmp$typeCols
  typeFull <- tmp$typeFull
  typeCode <- tmp$TYPES[typeCols]
  allTypes <- sort(unique(typeCols))
  names(typeCols) <- typeCode
  
  tvec  <- paste0(sort(unique(typeCode)),collapse=", ")
  
  facNames <- character(0)
  for(j in 1:ncol(xdata))if(is.factor(xdata[,j]))facNames <- c(facNames,colnames(xdata)[j])
  
  x <- model.frame(formula,data=xdata,na.action=NULL)
  x <- model.matrix(formula,data=x)
  colnames(x)[1] <- 'intercept'
  q <- Q <- ncol(x)
  
  tmp <- .gjamMissingValues(x,y)
    xmiss  <- tmp$xmiss; xbound <- tmp$xbound; missX <- tmp$missX
    missX2 <- tmp$missX2; ymiss <- tmp$ymiss; missY <- tmp$missY; xprior <- tmp$xprior
    yprior <- tmp$yprior
    nmiss  <- length(xmiss)
    mmiss  <- length(ymiss)
    if(nmiss > 0)x[xmiss] <- xprior
  
  tmp <- .gjamXY(x,y,typeCode,standardX,facNames)
    x <- tmp$x; y <- tmp$y; snames <- tmp$snames; xnames <- tmp$xnames
    isInt        <- tmp$isInt; intMat <- tmp$intMat; isSquare  <- tmp$isSquare
    factorList   <- tmp$factorList; isFactor <- tmp$isFactor; isNonLinX <- tmp$isNonLinX
    designTable  <- tmp$designTable; xscale <- tmp$xscale; predXcols <- tmp$predXcols
    modelSummary <- append(modelSummary,list(designTable = designTable, isFactor = isFactor))
  
  updateBeta <- .gjamUpdateBetaNoPrior
  loBeta <- hiBeta <- NULL
  
  if(!is.null(betaPrior)){
    loBeta     <- as.vector(betaPrior$lo)
    hiBeta     <- as.vector(betaPrior$hi)
    updateBeta <- .gjamUpdateBetaPrior
  }                 
  
  tmp <- .gjamHoldoutSetup(holdoutIndex,holdoutN,n)
    holdoutIndex <- tmp$holdoutIndex; holdoutN <- tmp$holdoutN
    inSamples <- tmp$inSamples; nIn <- tmp$nIn
  
  # w, z, cuts, cutLo, cutHi, plo, phi, ordCols, disCols, compCols, breakMat, minOrd, maxOrd
  tmp <- .gjamSetup(typeCols,x,y,breakList,holdoutN,holdoutIndex,censor=censor,effort=effort) 
    w <- tmp$w; z <- tmp$z; y <- tmp$y; other <- tmp$other; cuts <- tmp$cuts
    cutLo <- tmp$cutLo; cutHi <- tmp$cutHi; plo <- tmp$plo; phi <- tmp$phi
    ordCols <- tmp$ordCols; disCols <- tmp$disCols; compCols <- tmp$compCols 
    classBySpec <- tmp$classBySpec; breakMat <- tmp$breakMat
    minOrd <- tmp$minOrd; maxOrd <- tmp$maxOrd; censorCA <- tmp$censorCA
    censorDA <- tmp$censorDA; ncut <- ncol(cuts)
    
    modelSummary <- append(modelSummary,list(classBySpec = classBySpec))
  
  tg       <- cutg <- cuts
  sigmaDf  <- nIn - Q + S - 1
  
  if('x' %in% names(thetaPrior))Q1 <- ncol(thetaPrior$x)
  
  #initial values
  XX  <- crossprod(x[inSamples,])
  IXX <- solve(XX)
  WX  <- crossprod(x[inSamples,],w[inSamples,])
  WIX <- IXX%*%WX
  
  tmp <- updateBeta(WIX=WIX,IXX=IXX,sg=diag(.1,S),w=w,thetaPrior=thetaPrior,tau=rep(-10,S),
                    z = z,bzero=matrix(0,Q1,S),inSamples=inSamples,
                    alpha=matrix(0,Q,S),loBeta=loBeta,hiBeta=hiBeta)
  bg <- alpha <- tmp$bg
  betaZero <- tmp$bzero
  tau      <- tmp$tau
  sg       <- crossprod( w - x%*%bg )/(n - q)
  
  rownames(sg) <- colnames(sg) <- snames
  colnames(x)  <- xnames
  
  rgibbs <- matrix(0,ng,S*S)
  bgibbs <- matrix(0,ng,S*q)
  
  colnames(rgibbs) <- .multivarChainNames(snames,snames)
  sgibbs <- rgibbs
  colnames(bgibbs) <- .multivarChainNames(xnames,snames)
  dgibbs <- bgibbs

  pg <- .9
  
  sensGibbs <- matrix(0,ng,Q)   # sensitivity
  colnames(sensGibbs) <- xnames
  
  TRAITS <- F
  if(length(traitMat) > 0){
    TRAITS <- T
    traitMat <- traitMat[colnames(y),]
    tnames   <- colnames(traitMat)
    M        <- ncol(traitMat)
    traitMat <- t(traitMat)
    
    agibbs <- matrix(0,ng,M*q)
    mgibbs <- matrix(0,ng,M*M)
    tpred  <- tpred2 <- matrix(0,n,M)
    colnames(agibbs) <- .multivarChainNames(xnames,tnames)
    colnames(mgibbs) <- .multivarChainNames(tnames,tnames)
  }
  
  if('OC' %in% typeCode){
    cnames       <- paste('C',1:ncut,sep='-')
    nor          <- length(ordCols)
    cgibbs <- matrix(0,ng,(ncut-3)*nor)
    cutg <- cuts
    colnames(cgibbs) <- as.vector( outer(snames[ordCols],cnames[-c(1,2,ncut)],paste,sep='_') )
    tmp   <- .gjamGetCuts(z,ordCols)
    cutLo <- tmp$cutLo
    cutHi <- tmp$cutHi
    plo[,ordCols] <- tg[cutLo]                                        
    phi[,ordCols] <- tg[cutHi]
    lastOrd <- max(maxOrd) + 1
  }
  
  mm <- apply(z,2,range)
  rownames(mm) <- c('lo bin','hi bin')
 # print(mm)
  
  ycount <- rowSums(y)
  if('CC' %in% typeCode)ycount <- rowSums(y[,compCols])
  
  priorXIV <- diag(1e-5,ncol(x))
  priorX   <- colMeans(x)
  predx    <- predx2 <- x*0
  xpred    <- x
  xpred[,isNonLinX] <- 0            #all nonlinear terms (e.g., interactions)
  px <- 1:Q
  if(!is.null(isNonLinX))px <- px[-isNonLinX]
  px <- px[!xnames[px] %in% isFactor]
  px       <- px[px != 1]
  
  ypred  <- ypred2 <- wpred  <- wpred2 <- sumb <- ymissPred <- ymissPred2 <- y*0
  sumDev <- 0   #for DIC
  sMean  <- sg*0
  ntot   <- 0
  
  pbar <- txtProgressBar(min=1,max=ng,style=1)

  richness <- NULL
  notOther <- c(1:S)
  if(length(other) > 0){                     
    notOther <- notOther[notOther != other]
    sg[other,] <- sg[,other] <- 0
    sg[other,other] <- 1
  }
  
  if(is.null(pgPrior)){
    pgMu <- .97
    p1 <- n*S
    p2 <- p1*(1/pgMu - 1)
    pgPrior <- c(p1,p2)
  }

  for(g in 1:ng){
    
    XX  <- crossprod(x[inSamples,])
    IXX <- solve(XX)
    WX  <- crossprod(x[inSamples,],w[inSamples,notOther])                                #covariance
    WIX <- IXX%*%WX
    
    sg[notOther,notOther] <- .updateWishartNoPrior( x[inSamples,],w[inSamples,notOther],sigmaDf,
                                                   beta=bg[,notOther],
                                IXX=IXX,WX=WX,WIX=WIX,TRYPRIOR=T)$sigma
        
    tmp <- updateBeta(WIX=WIX,IXX=IXX,sg=sg[notOther,notOther],w=w[,notOther],tau=tau,
                      thetaPrior=thetaPrior,
                      z = z[,notOther],
                      bzero=betaZero[,notOther],inSamples=inSamples,
                      alpha=alpha[,notOther],loBeta=loBeta,hiBeta=hiBeta)
    bg[,notOther] <- alpha[,notOther] <- tmp$bg 
    betaZero[,notOther] <- tmp$bzero
    if(!is.null(tmp$tau))tau[,notOther]  <- tmp$tau
    
    muw   <- x%*%bg
    alpha <- .gjamA2B(bg,sg) # bg is covariance scale, alpha is correlation scale
    
    if( 'OC' %in% typeCode ){
      
      tg   <- .gjamUpdateTheta(w,cutg,cutLo,cutHi,ordCols,holdoutN,holdoutIndex,minOrd,maxOrd) # variance scale
      cutg <- .gjamCuts2theta(tg,ss = sg[ordCols,ordCols])                                     # correlation scale
      breakMat[ordCols,1:lastOrd] <- cutg
      cgibbs[g,] <- as.vector( cutg[,-c(1,2,ncut)] )
      
      plo[,ordCols] <- cutg[cutLo]                                        
      phi[,ordCols] <- cutg[cutHi]
    }
    
    tmp <- .gjamUpdateW(w,muw,x,y,sg,alpha,cutg,plo,phi,effort = effort,
                        typeNames,holdoutN, holdoutIndex, censorCA, censorDA,
                        notOther=notOther,pg=pg,pgPrior = pgPrior,breakMat = breakMat)
    w  <- tmp$w
    yp <- tmp$yp
    pg <- tmp$pg
    
    setTxtProgressBar(pbar,g)
  
    if(mmiss > 0)y[ymiss] <- yp[ymiss]
    
    if(nmiss > 0){
      sinv     <- solve(sg)
      x[xmiss] <- .imputX_MVN(x,w,bg,xmiss,sinv,xprior=xprior,xbound=xbound)[xmiss]
      XX  <- crossprod(x)
      IXX <- solve(XX)
      missX <- missX + x[xmiss]
      missX2 <- missX2 + x[xmiss]^2
    }
    
    cg <- .cov2Cor(sg)
    
    rgibbs[g,] <- cg
    bgibbs[g,] <- bg
    dgibbs[g,] <- alpha
    sgibbs[g,] <- sg
    sensGibbs[g,] <- diag( alpha[,notOther]%*%cg[notOther,notOther]%*%t(alpha[,notOther]) )
    
    if(TRAITS){
      Atrait <- bg%*%t(traitMat[,colnames(yp)])
      Strait <- traitMat[,colnames(yp)]%*%sg%*%t(traitMat[,colnames(yp)])
      agibbs[g,] <- Atrait
      mgibbs[g,] <- Strait
    }
    if(length(predXcols) > 0){
      if(!is.null(isNonLinX)){
        xpred <- .predictY2X_nonLinear(xpred,yy=.w2z(w[,notOther],sg[notOther,notOther]),
                                      bb=alpha[,notOther],ss=.cov2Cor(sg[notOther,notOther]),
                                      priorIV = priorXIV,priorX=priorX,
                                      predCols=isNonLinX,isInt,intMat,
                                      isFactor,factorList)$x
      }
      xpred[,px] <- .predictY2X_linear(xpred,yy=.w2z(w[,notOther],sg[notOther,notOther]),
                                      bb=alpha[,notOther],ss=.cov2Cor(sg[notOther,notOther]), 
                                      priorIV = priorXIV, 
                                      priorX=priorX,predCols=px)[,px]
    }
    
    if(g > burnin){

      ntot   <- ntot + 1
      ypred  <- ypred + yp
      ypred2 <- ypred2 + yp^2
      sumDev <- sumDev - 2*sum(.dMVN(w[,notOther],x%*%alpha[,notOther],
                                sg[notOther,notOther]) )
      sMean  <- sMean + sg
      
      wpred  <- wpred + w
      wpred2 <- wpred2 + w^2
      
      yy <- yp[,notOther]
      yy[yy > 0] <- 1
      yy[yy < 0] <- 0
      richness <- .add2matrix(rowSums(yy),richness)
      
      if(mmiss > 0){
        ymissPred[ymiss]  <- ymissPred[ymiss] + y[ymiss]
        ymissPred2[ymiss] <- ymissPred2[ymiss] + y[ymiss]^2
      }
      
      if(length(predXcols) > 0){
        predx  <- predx + xpred
        predx2 <- predx2 + xpred^2
      }
        
      if(TRAITS){
        yw     <- sweep(yp,1,rowSums(yp),'/')
        Ttrait <- .gjamPredictTraits(yw,traitMat[,colnames(yp)], traitTypes)
        tpred  <- tpred + Ttrait
        tpred2 <- tpred2 + Ttrait^2
      }
    }
  }
  
  richness <- richness/ntot
  
  if(mmiss > 0){
    ymissPred[ymiss]  <- ymissPred[ymiss]/ntot
    ymissPred2[ymiss] <- sqrt(ymissPred2[ymiss]/ntot - ymissPred[ymiss]^2)
  }
  
  x[xmiss] <- missX/ng
  xmissMu  <- missX/ng
  xmissSd  <- sqrt( missX2/ng - xmissMu^2 )
  
  sMean <- sMean/ntot
  
  chains <- list( rgibbs = rgibbs, sensGibbs = sensGibbs,
                  sgibbs = sgibbs, bgibbs = bgibbs) 
  
  xpredMu <- predx/ntot
  xpredSd <- sqrt(predx2/ntot - xpredMu^2)
  
  tmp <- .processPars(bgibbs)$summary
  bMu <- matrix(tmp[,'estimate'],q,S)
  bSe <- matrix(tmp[,'se'],q,S)
  aMu <- aSe <- bMu
  wz  <- which(tmp[,'0.025'] > 0 | tmp[,'0.975'] < 0)
  bwz <- matrix( unlist(strsplit(names(wz),'_')),ncol=2,byrow=T)[,1]
  
  tmp <- .processPars(rgibbs)$summary
  rMu <- matrix(tmp[,'estimate'],S,S)
  rSe <- matrix(tmp[,'se'],S,S)
  wz  <- which(tmp[,'0.025'] > 0 | tmp[,'0.975'] < 0)
  rwz <- matrix( unlist(strsplit(names(wz),'_')),ncol=2,byrow=T)[,1]
  
  tmp <- .processPars(sgibbs)$summary
  sMu <- matrix(tmp[,'estimate'],S,S)
  sSe <- matrix(tmp[,'se'],S,S)
  
  yallZero <- which(!snames %in% bwz & !snames %in% rwz)  #no significant coefficients
  
  yMu <- ypred/ntot
  ySd <- sqrt(ypred2/ntot - yMu^2)
  cMu <- cuts
  cSe <- numeric(0)
  
  wMu <- wpred/ntot
  wpp <- pmax(0,wpred2/ntot - wMu^2)
  wSd <- sqrt(wpp)
  
  tMu <- tSd <- tMuOrd <- btMu <- btSe <- stMu <- stSe <- numeric(0)
  
  if(TRAITS){
    
    tMu <- tpred/ntot
    tSd <- sqrt(tpred2/ntot - tMu^2)
    wo  <- which(traitTypes == 'OC')    #predict ordinal scores
    M   <- ncol(tMu)
    
    if(length(wo) > 0){
      tMuOrd <- tMu*0
      for(j in wo)tMuOrd[,j] <- findInterval(tMu[,j],traitPart[j,]) - 1
      tMuOrd <- tMuOrd[,wo]
    }
    
    tmp <- .processPars(agibbs)$summary
    btMu <- matrix(tmp[,'estimate'],q,M)
    btSe <- matrix(tmp[,'se'],q,M)

    tmp <- .processPars(mgibbs)$summary
    stMu <- matrix(tmp[,'estimate'],M,M)
    stSe <- matrix(tmp[,'se'],M,M)
    
    rownames(btMu) <- rownames(btSe) <- colnames(x)
    colnames(btMu) <- colnames(btSe) <- rownames(stMu) <- colnames(stMu) <- 
      rownames(stSe) <- colnames(stSe) <- tnames
      
    chains <- append( chains,list('agibbs' = agibbs))
    chains <- append( chains, list('mgibbs' = mgibbs) ) 
  }
  
  # note: on latent w scale
  meanDev <- sumDev/ntot
  pd  <- meanDev + 2*sum(.dMVN(yMu[,notOther],x%*%bMu[,notOther],
                                         sMean[notOther,notOther]) )
  DIC <- 2*pd + meanDev
  
  score <- mean( .getScoreNorm(y[,notOther],yMu[,notOther],ySd[,notOther]^2),na.rm=T )  # gaussian w
  
  if('OC' %in% typeNames){
    nk  <- length(ordCols)
    tmp <- .processPars(cgibbs)$summary
    cMu <- matrix(tmp[,'estimate'],nk,ncut-3)
    cSe <- matrix(tmp[,'se'],nk,ncut-3)
    colnames(cMu) <- colnames(cSe) <- cnames[-c(1,2,ncut)]
    tmp <- .processPars(dgibbs)$summary
    aMu <- matrix(tmp[,'estimate'],q,S)
    aSe <- matrix(tmp[,'se'],q,S)
    colnames(aMu) <- snames
    
    rownames(cMu) <- rownames(cSe) <- snames[ordCols]
    chains <- c(chains,list(cgibbs = cgibbs, dgibbs = dgibbs))
  }
  
  if('PA' %in% typeNames){
    zMu <- yMu
    zSd <- ySd
  }
  
  prAbs <- betaZeroMu <- betaZeroSe <- numeric(0)
  
  if(Q == 2)xscore <- mean( .getScoreNorm(x[,2],xpredMu[,2],xpredSd[,2]^2) )
  if(Q > 2) xscore <- colMeans( .getScoreNorm(x[,-1],xpredMu[,-1],xpredSd[,-1]^2) )
  
  colnames(bMu)   <- colnames(bSe) <-
    rownames(rMu) <- colnames(rMu) <- rownames(rSe) <- 
    colnames(rSe) <- snames
  rownames(bMu)   <- rownames(bSe) <- xnames
  
  modelSummary <- append(modelSummary,
                         list(typeNames = typeNames,DIC = DIC, score = score, xscore = xscore,
                              tMuOrd = tMuOrd,tMu = tMu,tSd = tSd,cutMu = cMu, cutSe = cSe, 
                              betaMu = bMu, betaSe = bSe,corMu = rMu, corSe = rSe, 
                              sigMu = sMu, sigSe = sSe,
                              betaTraitMu = btMu, betaTraitSe = btSe, 
                              betaSigMu = stMu, betaSigSe = stSe,
                              xpredMu = xpredMu,xpredSd = xpredSd,
                              yMu = yMu, ySd = ySd, wMu = wMu, wSd = wSd, prAbs = prAbs))
  
  list(burnin=burnin,missingIndex = xmiss, missingX = xmissMu, missingXSd = xmissSd,yallZero = yallZero,
       chains = chains, x = x, y = y, holdoutIndex = holdoutIndex, richness = richness,
       yMissMu = ymissPred, yMissSd = ymissPred2, ymiss = ymiss, modelSummary = modelSummary,
       censor = censor, TRAITS = TRAITS, traitList = traitList)
}
