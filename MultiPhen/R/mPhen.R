mPhen.options<-function(type=c("regression","plot","geno.input","pheno.input","meta.analysis","misc"), descr = FALSE){
  regoptsall = list()
  descriptionall = list()
  if(type[1] == "all"){
    type = c("regression","plot","geno.input","pheno.input","meta.analysis","misc")
  }
  for(k in 1:length(type)){
    if(length(grep("regr",type[k]))>0){
      opts =list(
                 "mPhen.inverseRegress" ,getOption("mPhen.inverseRegress",TRUE),"whether to use genotypes as outcome",
                 "mPhen.JointModel" , getOption("mPhen.JointModel",TRUE), "whether to use a multivariate model", 
                 "mPhen.scoreTest" , getOption("mPhen.scoreTest",FALSE), "whether to use a score test rather than likelihood test, for ordinal inverse regression only",
                 "mPhen.calcHWE",getOption("mPhen.calcHWE",FALSE), "whether to calculate HWE values",
                 "mPhen.exactMethod",getOption("mPhen.exactMethod",NULL), "if non-null, then which exact method to use. Can be wald,midp,fisher,small, see ?oddsratio()",
                 "mPhen.rescale" , getOption("mPhen.rescale",100), "when using imputed values and ordinal regression, how to scale probabilities to get integer values",
                 "mPhen.variable.selection",getOption("mPhen.variable.selection",FALSE), "whether backward selection should be used to reduce number of associated variables",
                 "mPhen.link.geno" , getOption("mPhen.link.geno",  "ordinal" ), "the link function when using genotype as outcome, can be ordinal or gaussian",
                 "mPhen.adjustSinglePv" , getOption("mPhen.adjustSinglePv",FALSE),
                 "whether to use a Nyoldt-Sidak correction on single-trait pvalue associations based on effective number of phenotypes tested",
                 "mPhen.useExpected" , getOption("mPhen.useExpected",FALSE), "whether to converted imputed data to dosage data. Should set link.geno to gaussian in this case",
                 "mPhen.link.pheno" , getOption("mPhen.link.pheno",  NULL ), "link function when using phenotpe as outcome.  If it is NULL, then automatically determined based on number of unique levels in each phenotype",
                 "mPhen.ccarep" , getOption("mPhen.ccarep",5), "Number of repetitions when using mPhen.cca",
                 "mPhen.thresh.diff" , getOption("mPhen.thresh.diff", 1e-3), "Threshold on euclidian distances (L2) between genotype weights from one iteration to the next, and also between phenotype weights from one iteration and the next, when using mPhen.CCA",  
                 "mPhen.orthPheno" , getOption("mPhen.orthPheno","none"),"Whether to orthogonalise phenotypes prior to association. Can be 'none', 'PCA', or 'GramSchmidt'"
                 )

    } else if(type[k]=="plot"){
      opts = list(
                  "mPhen.sigThresh",getOption("mPhen.sigThresh",1.0),"only plot pvalues more significant than this",
                  "mPhen.noPhenosPerPlot",getOption("mPhen.noPhenosPerPlot",20), "maximum number of phenotypes to include per plot",
                  "mPhen.limitP" ,getOption("mPhen.limitP",1e-200), "the smallest pvalue to plot in manhattan (smaller values will be set to this value)", "mPhen.title" , "all", "base title for all plots",
                  "mPhen.minlogpv",getOption("mPhen.minLogPv",20),"the smallest pvalue to show on qqplot",
                  "mPhen.colourByChroms" , getOption("mPhen.colourByChroms",FALSE),"whether to colour different chroms on manhattan by colour, or to colour different phenotypes by different colours.  Phenotypes are already assigned different symbols",
                  "mPhen.cex",getOption("mPhen.cex",0.5),"scaling for x-axis in plots",
                  "mPhen.pool",getOption("mPhen.pool",FALSE),"in drawing qq plot, whether to pool pvalues across phenotypes before drawing plot",
                  "mPhen.Colv",getOption("mPhen.Colv", TRUE),"whether to use a dendrogram in heatmap on phenotypes.  In some cases making the dendrogram will throw an error, in which case set this to FALSE",
                  "mPhen.indexMatch",getOption("mPhen.indexMatch",NULL),"this is the co-ordinates of a locus you wish to highlight on heatmap",
                  "mPhen.onlyShowJoint",getOption("mPhen.onlyShowJoint",NULL),"whether to only show Joint model pvalues whenever a joint model is used.  If false, then joint model not plotted, if TRUE, only joint model plotted, and if NULL, everything is plotted",
                  "mPhen.nacolor",getOption("mPhen.nacolor","black"),"what colour to use on heatmap for NA"
                  )
    }else if(type[k]=="geno.input"){
      opts = list(

                  "mPhen.batch",getOption("mPhen.batch",100000),"the number of snp genotypes to read into memory",
                  "mPhen.baseValue",getOption("mPhen.baseValue",0), "Default genotype, use 0 for SNP genotypes, and 2 for CN genotypes",
                  "mPhen.format",getOption("mPhen.format","GT"), "Which format field to read in from VCF file",
                  "mPhen.starting" , getOption("mPhen.starting",0), "Defines start position on chromosome for analysis",
                  "mPhen.ending" , getOption("mPhen.ending",Inf), "Defines end position on chromosome for analysis",
                  "mPhen.thin" , getOption("mPhen.thin",1),  "Allows taking every 10th snp, default is 1",
                  "mPhen.rsidToDo",getOption("mPhen.rsidToDo",NULL), "Which rs identifiers to run the analysis on",
                  "mPhen.matchPosition",getOption("mPhen.matchPosition",FALSE),"whether to match mPhen.rsidToDo to position rather than id" ,
                  "mPhen.numGenoPCs",getOption("mPhen.numGenoPCs",0), "How many geno pcs to calculate.  Note: the relatedness matrix is updated if this value is greater than zero, but the PCs (eigenvectors of this matrix are not calculated until closeConnection = TRUE",
                  "mPhen.onlyCommon", getOption("mPhen.onlyCommon" , FALSE),"whether to restrict to only positions in common when merging multiple cohorts"        
                  )
    }else if(type[k]=="pheno.input"){
      opts = list(
                  "mPhen.sep.pheno" ,getOption("mPhen.sep.pheno","\t"), "separator between cols in phenotype file",
                  "mPhen.onlyCommonPheno" ,getOption("mPhen.onlyCommonPheno",TRUE), "whether to include only common phenotypes",
                  "mPhen.numHeaderRows.pheno" ,getOption("mPhen.numHeaderRows.pheno",1), "how many header rows in phenotype file",
                  "mPhen.naRowThresh" ,getOption("mPhen.naRowThresh",0.1), "exclude samples with more than this fraction of missing genotypes",
                  "mPhen.naColThresh" , getOption("mPhen.naColThresh",0.1), "exclude phenotypes with more than this fraction of missing values",
                  "mPhen.fillMissingPhens" , getOption("mPhen.fillMissingPhens",FALSE),"whether to impute missing phenotype values when reading in a phenotype file",
                  "mPhen.quantileThresh",getOption("mPhen.quantileThresh",0.9),"If  excludeFile has a second column, indicating quality, then include the top percentage only in analysis"
                  )
    }else if(type[k]=="misc"){
      opts = list(
                  "mPhen.log10p",getOption("mPhen.log10p",FALSE),"whether to calculate probabilities in log10 space",
                  "mPhen.defaultBackwardFrac",getOption("mPhen.defaultBackwardFrac" ,0.1 ),"what proportion of variables to remove at each iteration",
                  "mPhen.defaultBackwardThresh", getOption("mPhen.defaultBackwardThresh" , 0.01), "threshold for backward selection significance",
                  "mPhen.defaultUseBF",getOption("mPhen.defaultUseBF", FALSE),"whether to use bayes factor for backward selection",
                  "mPhen.resultsDir",getOption("mPhen.resultsDir", "results/"),"output directory",
                  "mPhen.metaMinCohortNum",getOption("mPhen.metaMinCohortNum" , 2),"the minimum number of cohorts which have a beta/se in order to report a meta-analysis result",

                  "mPhen.useFisherIfAllBetaNA",getOption("mPhen.useFisherIfAllBetaNA",TRUE),"if beta is not available, whether to revert to Fishers method for combining pvalues"
                  )
    }
    seqs = seq(1,length(opts),3)
    regopts = opts[seqs+1]
    descr1 = opts[seqs+2]
    names(regopts) = unlist(opts[seqs])
    names(descr1) = names(regopts)
    regoptsall = c(regoptsall,regopts)
    descriptionall = c(descriptionall,descr1)
  }
  if(descr) return(descriptionall)
  else return(regoptsall)
}




mPhen <-function(genoData, phenoData, 
                 phenotypes ="all",
                 covariates = NULL, resids = NULL, strats = NULL,
                 opts = mPhen.options(c("regression","pheno.input"))){
  if(opts$mPhen.fillMissingPhens) phenoData = .fillMissing(phenoData)
  limit = list(phenotypes = "all",covariates = covariates, resids = resids, strats = strats)
  phenoObject = mPhen.preparePheno(list(pheno = phenoData,limit = limit),  indiv = dimnames(genoData)[[1]],  opts = opts)
  results = mPhen.assoc(genoData,phenoObject,opts = opts)
  results  
}

mPhen.assoc <-function(genoData, 
                       phenoObject,  
                       opts = mPhen.options("regression"), subinds = 1:(dim(genoData)[1])
                       ){
  if(is.null(genoData)) return (NULL)
  testgeno = !is.array(genoData)
  if(!is.null(opts$mPhen.exactMethod))stop('ERROR, need to re-implement exact method')
  if(testgeno) {
    stop('ERROR! The genetic data is NOT in matrix format. Please reformat the genetic data as a matrix and rerun the test (see ?as.matrix)')
  }
  if(opts$mPhen.scoreTest & opts$mPhen.variable.selection) stop('score test not implemented for variable selection')
  ordinal = opts$mPhen.link.geno=="ordinal"
  isImputed = length(dim(genoData))==3
  if(isImputed && opts$mPhen.useExpected) {
    genoData <- .calcExpectedGeno(genoData)
    isImputed <- FALSE
    if(ordinal & opts$mPhen.inverseRegress) stop('ordinal not consistent with useExpected=TRUE for inverse regression')
  }

  if(length(which(dimnames(phenoObject$phendata)[[3]] !=dimnames(genoData)[[1]]))>0) stop('sample id names do not match between geno and pheno')
  if(is.null(dimnames(genoData)[[1]])){
    dimnames(genoData)[[1]] = 1:(dim(genoData)[1])
  }
  dimg = dim(genoData)
  if(opts$mPhen.JointModel & !opts$mPhen.inverseRegress & isImputed){
    warning('cannot run multi-geno regression on imputed data, converting to dosage')
    genoData = .calcExpectedGeno(genoData)
    isImputed = FALSE
    dimg = dim(genoData)
  }
  if(is.null(dimnames(genoData)[[2]])) dimnames(genoData)[[2]] = 1:(dim(genoData)[[2]])

  if(isImputed){
    if(opts$mPhen.scoreTest) stop('score test does not work with imputed data')
    phenoObject = .expandDat(phenoObject,dimg[3])
    subinds1 =  .subSample(subinds,dimg[3])
  }else{
    subinds1 = subinds
    rescale = 1
  }
  #if(opts$mPhen.scoreTest & !ordinal) stop('score test only coded for ordinal regression')
  x = .makeGenoWeights(genoData,rescale=opts$mPhen.rescale,ordinal=ordinal)
  res = .mPhen(x, phenoObject,subinds =subinds1,  opts=opts)

  if(phenoObject$effPhe>1 & opts$mPhen.adjustSinglePv) {
    for(k in 1:dim(res$Res)[1]){
      for(j in 1:dim(res$Res)[2]){
        phenNames = dimnames(res$Res)[[3]]
        for(l in 1:length(phenNames)){
          if(phenNames[l]!="JointModel"){
            res$Res[k,j,l,2] <-  .applySidakCorrection(res$Res[k,j,l,2],phenoObject$effPhe)
          }
        }
      }
    }
  }

  attr(res,"opts")<-opts
  res
}






## this calculates best linear combination of snps and phenotypes
mPhen.cca<-function(genoData, phenoObject,
                    opts = mPhen.options("regression"), 
                    subinds = 1:(dim(genoData)[1]),
                    vs.G =opts$mPhen.variable.selection,
                    vs.P = opts$mPhen.variable.selection
                    ){
  if(is.null(genoData)) return (NULL)
  isImputed = length(dim(genoData))==3
  if(isImputed) genoData = .calcExpectedGeno(genoData)
  opts$mPhen.scoreTest = FALSE; opts$mPhen.calcHWE=FALSE; 
  noPhenos = length(phenoObject$phenN)
  if(noPhenos<=1) stop('need more than one pheno, otherwise just use mPhen.assoc')
  # if(thresh_initial < 1.0){
  #      opts$mPhen.inverseRegress = TRUE; opts$mPhen.JointModel = TRUE;       opts$mPhen.variable.selection = vs.P
  #      resultsJoint = mPhen.assoc(genoData, phenoObject,  opts = opts,  subinds = subinds)
  #      pvjoint = resultsJoint$Res[,,(noPhenos+1),2]
  #      pvinds = which(pvjoint < opts$mPhen.thresh.initial)
  #  }else{
  pvinds = 1:dim(genoData)[[2]]
  #  }

  geno1 = genoData[,pvinds,drop=F]   
  noSnps = dim(geno1)[2]
  betasg = rep(1,noSnps)
  phend = phenoObject$phendata[,1,]
  phend2 = phenoObject$phendata[1,,,drop=F]
  phenoObject2 = list(phendata = phend2, index = phenoObject$index[1],
                      index_strat = phenoObject$index_strat, index_resid = phenoObject$index_resid,
                      pheno_cov1 = phenoObject$pheno_cov1, pheno_resid1= phenoObject$pheno_resid1, stratMatrix1 =phenoObject$stratMatrix1, 
                      families = "gaussian", phenN="combined", effPhe = 1)
  betasp_prev = rep(1,noPhenos)
  betasg_prev = betasg
  delta = 1
  deltag = 1
  cntk =0
  while(delta > opts$mPhen.thresh.diff | deltag > opts$mPhen.thresh.diff){
    opts$mPhen.inverseRegress = TRUE; opts$mPhen.JointModel = TRUE;       opts$mPhen.variable.selection = vs.P
    geno2 = geno1 %*% betasg
    dimnames(geno2)[[2]] = "-1_-1"
    resultsJoint = mPhen.assoc(geno2, phenoObject,  opts = opts,  subinds = subinds)
    betasp = resultsJoint$Res[,,1:noPhenos,1]
    nap = is.na(betasp)
    normp = sqrt(betasp[!nap] %*% betasp[!nap])
    if(normp==0){
      opts$mPhen.variable.selection = FALSE

      resultsJoint = mPhen.assoc(geno2, phenoObject,  opts = opts, subinds = subinds)
      betasp = resultsJoint$Res[,,1:noPhenos,1]
      nap = is.na(betasp)
      normp = sqrt(betasp[!nap] %*% betasp[!nap])
    }
    betasp = betasp / normp
    phend2[,1,] = betasp %*% phend
    phenoObject2[[1]] <- phend2
    opts$mPhen.inverseRegress = FALSE; opts$mPhen.JointModel = TRUE;       opts$mPhen.variable.selection = vs.G
    resultsSingle1 = mPhen.assoc(geno1, phenoObject2,  opts = opts,subinds = subinds)
    betasg = resultsSingle1$Res[,1:noSnps,,1]
    nag = is.na(betasg)
    normg = sqrt(betasg[!nag] %*% betasg[!nag])
    if(normg==0){
      opts$mPhen.variable.selection = FALSE
      resultsSingle1 = mPhen.assoc(geno1, phenoObject2, opts = opts,  subinds = subinds)
      betasg = resultsSingle1$Res[,1:noSnps,,1]
      nag = is.na(betasg)
    }
    betasg = betasg / normg

    delta = (betasp_prev[!nap]-betasp[!nap]) %*% (betasp_prev[!nap]-betasp[!nap])
    deltag = (betasg_prev[!nag]-betasg[!nag]) %*% (betasg_prev[!nag]-betasg[!nag])
    print(paste("delta",delta, deltag)) 
    betasp_prev = betasp
    betasg_prev = betasg
    cntk = cntk+1
    if(cntk>opts$mPhen.ccarep) break
  }
  betasg1 =rep(0, dim(genoData)[2])
  betasg1[pvinds] = betasg
  list(betasp = betasp, betasg = betasg1, resultsGeno = resultsSingle1, resultsPheno = resultsJoint, gComb = geno1 %*% betasg, pComb =betasp %*% phend)
}


#

mPhen.preparePheno<-function(phenoData, 
                             pcs = NULL,        
                             indiv=if(is.null(pcs)) rownames(phenoData) else rownames(pcs),
                             opts = mPhen.options("regression")
                             ){
  fillMissingPhens = FALSE
  orthPheno = opts$mPhen.orthPheno=="PCA"  
  usePCA=orthPheno=="PCA"
  orthogonalise = orthPheno=="GramSchmidt"
  link.pheno = opts$mPhen.link.pheno
  if(is.list(phenoData)){
    limit = phenoData$limit
    phenoData = phenoData$pheno
  }else{
    limit = list(phenotypes = "all")  #,covariates = NULL, resids = NULL, strats = NULL, excls = NULL)  
  }
  if(is.null(indiv)) stop('this cannot be null')
  if(!is.null(pcs)){
    if(length(which(indiv!=dimnames(pcs)[[1]]))>0) stop("if pcs are specified, they should be in same order as indiv")
  }
  testpheno = !is.matrix(phenoData) || is.null(dimnames(phenoData)[[1]]) ||  is.null(dimnames(phenoData)[[2]])
  if(testpheno){
    stop('ERROR! the phenotype data is NOT in matrix format. Please reformat. 
         NOTE, if testing one single phenotype remember to add drop = F (see ?Extract).  
         Also, dimnames(phenoData) must be defined.  
         Try pheno = as.matrix(read.table(phenoFile,header=T,row.names=1,as.is=T,sep="\t",comment=""))')
  }
  mat = match(indiv,dimnames(phenoData)[[1]])
  if(length(mat[!is.na(mat)])==0) stop('probably matching pheno and geno names')
  pheno1 = phenoData[mat,,drop=FALSE]
  if(!is.null(pcs)){
    pheno1 = cbind(pheno1,pcs)
  }
  if(is.list(limit)){
    limit = .getLimitMatrixFromList(limit)
  }
  {
    if(is.null(dim(limit)) || dim(limit)[[2]]<2 ||  min(apply(limit,c(1,2),is.character))==0)
      stop('ERROR limit should be in matrix format with more than 2 cols and in string format.  
           Try read.table(limitfile,as.is=T,fill=T,header=F)')
    if(length(grep('#',limit[,1])>0)) limit = limit[-grep('#',limit[,1])>0,]
    ph_ind = grep("^pheno",limit[,1])
    if(length(ph_ind)==0 ){
      lmt1 = cbind(rep('pheno', dim(pheno1)[[2]]), dimnames(pheno1)[[2]])
      limit = rbind(limit,lmt1)
    }else if(length(ph_ind)==1 & limit[ph_ind[1],2]=="all"){
      #print(limit)
      toexcl_ = c(limit[grep("^covar",limit[,1]),2],
                  limit[grep("^resid",limit[,1]),2],
                  limit[grep("^resid",limit[,1]),2],
                  limit[grep("^nopheno",limit[,1]),2] )
      if(!is.null(pcs)) toexcl_ = c(toexcl_,dimnames(pcs)[[2]])
      phensI = dimnames(pheno1)[[2]]
      if(length(toexcl_)>0) phensI = phensI[!(phensI %in% toexcl_)]
      lmt1 = cbind(rep('pheno', length(phensI)),phensI)
      if(dim(limit)[2]>2) lmt1 = cbind(lmt1,rep(limit[ph_ind,3], length(phensI)))
      dimnames(lmt1)[[2]] = dimnames(limit)[[2]]
      limit = rbind(limit[-ph_ind[1],],lmt1)
    }
    splitind = grep(";",limit[ph_ind,2])
    if(length(splitind)>0){
      for(k in splitind){
        splitstr = strsplit(limit[ph_ind[k],2],";")[[1]]
        lmt1 = cbind(rep('pheno', length(splitstr)), splitstr)
        if(dim(limit)[2]>2) lmt1 = cbind(lmt1,rep(limit[k,3],length(splitstr)))
        dimnames(lmt1)=NULL
        limit = rbind(limit,lmt1)
      }
      limit = limit[-ph_ind[splitind],]
    }
    if(length(which(is.na(match(limit[,1],c("pheno","covar","strat","resid", "excl","nopheno")))))>0) 
      stop('should match pheno,covar,resid or strat')
    phenotypes = limit[grep("^pheno",limit[,1]),2]
  }
  transformed = .applyTrans(pheno1,limit)
  pheno1 = transformed$pheno1
  #if(usePCA) fillMissingPhens = TRUE
  if(fillMissingPhens  &   length(which( apply(is.na(pheno1),1,sum)>0))>0){
    #     pheno1 = kNNImpute(pheno1,3)$x
    pheno1 = .fillMissing(pheno1)
  }
  limit = transformed$limit
  pheno1 = .fixNonNumeric(pheno1)
  phenNames = dimnames(pheno1)[[2]]
  todo_ = limit
  todo = sort(todo_[grep("^pheno",todo_[,1]),2])
  #nottodo = sort(todo_[grep("^nopheno",todo_[,1]),2])
  excl = todo_[grep("^excl",todo_[,1]),2]
  covar = todo_[grep("^covar",todo_[,1]),2]
  stratify=todo_[grep("^strat",todo_[,1]),2]
  resid = todo_[grep("^resid",todo_[,1]),2]
  index_cov = NULL
  index = NULL
  index_resid = NULL
  index_strat = NULL
  exclindex = NULL
  for ( ik in todo) index =c(index,which(phenNames==ik))
  for ( ik in excl) exclindex =c(exclindex,grep(ik,phenNames))
  for ( ik in covar) index_cov = c(index_cov,which(phenNames==ik))
  for ( ik in resid) index_resid = c(index_resid,which(phenNames==ik))
  for ( ik in stratify) index_strat = c(index_strat,which(phenNames==ik))
  index = unique(index)
  exclindex = unique(exclindex)
  toexcl=c()
  for(i in exclindex){
    toexcl=unique(c(toexcl,which(pheno1[,i]==1 | is.na(pheno1[,i]))))
  }
  print(paste("excluding",length(toexcl),"samples based on exclusion criteria"))
  if(length(toexcl)>0) pheno1 =pheno1[-toexcl,]
  index_cov = unique(index_cov)
  index_resid = unique(index_resid)
  index_strat = unique(index_strat)
  toStrat = pheno1[,index_strat,drop=FALSE]
  stratNames = dimnames(pheno1)[[2]][index_strat]
  pheno_cov = pheno1[,index_cov,drop=FALSE]
  pheno_resid =pheno1[,index_resid,drop=FALSE]
  pheno1 = pheno1[,index,drop=FALSE]
  loadings = NULL
  if(usePCA){
    narow = apply(is.na(pheno1),1,sum)>0
    sv = svd(apply(pheno1[!narow,],2,.centralise))
    pheno2 = array(NA, dim = c(dim(pheno1)[1],dim(sv$u)[2]), dimnames =  list(dimnames(pheno1)[[1]],paste("PC",1:length(sv$d),sep="")))
    pheno2[!narow,] = sv$u
    loadings = solve(diag(sv$d) %*% t(sv$v))
    pheno1 = pheno2
  }
  else if(orthogonalise){
    narow = apply(is.na(pheno1),1,sum)>0
    orth = .orthogonalise(apply(pheno1[!narow,],2,.centralise))
    pheno2 = array(NA, dim = c(dim(pheno1)[1],dim(orth)[2]), dimnames =  list(dimnames(pheno1)[[1]],paste("ORTH",1:length(orth),sep="")))
    pheno2[!narow,] = orth$pheno 
    pheno1 = pheno2
    loadings = orth$W
  }
  #  print(max(mat,na.rm=T))
  #  print(dim(pheno_cov))
  pheno_cov1 = pheno_cov
  pheno_resid1 = pheno_resid
  stratMatrix1 = .makeFactor(toStrat,stratNames)
  if(dim(pheno_cov1)[2]==0) phenoCovna = rep(0,length(indiv)) else phenoCovna = is.na(apply(pheno_cov1,1,sum))
  if(dim(pheno_resid1)[2]==0) phenoResna = rep(0,length(indiv)) else phenoResna = is.na(apply(pheno_resid1,1,sum))
  families = rep(link.pheno, dim(pheno1)[2])
  inclM = as.matrix(apply(pheno1,2,.getIncl, phenoCovna | phenoResna))
  if(is.null(link.pheno)){
    for(i in 1:(dim(pheno1)[2])) families[i] = .getFamily(pheno1[,i],inclM[,i])
  }
  binom = families=="binomial"
  caseInd = apply(pheno1,2,.iscase)
  for(i in 1:(dim(pheno1)[2])) if(binom[i]) pheno1[,i] = pheno1[,i]-min(pheno1[,i],na.rm=TRUE)
  ##if(dim(pheno_cov1)[2]>0)for(i in 1:(dim(pheno_cov1)[2])) pheno_cov1[,i] = .standardise(pheno_cov1[,i])
  typenames = c("pheno","incl","caseInd","offsets")
  phenN = dimnames(pheno1)[[2]]
  arraynames = list(phenN, typenames,indiv )
  phendata = .getArray(arraynames)
  phendata[,1,] = t(pheno1)
  phendata[,2,] = t(inclM)
  phendata[,3,] = t(caseInd)
  if(length(index)>0) for(i in 1:(length(index))) phendata[i,4,] = .calcOffset(phendata[i,,],families[i], pheno_resid1)
  list(index=index,index_strat=index_strat,index_resid=index_resid,phendata= phendata,pheno_cov1=pheno_cov1,
       ### pheno_resid1=pheno_resid1,
       stratMatrix1=stratMatrix1,families=families,phenN = phenN, effPhe=.nyholdtSidak(pheno1), loadings = loadings,limit = limit)
}






## can be used to check required parameters are defined
##tocheck is a list of params which should be defined
##allnames can be obtained from  ls(all.names=T)
## if some parameters are not defined will throw an error
## This function will also 'fix' assignment. In particular, if there is a string with a $WORK, then
## the function will replace $WORK with the system defined variable value
mPhen.defineOptions<-function(file = NULL, getOptionsFromCommandLine = TRUE){
  if(!is.null(file)){
    if(file.exists(file))	 source(file)
    else warning(paste(file,'does not exist, not sourcing options'))
  }
  if(getOptionsFromCommandLine) .parseOpts()
  outDir=getOption("mPhen.resultsDir",NULL)
  tocheck = .Options[grep("mPhen",names(.Options))]
  outfile = NULL
  if(!is.null(outDir)) outfile = paste(outDir,"options.txt",sep="/")
  for(i in 1:length(tocheck)){
    if(!is.null(outfile)){
      write(names(tocheck)[i],file=outfile,append=i>1)
    }
    v = tocheck[[i]]
    if(is.list(v)){
      if(!is.null(outfile)){
        for(j in 1:length(v)){
          write(v[[j]],file=outfile,append=i>1)
        }
      }
      starind = grep('\\*',v)
      if(length(starind)>0){
        v2 = list()
        for(kk1 in 1:length(starind)){
          kk = starind[kk1]
          v1 <- v[[kk]]
          inds <- grep('\\*',v1)
          if(length(inds)>0){
            toadd <- c()
            for(jj in 1:length(inds)){
              line <- paste("ls",v1[inds[jj]])
              toadd <- c(toadd,system(line,intern=TRUE))
            }
            v2[[kk1]] <- c(v1[-inds],toadd)
          }
        }
        v[starind] <- v2
        tocheck[[i]] = v
      }
    }
    else{
      if(!is.null(outfile)){
        write(v,file=outfile,append=i>1)
      }
    }
    if(is.character(v)){
      for(kk in 1:length(v)){
        v1 = strsplit(v[kk],"/")[[1]]
        subinds = grep('\\$',v1)
        for(k in subinds){
          line = paste("echo",v1[k])
          v1[k] = system(line,intern=TRUE)
        }
        v[kk] = paste(v1,collapse="/",sep="")
      }
      v = sub('gb','000000000',v)
      v = sub('Gb','000000000',v)
      v = sub('mb','000000',v)
      v = sub('Mb','000000',v)
      v = sub('kb','000',v)
      v = sub('Kb','000',v)
      #print(paste("assign",allnames[index],v))
      #assign(allnames[index],v,inherits =TRUE)
      tocheck[[i]] = v
    }
  }
  options(tocheck)
}










