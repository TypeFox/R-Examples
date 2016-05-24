anova.manyany = function(object, ..., nBoot=99, p.uni="none", block = object1$block, nCores = 1, bootID=NULL, replace=TRUE)
{
  #set default na.action to exclude in order to not change dimensions of anything when NA's are present
  naOptInit = getOption("na.action")
  options(na.action="na.exclude")
  
  object1 = object 
  # get object 2
    dots <- list(...)
    ndots <- length(dots)
    fndObj2 <- FALSE
    if (ndots==0) {
       stop("missing a second manyany object")
    }else {
       if (ndots>1)
           warning("current version only compares two manyany objects")
       for (i in 1:ndots) {
           if (any(class(dots[[i]])=="manyany")){
              object2 <- dots[[i]]
              fndObj2 <- TRUE
              break
           }
       }
       if (!fndObj2) stop("cannot find object 2")
    }

  if(any(names(object1$call)=="composition"))
  {
    if(object1$call$composition==TRUE) #recode so that it fits compositional models as univariate, to save time and fuss/bother.
    {
      object1$call$formula = object1$formula
      object2$call$formula = object2$formula
      object1$call$data = object1$model
      object2$call$data = object2$model
      object1$residuals = as.matrix(c(object1$residuals))
      object1$call$composition=FALSE
      object2$call$composition=FALSE
      assign(as.character(object1$call[[3]]),object1$model$y) 
      assign(as.character(object2$call[[3]]),object2$model$y) 
    }
  }
  n.rows = dim(object1$resid)[1]
  n.vars = dim(object1$resid)[2]
  
  qfn = rep(NA,n.vars)
  for(i.var in 1:n.vars)
  {
    if(grepl("egative",object1$family[[i.var]]$family) || object1$family[[i.var]]$family == "negbinomial")
      qfn[i.var] = "qnbinom"
    if(object1$family[[i.var]]$family=="poisson")
      qfn[i.var] = "qpois"
    if(object1$family[[i.var]]$family=="binomial")
    {
      qfn[i.var] = "qbinom"
      warning("The binomial option of manyany currently assumes you have presence/absence data")
    } 
    if(object1$family[[i.var]]$family=="gaussian")
      qfn[i.var] = "qnorm"  
    if(object1$family[[i.var]]$family=="Tweedie")
      qfn[i.var] = "qtweedie"
    if(object1$family[[i.var]]$family=="ordinal")
      qfn[i.var] = "qordinal"
  }

  if(is.null(bootID)==FALSE)
  {
    bootID = as.matrix(bootID)
    if(dim(bootID)[2]!=n.rows)
      stop("Number of columns of bootID must match number of rows in data")
    nBoot = dim(bootID)[1] #overwriting nBoot with value implied by user-entered ID matrix
    block = NULL #overwriting previous value for block
    print("User-entered bootID matrix will be used to generate bootstrap samples")
  }  
  n.levels  = n.rows; unlistIDs = NULL #objects needed for block resampling otherwise ignorable 
  blockIDs = NULL
  if(is.null(block)==FALSE)
  {
    tb=table(block)
    n.levels = length(tb)
    if(any(tb!=n.rows/n.levels))
    {   
      print(tb) 
      stop("Sorry, block needs to be a balanced factor - same number of rows for each level")
    }
    else
    {
      blockIDs = vector("list",n.levels)
      for(i.level in 1:n.levels)
        blockIDs[[i.level]] = which(block==names(tb)[i.level])
      unlistIDs = unlist(blockIDs) #needed to match each resampled observation with its correct location
    }
  }
  #get observed test stat
#  ft.1i=eval(object1$call) #this call seems unnecessary
#  ft.2i=eval(object2$call) #this call seems unnecessary
  statj = 2 * ( logLik(object2)-logLik(object1) )
  stat = sum(statj)
  

  if(nCores>1)
  {
    nBooti = ceiling(nBoot/nCores)

    # construct a list which says which rows of bootID to use in which cluster: only needed when bootID provided
    bootRows = vector(length=nCores,mode="list")
    for(iCore in 1:nCores)
      bootRows[[iCore]] = 1:nBooti + nBooti*(iCore-1)
    bootRows[[nCores]]  = pmin( bootRows[[nCores]], nBoot )

    # set up clusters, pass through arguments
    cl=makeCluster(nCores)
    argList = list(bootID=bootID, block=block, blockIDs = blockIDs, n.rows=n.rows, n.vars=n.vars, replace=replace, unlistIDs=unlistIDs, n.levels=n.levels, object1=object1, object2=object2, qfn=qfn)
    clusterExport(cl,"argList", envir=environment())
    #clusterExport(cl,c("nBooti","bootID","block","n.rows","n.vars","replace","unlistIDs","n.levels","object1","object2","qfn"), envir=environment())
    out=parLapply(cl,bootRows,bootAnova)
    #why not clusterapply??
    stopCluster(cl)

    # store results in vectors/matrices not lists
    stat.i  = rep(NA,nBooti*nCores)
    statj.i = matrix(NA,n.vars,nBooti*nCores)
    for(i.core in 1:nCores)
    {
      stat.i[(i.core-1)*nBooti+1:nBooti]     = out[[i.core]]$stati.i
      statj.i[,(i.core-1)*nBooti+1:nBooti] = out[[i.core]]$statj.ii
    }  
    stat.i = stat.i[1:nBoot]
    statj.i = statj.i[,1:nBoot]
  }
  else
  {
    out = bootAnova(bootRows=1:nBoot,bootID=bootID,block=block, blockIDs=blockIDs, n.rows=n.rows,n.vars=n.vars,replace=replace,unlistIDs=unlistIDs,n.levels=n.levels,object1=object1,object2=object2,qfn=qfn,nCores=1)
    stat.i=out$stati.i
    statj.i = out$statj.ii
  }
  if(n.vars>1)
    dimnames(statj.i)[[1]] = dimnames(object1$residuals)[[2]]

  p = ( 1 + sum(stat.i>stat-1.e-8) ) / (nBoot + 1)
  pj = ( 1 + apply(statj.i>statj-1.e-8,1,sum) ) / ( nBoot + 1)

  class(stat.i) = "numeric"
  if(p.uni=="unadjusted")
    result = list(stat=stat,p=p,uni.test=statj,uni.p=pj,stat.i=stat.i,statj.i=statj.i,p.uni=p.uni,nBoot=nBoot) 
  if(p.uni=="none")
    result = list(stat=stat,p=p,stat.i=stat.i,p.uni=p.uni,nBoot=nBoot) 
  options(na.action=naOptInit) #restore previous default for na.action

  class(result) = "anova.manyany"
  return(result)  
}

print.anova.manyany=function(x, ...)
{
  #get overall results in a table
  table=matrix(c(x$stat,x$p),1,2)
  dimnames(table)[[2]]=c("LR","Pr(>LR)")
  dimnames(table)[[1]]=c("sum-of-LR")

  allargs <- match.call(expand.dots = FALSE)
  dots <- allargs$...
  s.legend = TRUE
  if(length(dots)>1)
  {
    if("signif.legend" %in% dots)
      s.legend = signif.legend
  }
  if(x$p.uni=="none")
    signif.legend = s.legend
  else
    signif.legend = FALSE
  
  #print overall results
  cat("\n")
  printCoefmat(table, tst.ind=1, P.values=TRUE, has.Pvalue=TRUE, signif.legend=signif.legend, eps.Pvalue=1/(x$nBoot+1-1.e-8),...)
  cat("\n")
  #print univariate results in a table, if required
  if(x$p.uni!="none")
  {
    signif.legend = s.legend
    tablej=cbind(x$uni.test,x$uni.p)
    dimnames(tablej)[[2]]=c("LR","P(>LR)")
    printCoefmat(tablej, tst.ind=1, P.values=TRUE, has.Pvalue=TRUE, signif.legend=signif.legend, eps.Pvalue=1/(x$nBoot+1-1.e-8), ...)
  }
}

bootAnova = function(bootRows,...)
{
  nBooti = length(bootRows)
  dots = list(...)
  if ( any(names(dots)=="nCores") ) # if nCores=1, take ... and call it argList, to match parLapply call
  {
    if(dots$nCores==1)
      argList=list(...)
  }
  #initialise parameters for bootstrapping

  require(mvabund)
  yMat = matrix(NA,argList$n.rows,argList$n.vars)
  if(argList$object1$family[[1]]$family=="ordinal")
    yMat=data.frame(yMat)
  argList$object1$call$get.what="none" #to avoid wasting time computing residuals etc when resampling
  argList$object2$call$get.what="none" #to avoid wasting time computing residuals etc when resampling

  stati.i  = rep(NA,nBooti)
  statj.ii = matrix(NA,argList$n.vars,nBooti)
  
  if(is.null(argList$bootID))
    boot.Resamp = rep(NA,argList$n.rows)

  # need to find data frame and call it what it the same as in the original call for analysis
  whichData = which(names(argList$object2$call)=="data")
  assign(as.character(argList$object2$call[[whichData]]),argList$object2$model) 
    
  #now do the bootstrap
  for(iBoot in 1:length(bootRows))
  {
    # first get ID entries for this resample
    if(is.null(argList$bootID)==FALSE)
      boot.Resamp = argList$bootID[bootRows[iBoot],]
    else      
    {
      if(is.null(argList$block))
        boot.Resamp = sample(1:argList$n.rows,replace=argList$replace)
      else
        boot.Resamp[argList$unlistIDs] = unlist(argList$blockIDs[sample(argList$n.levels,replace=argList$replace)]) #unlistIDs is needed to make sure each unlisted blockID ends up in the right place
    }

    # resample PIT residuals
    resid.i = as.matrix(argList$object1$residuals[boot.Resamp,])
  
    # now use PIT-transform to get resampled yMat
    for(i.var in 1:argList$n.vars)
    {
      qparams = argList$object1$params[[i.var]]
      qparams[[1]]=resid.i[,i.var]
      names(qparams)[1]="p"
      yMat[,i.var] = do.call(argList$qfn[i.var], qparams)
    }
    
    #save resampled yMat as whatever the original yMat was called in workspace - but without zerotons
    if(argList$object1$family[[1]]$family=="ordinal")
      is.zeroton = apply(yMat,2,function(x) length(table(x)))==1
    else
      is.zeroton = apply(yMat,2,sum,na.rm=TRUE)==0
    assign(as.character(argList$object1$call[[3]]),yMat[,is.zeroton==FALSE]) 
    assign(as.character(argList$object2$call[[3]]),yMat[,is.zeroton==FALSE]) 

    #re-fit manyany functions and calculate test stats using the resampled yMat:
    if(sum(is.zeroton==FALSE)>0)
    {
#      recover()  
      ft.1i=eval(argList$object1$call)
      ft.2i=eval(argList$object2$call)
      statj.ii[is.zeroton==FALSE,iBoot]=2 * ( logLik(ft.2i)-logLik(ft.1i) )
      stati.i[iBoot] = sum(statj.ii[,iBoot], na.rm=TRUE)      
    }
    else
      stati.i[iBoot] = 0
  }
  return(list(stati.i=stati.i,statj.ii=statj.ii))
} 
