#function for message on load
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("HHG package for non parametric tests of independence and equality of distributions.")
  packageStartupMessage("type vignette(\'HHG\') or ?HHG for documentation, examples and a quickstart guide.")
  packageStartupMessage("use suppressPackageStartupMessages(library(HHG)) to suppress this message.")
}

# function for computing the statistic of the univariate df test of independence
hhg.univariate.ind.stat = function(x, y, variant = 'ADP',aggregation.type='sum',score.type='LikelihoodRatio', mmax = max(floor(sqrt(length(x))/2),2), mmin =2, w.sum = 0, w.max = 2){
  correct.mi.bias = F # currently mi bias correction is not supported in interface
  type='Independence'
  flag_perform_ADP_multiple_partitions = F #flag used to check whether all partition sizes can be computed at single call
  
  #input checks:
  if(is.null(x) |is.null(y)){
    stop("x & y should be vectors of doubles." )
  }
  if(!(correct.mi.bias==TRUE || correct.mi.bias==FALSE) || w.sum!=as.integer(w.sum) || w.max!=as.integer(w.max) ){
    stop("Correct bias should be boolean, w.sum & w.max should be integers.")
  }
  if(length(x)!= length(y)){
    stop("X & Y not the same length")
  }
  if(!is.numeric(mmin) | !is.numeric(mmax)){
    stop(" parameters should be 2<= mmin<= mmax, integers")  
  }
  if(mmin !=round(mmin) | mmax !=round(mmax)){
    stop(" parameters should be 2<= mmin<= mmax, integers")  
  }
  if(mmin<2 | mmax<mmin){
    stop(" parameters should be 2<= mmin<= mmax, integers")  
  }
  
  .hhg.univariate.check.inputs(type,variant,length(x),score.type,aggregation.type)
  
  max_not_available=FALSE
  if (((variant == 'DDP') && ((mmax > 4))) || ((variant == 'ADP') && ((mmax > 2)))) {
    max_not_available=TRUE
  }
  if(max_not_available == TRUE  & is.element(aggregation.type,c('max'))){
    stop(" Maximum based scores only available for ADP with m=2 and DDP with m=2,3,4")
  }
  if(max_not_available == TRUE  & is.element(aggregation.type,c('both'))){
    warning(" Maximum based scores only available for ADP with m=2 and DDP with m=2,3,4")
  }
  if(variant == 'ADP' & correct.mi.bias == F & aggregation.type == 'sum'){
    flag_perform_ADP_multiple_partitions = T
  }
  
  #statistics computation:
  stat.sl = rep(NA,(mmax)-mmin+1)
  stat.sc = rep(NA,(mmax)-mmin+1)
  stat.ml = rep(NA,(mmax)-mmin+1)
  stat.mc = rep(NA,(mmax)-mmin+1)
  
  if(flag_perform_ADP_multiple_partitions){ 
    res = .hhg.test.adp_mk(x,y,M=(mmin:mmax))
    stat.sl = res$sum.lr
    stat.sc = res$sum.chisq
  }else{
    ms_vector=(mmin:mmax)
   for(i in 1:length(ms_vector)){
     m=ms_vector[i]
     #call by different variant an m to specific functions (some of which can compute max variants or more efficiently for small m)
     if (variant == 'DDP') {
       if (m != as.integer(m) || m < 2) {
         stop('m, number of partitions,  must be an integer greater than 1')
       } else if (m == 2) {
         x.variant = 'spr.obs'
       } else if (m == 3) {
         x.variant = 'ppr.33.obs'
       } else if (m == 4) {
         x.variant = 'tpr.obs'
       } else {
         x.variant = 'ddp.obs'
       }
     } else if (variant == 'ADP') {
       if (m != as.integer(m) || m < 2) {
         stop('m, number of partitions,  must be an integer greater than 1')
       } else if (m == 2) {
         x.variant = 'spr.all'
       } else if (m == 3) {
         # One could use 'ppr.33.all' that has the same complexity, but in practice it is much slower
         x.variant = 'ddp.all'
       } else if (m == 4) {
         # tpr.all is too time consuming
         x.variant = 'ddp.all'
       } else {
         x.variant = 'ddp.all'
       }
     }else{stop("Unkown Variant, should be \'ADP\' or \'DDP\' " )}
     
     
     ret_raw = .hhg.test.udfree(x = x, y = y, variant = x.variant, K = m, correct.mi.bias = correct.mi.bias,w.sum = w.sum,w.max = w.max)
     if (max_not_available) {
       ret_raw$max.chisq = NA
       ret_raw$max.lr    = NA
     }
     
     stat.sl[i] = ret_raw$sum.lr
     stat.sc[i] = ret_raw$sum.chisq
     stat.ml[i] = ret_raw$max.lr
     stat.mc[i] = ret_raw$max.chisq
   } 
  }
  
  #arrange output:
  names(stat.sl) = paste0('m.',(mmin:mmax))
  names(stat.sc) = paste0('m.',(mmin:mmax))
  
  names(stat.ml) = paste0('m.',(mmin:mmax))
  names(stat.mc) = paste0('m.',(mmin:mmax))
  
  flag_SL = (score.type=='LikelihoodRatio' || score.type=='both') & ( aggregation.type=='sum' || aggregation.type=='both')
  flag_SC = (score.type=='Pearson' || score.type=='both') & ( aggregation.type=='sum' || aggregation.type=='both')
  flag_ML = (score.type=='LikelihoodRatio' || score.type=='both') & ( aggregation.type=='max' || aggregation.type=='both')
  flag_MC = (score.type=='Pearson' || score.type=='both') & ( aggregation.type=='max' || aggregation.type=='both')
  ret=list()
  if(flag_SL + flag_SC + flag_ML + flag_MC==1){
      if(flag_SL){
        ret$statistic = stat.sl
      }else if(flag_SC){
        ret$statistic = stat.sc
      }else if(flag_ML){
        ret$statistic = stat.ml
      }else if(flag_MC){
        ret$statistic = stat.mc
      }
  }else{
    if(flag_SL){
      ret$sum.lr = stat.sl
    }
    if(flag_SC){
      ret$sum.chisq = stat.sc
    }
    if(flag_ML){
      ret$max.lr = stat.ml
    }
    if(flag_MC){
      ret$max.chisq = stat.mc
    }
  }
  ret$type = type
  ret$stat.type = 'Independence-Stat'
  ret$size = length(x)
  ret$variant = variant
  ret$score.type = score.type
  ret$aggregation.type = aggregation.type
  ret$mmin = mmin
  ret$mmax = mmax
  ret$additional = c(correct.mi.bias , w.sum,w.max)
  class(ret)='UnivariateStatistic'
  return (ret)
}

#function for performing validation checks on inputs. the logic of this function is being shared by many other functions,
#but other functions also have specific checks
.hhg.univariate.check.inputs=function(type,variant,size,score.type,aggregation.type){
  if(is.null(variant)){
    stop("unknown variant - for independence should be 'ADP' or 'DDP'")}
  
  if(is.null(score.type)){
    stop("score.type should be \'LikelihoodRatio\' , \'Pearson\'  or \'both\' ")
    
  }
  if(is.null(aggregation.type)){
    stop("aggregation.type should be \'sum\' , \'max\'  or \'both\' ")
  }
  if(type=='Independence'){
    if( length(size)>1  ){stop("For independence - group size should be single number")}
    
    if(!is.element(score.type,c('LikelihoodRatio','Pearson','both'))){
      stop("score.type should be \'LikelihoodRatio\' , \'Pearson\'  or \'both\' ")
    }
    
    if(!is.element(aggregation.type,c('sum','max','both'))){
      stop("aggregation.type should be \'sum\' , \'max\'  or \'both\' ")
    }
    
    if(!is.element(variant,c('DDP','ADP'))){
      stop("unknown variant - for independence should be 'ADP' or 'DDP'")}
      
  }else if(type=='KSample'){
    if(!is.element(variant,c('KSample-Variant','KSample-Equipartition'))){
      stop("unknown variant - for independence should be 'KSample-Variant' or 'KSample-Equipartition'")}
    if(length(size)<2 ){
      stop("For k-sample - data should have at least two samples")}
    if(!is.element(score.type,c('LikelihoodRatio','Pearson','both'))){
      stop("score.type should be \'LikelihoodRatio\' , \'Pearson\'  or \'both\' ")
    }
    if(!is.element(aggregation.type,c('sum','max','both'))){
      stop("aggregation.type should be \'sum\' , \'max\'  or \'both\' ")
    }
    
  }else{
    stop('Unkown type, should be either \'Independence\' or \'KSample\' ')
  }
}

# MinP Object constructor from m.stats. Function constructs a univariate object, with ecdfs over the specific m distributions,
#but also ecdf for minp and fisher statistics with the specific m
#in order to handle ties in distributions, pvalues are computed using ecdf over the negative values (of single m stattistics, and Fisher)
# with an addition of a value of -Inf. MinP ecdf is attached a value of 0.
.hhg.univariate.object.constructor.from.mstats=function(m.stats,minm,maxm,type,variant,size,score.type,aggregation.type,additional=NULL){
  .hhg.univariate.check.inputs(type,variant,size,score.type,aggregation.type)
  MESSAGES=F #flag used to show messages while computing the null distribution. currently has a value of F hardcoded.
  if(min(maxm)-minm+1 > dim(m.stats)[2] | min(maxm) <2  | min(maxm)<minm){
    stop("minm should be at least 2,  columns in m.stats should be the statistics for null table from, minm to maxm")
  }
  
  if(type =='Independence'){
     #additional checks required? - not currently, all implemented in the checks function
    }
  
  if(type =='KSample'){
    if(sum(size) <min(maxm)){
      stop("size should be larger than maximum number of partitions")
    }
  }
  
  #mode switches - shouldn't be touched
  
  lookup.table=m.stats
  ret=list()
  
  m.bound=max(maxm)
  nr.mmax=length(maxm)
  
  # used to sort and count the ecdf by m (number of partitions)
  m.stat.negative.ecdf=list()
  m.stat.size = list()
  
  #These are matrices for holding the results for many computed m!
  MinP=matrix(rep(NA,nr.mmax*nrow(lookup.table)),nrow =nrow(lookup.table),ncol = nr.mmax)
  m_chosen=matrix(rep(NA,nr.mmax*nrow(lookup.table)),nrow =nrow(lookup.table),ncol = nr.mmax)
  Fisher=matrix(rep(0,nr.mmax*nrow(lookup.table)),nrow =nrow(lookup.table),ncol = nr.mmax)
  
  MinP_ECDF_LIST=list()
  FISHER_NEGATIVE_ECDF_LIST=list() #for higher alternative , one should use an ecdf on the negative to handle ties (conservativly).
  
  
  p.val=matrix(rep(NA,nrow(lookup.table)*(m.bound-minm+1)),nrow = nrow(lookup.table),ncol = m.bound - minm + 1)
  
  for(ci in 1:(m.bound - minm + 1)){ 
    if(MESSAGES){
      print(paste("Sorting statistics, m-wise:",as.character(ci+minm-1)))
    }
    #computing the ecdf over the unique values
    m.stat.negative.ecdf[[ci]] = ecdf(-1*c(lookup.table[,ci],Inf))
    m.stat.size[[ci]] = length(lookup.table[,ci] +1)
    
    #compte the null ranks
    p.val[,ci]=m.stat.negative.ecdf[[ci]](-1* lookup.table[,ci])
    
  }
  #Compute the actual statistics
  for(mi in 1:nr.mmax){
      mim=maxm[mi] # the current k
      if(MESSAGES){      print(paste("Computing Statistic for M.max:",mim))    }
      
      if(MESSAGES){     print("Finding fields with NA")   }
      ind.na=is.na(lookup.table[,(mim-minm+1)]) #find N.A.'s
      if(MESSAGES){     print("Calling quick compute function")   }
      qc=.hhg.univariate.object.fast.compute.statistic(p.val,mmax=mim,mmin=minm,ind.na,MESSAGES)
      MinP[,mi]=qc$MinP
      m_chosen[,mi]=qc$m_chosen
      Fisher[,mi]=qc$Fisher
      MinP_ECDF_LIST[[mi]]=ecdf(c(qc$MinP,0))
      FISHER_NEGATIVE_ECDF_LIST[[mi]]=ecdf(c(-qc$Fisher,-Inf))
  }    
  
  #create the object
  ret$m.stat.negative.ecdf = m.stat.negative.ecdf
  ret$m.stat.size = m.stat.size
  
  
  ret$MinP=MinP
  ret$m_chosen=m_chosen
  ret$Fisher=Fisher
  ret$MinP_ECDF_LIST=MinP_ECDF_LIST
  ret$FISHER_NEGATIVE_ECDF_LIST=FISHER_NEGATIVE_ECDF_LIST
    
  ret$type=type
  ret$variant=variant
  ret$size=size
  ret$score.type=score.type
  ret$aggregation.type=aggregation.type
  ret$minm = minm 
  ret$maxm = maxm
  ret$additional=additional
  class(ret)='UnivariateObject'
  
  return(ret)
}

#function for printing the univariate object
print.UnivariateObject = function(x,...){
  cat(paste0('Univariate Null Table Object \n'))
  cat(paste0('Type: ',x$type,'\n'))
  cat(paste0('Variant: ',x$variant,'\n'))
  if(x$type == 'Independence'){
    cat(paste0('Sample Size: ',x$size,'\n'))  
  }else if(x$type == 'KSample'){
    cat(paste0('Group Sizes: ','\n'))  
    cat(paste0(x$size))
    cat('\n')
  }
  cat(paste0('Score Type: ',x$score.type,'\n'))
  cat(paste0('Aggregation Type: ',x$aggregation.type,'\n'))
  cat(paste0('mmin: ',x$minm,'\n'))
  cat(paste0('mmax: ',x$maxm,'\n'))
  cat(paste0('additional: ','\n'))
  cat((x$additional))
}

#function for printing the univariate statistic object.
print.UnivariateStatistic = function(x,...){
  if(x$stat.type == 'Independence-Stat'){
    .print.ind.stat(x)
  }else if(x$stat.type == 'KSample-Stat'){
    .print.ks.stat(x)
  }else if(x$stat.type == 'Independence-Combined'){
    .print.ind.combined(x)
  }else if(x$stat.type == 'KSample-Combined'){
    .print.ks.combined(x)
  }

}

#auxliary function used for printing hhg.univariate.ind.stat
.print.ind.stat=function(x){
  cat(paste0('HHG univariate independence statistic of type: \n'))
  agg.text = x$aggregation.type
  if(agg.text == 'both'){
    agg.text = 'Both sum & max'
  }
  score.text = x$score.type
  if(score.text == 'both'){
    score.text = 'both Likelihood Ratio & Pearson'
  }
  if(score.text == 'LikelihoodRatio'){
    score.text = 'Likelihood Ratio'
  }
  cat(paste0(agg.text,' of ',x$variant, ' on ',score.text, ' scores.\n\n'))  
  cat(paste0('Minimum partition size: ',x$mmin,'  Maximum partition size: ',x$mmax,' \n\n'))
  cat(paste0('Sample size: ',x$size,' \n\n'))
  if('statistic' %in% names(x)){
    cat(paste0('Statistics, by partition size: \n'))
    print(x$statistic)
  }
  if('sum.lr' %in% names(x)){
    cat(paste0('Sum over Likelihood Ratio scores of Partitions, by partition size: \n'))
    print(round(x$sum.lr,digits = 3))
    cat('\n')
    
  }
  if('max.lr' %in% names(x)){
    cat(paste0('Maximum over Likelihood Ratio scores of Partitions, by partition size: \n'))
    print(round(x$max.lr,digits = 3))
    cat('\n')
    
  }
  if('sum.chisq' %in% names(x)){
    cat(paste0('Sum over Pearson Chi Square scores of Partitions, by partition size: \n'))
    print(round(x$sum.chisq,digits = 3))
    cat('\n')
    
  }
  if('max.chisq' %in% names(x)){
    cat(paste0('Maximum over Pearson Chi Square scores of Partitions, by partition size:  \n'))
    print(round(x$max.chisq,digits = 3))
    cat('\n')
    
  }
}

#auxliary function used for printing hhg.univariate.ks.stat
.print.ks.stat = function(x){
  cat(paste0('HHG univariate ksample statistic of type: \n'))
  agg.text = x$aggregation.type
  if(agg.text == 'both'){
    agg.text = 'Both sum & max'
  }
  score.text = x$score.type
  if(score.text == 'both'){
    score.text = 'both Likelihood Ratio & Pearson'
  }
  if(score.text == 'LikelihoodRatio'){
    score.text = 'Likelihood Ratio'
  }
  cat(paste0(agg.text,' of ',score.text, ' scores over possible partitions.\n\n'))  
  cat(paste0('Minimum partition size: ',x$mmin,'  Maximum partition size: ',x$mmax,' \n\n'))
  cat(paste0('Sample size, by groups: ',' \n'))
  cat(paste0(x$size))
  cat('\n\n')
  if('statistic' %in% names(x)){
    cat(paste0('Statistics, by partition size: \n'))
    print(x$statistic)
  }
  if('sum.lr' %in% names(x)){
    cat(paste0('Sum over Likelihood Ratio scores of Partitions, by partition size: \n'))
    print(round(x$sum.lr,digits = 3))
    cat('\n')
    
  }
  if('max.lr' %in% names(x)){
    cat(paste0('Maximum over Likelihood Ratio scores of Partitions, by partition size: \n'))
    print(round(x$max.lr,digits = 3))
    cat('\n')
    
  }
  if('sum.chisq' %in% names(x)){
    cat(paste0('Sum over Pearson Chi Square scores of Partitions, by partition size: \n'))
    print(round(x$sum.chisq,digits = 3))
    cat('\n')
    
  }
  if('max.chisq' %in% names(x)){
    cat(paste0('Maximum over Pearson Chi Square scores of Partitions, by partition size:  \n'))
    print(round(x$max.chisq,digits = 3))
    cat('\n')
    
  }
}

#auxliary function used for printing hhg.univariate.ind.combined.test
.print.ind.combined = function(x){
  cat(paste0('HHG univariate combined independence statistic\n'))
  cat('Statistics type combined: \n')
  agg.text = x$aggregation.type
  if(agg.text == 'both'){
    agg.text = 'Both sum & max'
  }
  score.text = x$score.type
  if(score.text == 'both'){
    score.text = 'both Likelihood Ratio & Pearson'
  }
  if(score.text == 'LikelihoodRatio'){
    score.text = 'Likelihood Ratio'
  }
  cat(paste0(agg.text,' of ',x$variant, ' on ',score.text, ' scores.\n\n'))  
  cat(paste0('Minimum partition size: ',x$mmin,'  Maximum partition size: ',x$mmax,' \n\n'))
  cat(paste0('Sample size: ',x$size,' \n\n'))
  
  cat('Single m (partition size) statistics:\n')
  print(round(x$m.stats,digits = 4))
  cat('\n')
  
  cat('Single m (partition size) pvalues:\n')
  cat(round(x$pvalues.of.single.m,digits = 4))
  cat('\n\n')
  
  if('MinP' %in% names(x)){
    cat('MinP Statistic - Test statistic is minimum of above single m pvalues: ')
    cat(round(x$MinP,digits = 4))
    cat('\n\n')
    
    cat('Partition size with minimum p-value (# of cells):\n')
    cat(round(x$MinP.m.chosen,digits = 4))
    cat('\n\n')
    
    cat('p-value for MinP test:')
    cat(round(x$MinP.pvalue,digits = 4))
    cat('\n\n')
  }
  
  if('Fisher' %in% names(x)){
    cat('Fisher Statistic - Test statistic is -sum(log(pval)), over  single m pvalues: ')
    cat(round(x$Fisher,digits = 4))
    cat('\n\n')
    
    cat('p-value for Fisher test:')
    cat(round(x$Fisher.pvalue,digits = 4))
    cat('\n\n')
  }
}

#auxliary function used for printing hhg.univariate.ind.ks.test
.print.ks.combined = function(x){
  cat(paste0('HHG univariate combined K-sample statistic\n'))
  cat('Statistics type combined: \n')
  agg.text = x$aggregation.type
  if(agg.text == 'both'){
    agg.text = 'Both sum & max'
  }
  score.text = x$score.type
  if(score.text == 'both'){
    score.text = 'both Likelihood Ratio & Pearson'
  }
  if(score.text == 'LikelihoodRatio'){
    score.text = 'Likelihood Ratio'
  }
  cat(paste0(agg.text,' of ',score.text, ' scores over possible partitions.\n\n'))  
  cat(paste0('Minimum partition size: ',x$mmin,'  Maximum partition size: ',x$mmax,' \n\n'))
  cat(paste0('Sample size, by groups: ',' \n'))
  cat(paste0(x$size))
  cat('\n\n')
  
  cat('Single m (partition size) statistics:\n')
  print(round(x$m.stats,digits = 4))
  cat('\n')
  
  cat('Single m (partition size) pvalues:\n')
  cat(round(x$pvalues.of.single.m,digits = 4))
  cat('\n\n')
  
  if('MinP' %in% names(x)){
    cat('MinP Statistic - Test statistic is minimum of above single m pvalues: ')
    cat(round(x$MinP,digits = 4))
    cat('\n\n')
    
    cat('Partition size with minimum p-value(# of cells):\n')
    cat(round(x$MinP.m.chosen,digits = 4))
    cat('\n\n')
    
    cat('p-value for MinP test:')
    cat(round(x$MinP.pvalue,digits = 4))
    cat('\n\n')
  }
  
  if('Fisher' %in% names(x)){
    cat('Fisher Statistic - Test statistic is -sum(log(pval)), over  single m pvalues: ')
    cat(round(x$Fisher,digits = 4))
    cat('\n\n')
    
    cat('p-value for Fisher test:')
    cat(round(x$Fisher.pvalue,digits = 4))
    cat('\n\n')
  }
}

#Internal Function for fast computing MinP & Fisher statistics by p.values (works on a matrix of values, not just a single observations)
.hhg.univariate.object.fast.compute.statistic=function(p.val,mmax,mmin=2,ind.na,MESSAGES=F){
  nr.row=nrow(p.val)
  if(MESSAGES){ print(paste("Fast Compute, # of rows is: ",nr.row,sep=''))}
  m.bound=max(mmax)
  nr.kmax=length(mmax)
  
  m_chosen=rep(NA,nr.row)
  Statistic=rep(NA,nr.row)
  Fisher=rep(NA,nr.row)
  
  if(MESSAGES){      print("Generating : finding indexes of mins")    }
  if(mmax-mmin>0){
    ind=apply(p.val[,1:(mmax-mmin+1)],1,which.min) # before na check
  }else{
    ind=rep(2,nrow(p.val))
  }
  if(MESSAGES){     print("Generating : writing results")   }
  for(ri in 1:nr.row){
    
    if(mmax>2){
      Statistic[ri]=p.val[ri,ind[ri]]
      Fisher[ri]=(-1)*sum(log(p.val[ri,(1 ):(mmax-mmin+1)]))
    }else{
      Statistic[ri]=p.val[ri]
      Fisher[ri]=(-1)*sum(log(p.val[ri]))
    }
  }
  m_chosen=ind+1
  Statistic[ind.na]=NA #handle na's
  m_chosen[ind.na]=NA
  Fisher[ind.na]=NA
  
  ret=list()
  ret$MinP=Statistic
  ret$Fisher=Fisher
  ret$m_chosen=m_chosen
  return(ret)
}

#function for constructing a null table object for independece statistics
hhg.univariate.ind.nulltable=function(size,mmin=2,mmax = max(floor(sqrt(size)/2),2),variant = 'ADP',aggregation.type = 'sum',score.type='LikelihoodRatio', w.sum = 0, w.max = 2,nr.replicates=1000,keep.simulation.data=F){  
  correct.mi.bias = F #mi estimation correction not supported yet
  type='Independence'
  
  #input checks:
  
  if(is.null(size) |any(is.null(mmax)) | is.null(nr.replicates) | is.null(variant) ){
    stop("size, mmax, nr.replicates should be integers. mmax and mmin should be less than size and >=2")
  }
  
  if(is.na(size) | is.na(nr.replicates) |is.na(variant) |any(is.na(mmax))){
    stop("size, Max.m, nr.replicates should be integers. mmax and mmin should be less than size and >=2")
  }
  
  
  for(i in 1:length(mmax)){
    if(mmax[i]!=as.integer(mmax[i]) & !is.na(mmax[i])){
      stop("Max.m must be vector of integers")
    }
  }  
  
  .hhg.univariate.check.inputs(type,variant,size,score.type,aggregation.type)
  
  if(size!= as.integer(size) |  nr.replicates!=as.integer(nr.replicates) | size<1 | max(mmax) > size | min(mmax) <2 | mmin<2| nr.replicates < 1){
    stop("size, mmax, nr.replicates should be integers. mmax shoudld be less than size and >=2")
  }
  
  if(score.type=='both' | aggregation.type=='both'){
    stop('null table can only be constructed of a single statistic, \'both\' not allowed in score.type or aggregation.type')
  }
  
  #create null table statistics
  m.stats=matrix(rep(NA,(max(mmax)-mmin+1)*nr.replicates),nrow = nr.replicates,ncol = (max(mmax)-mmin+1))
  x=(1:size)
  y_samp=(1:size)
  for(r in 1:nr.replicates){
    y=sample(y_samp)
    current_test=hhg.univariate.ind.stat(x,y,variant,aggregation.type,score.type,mmax = max(mmax),mmin = mmin,w.sum = w.sum,w.max = w.max)
    m.stats[r,]=current_test$statistic
  }
  
  #create objects
  colnames(m.stats)=paste0('m.',(mmin:max(mmax)))
  univariate.object = .hhg.univariate.object.constructor.from.mstats(m.stats,mmin,mmax,type,variant,size,score.type,aggregation.type,c(correct.mi.bias,w.sum,w.max))
  ret=list()
  if(keep.simulation.data){
    ret$m.stats = m.stats  
  }
  ret$univariate.object = univariate.object
  return(ret)
}

# Perform MinP/Fisher Test from m.stats
.hhg.univariate.compute.from.mstats=function(m.stats,univariate.object,minm,maxm){
  
  #inputs checks:
  
  if( minm!=as.integer(minm) | maxm!=as.integer(maxm)){
    stop("maxm and minm must be integers must be an integer")
  }
  if(maxm-minm+1>length(m.stats)){
    s=paste("m.stats provided is shorter than required maxm. # m.stats provided:",length(m.stats), "(maybe due to clumps) maxm:",maxm,' minm:',minm,sep='')
    stop(s)
  }
  if(max(univariate.object$maxm)<maxm){
    stop("Univariate object's maximum number of partitions computed is not sufficient")
  }
  if((univariate.object$minm)!=minm){
    stop("Univariate object's minimum number of partitions computed is not same as required in statistic computation")
  }
  
  #compute:
  ret=list()
  pm=rep(NA,maxm-minm+1)
  if(is.na(m.stats[length(m.stats)])){ #cannot perform test 
    stop("Missing values in m.stats")
  }
  for(pi in 1:(maxm-minm+1)){
    #compute the p.vals of the m.stats by the MinP.Object
    pm[pi]=univariate.object$m.stat.negative.ecdf[[pi]](-1*m.stats[pi])
  }
  
  #create object and return
  ret$MinP.m.chosen=which.min(pm)+1 #find the minimum
  ret$MinP=pm[ret$MinP.m.chosen-1]
  ret$pvalues.of.single.m=pm
  ret$Fisher=(-1)*sum(log(pm[1:(maxm-minm+1)]))
  if(sum(univariate.object$maxm==maxm)>0){
    current_m=which(univariate.object$maxm==maxm)
    ret$MinP.pvalue=univariate.object$MinP_ECDF_LIST[[current_m]](ret$MinP)
    ret$Fisher.pvalue=univariate.object$FISHER_NEGATIVE_ECDF_LIST[[current_m]](-ret$Fisher)
  }
  return(ret)
}

#function for compiling a null table object from existing m.stats
#(like the ones given on the site, or the ones used in the code example for function and vignette, for computing a null table using multiple cores)
hhg.univariate.nulltable.from.mstats=function(m.stats,minm,maxm,type,variant,size,score.type,aggregation.type, w.sum = 0, w.max = 2, keep.simulation.data=F){
  correct.mi.bias = F #not supported yet
  ret=list()
  #call the private function to create and object, and return it. all inputs checks are to be dont inside
  ret$univariate.object = .hhg.univariate.object.constructor.from.mstats(m.stats,minm,maxm,type,variant,size,score.type,aggregation.type,additional =c(correct.mi.bias,w.sum,w.max))
  if(keep.simulation.data){
    ret$m.stats=m.stats
  }
  return(ret)
}

#function for performing the combined independence test, combining different partition sizes.
hhg.univariate.ind.combined.test=function(X,Y=NULL,NullTable=NULL,mmin=2,mmax=max(floor(sqrt(length(X))/2),2),variant='ADP',aggregation.type='sum',score.type='LikelihoodRatio', w.sum = 0, w.max = 2 ,combining.type='MinP',nr.perm=100){
  
  correct.mi.bias = F #not supported yet
  # function gets either a vector at X and Y, or the result of hhg.univariate.ind.stat at X and null Y. function
  # may also receive null table object or create null table by parameters. since exiting data and null table need to work toghether,
  # or null table and hhg.univariate.ind.stat result, the function has the following priorities:
  
  #  if no null table
  #     if statisticc is found
  #         take statistic parameters
  #     else
  #         take parameters from function inputs
  #  else (there is a null table)
  #     if there is a statistic, 
  #       take its parameters.
  #       try to find an mmax in the null table, that fits statistic.
  #       if not found, take the min mmax (it could be the statistic will still not have enough terms, or it wont use all   		terms)
  #     else
  #       take null table parameters
  
    #these variables house the parameters chosen by the above logic.
    current_null_table=NULL
    current_mmax=NULL
    current_mmin=NULL
    current_variant=NULL
    current_aggregation.type=NULL
    current_score.type=NULL
    current_correct.mi.bias = NULL
    current_size = NULL
    current_w.sum = NULL
    current_w.max = NULL
    XY_is_stat=F
    m.stats=rep(NA,mmax-mmin+1) #if Y is null and X is statistic, this doesn't matter, since it is overridden (it is of a false dimensionality)
    
    #input checks
    if(!is.element(combining.type,c('MinP','Fisher','Both'))){
      stop('combining.type should be \'MinP\' , \'Fisher\'  or \'Both\' ')  
    }
    if(is.null(Y)){ if(is.list(X)){ if( class(X) =='UnivariateStatistic'){ if(X$stat.type=='Independence-Stat'){
      XY_is_stat=T #we already acctually have the statistics
    }}}}
    if(XY_is_stat == F & is.null(Y)){
      stop('Y not supplied or X is not valid result of hhg.univariate.ind.stat.')
    }
    if(is.null(NullTable)){ #handle the case of no null table found, generate, according to parameters. remember to take parameters from object if found.
      current_mmax = mmax
      current_mmin = mmin
      current_size = length(X)
      current_variant = variant
      current_aggregation.type = aggregation.type
      current_score.type = score.type
      current_correct.mi.bias = correct.mi.bias
      current_w.sum = w.sum
      current_w.max = w.max
      current_size = length(X)
      if(XY_is_stat){ #handle object
        current_mmax=X$mmax
        current_mmin=X$mmin
        current_variant=X$variant
        current_size = X$size
        current_aggregation.type = X$aggregation.type
        current_score.type = X$score.type
        current_correct.mi.bias = X$additional[1]
        current_w.sum = X$additional[2]
        current_w.max = X$additional[3]
      }else{
        if(length(X) != length(Y)){
          stop('X and Y not of the same length')
        }
      }
      .hhg.univariate.check.inputs('Independence',current_variant,current_size,current_score.type,current_aggregation.type)
      
      current_null_table = hhg.univariate.ind.nulltable(size= current_size,mmin,mmax,variant,aggregation.type,score.type,nr.replicates = nr.perm, w.sum = w.sum, w.max = w.max,keep.simulation.data = T)
      
    }else{ #null table is found
      current_null_table = NullTable
      #checking its of the right type
      if(!is.list(current_null_table) | !('univariate.object' %in% names(current_null_table))){
        stop('NullTable supplied is not valid. construct null table using hhg.univariate.ind.nulltable')
      }
      if(current_null_table$univariate.object$type!='Independence'){
        stop('null table type is not \'Independence\' ')
      }
      
      current_variant = NullTable$univariate.object$variant
      current_aggregation.type = NullTable$univariate.object$aggregation.type
      current_score.type = NullTable$univariate.object$score.type
      current_mmin = NullTable$univariate.object$minm
      
      
      for(m in current_null_table$univariate.object$maxm){
        if(m == mmax){current_mmax=m}
      }
      if(is.null(current_mmax)){ current_mmax = min(current_null_table$univariate.object$maxm)}

      current_correct.mi.bias = NullTable$univariate.object$additional[1]
      current_w.sum = NullTable$univariate.object$additional[2]
      current_w.max = NullTable$univariate.object$additional[3]
      current_size = NullTable$univariate.object$size
      if(XY_is_stat){ #handle the case we have both a null table and an object
        
        current_mmax=NULL
        for(m in current_null_table$univariate.object$maxm){
          if(m == mmax){current_mmax=m}
        }
        if(is.null(current_mmax)){ current_mmax = min(current_null_table$univariate.object$maxm)}
        
        
        
        if(current_mmax > X$mmax){
          stop('Null Table mmax bigger than computed stat m')
        }
        if(current_mmin != X$mmin){
          stop('Null Table mmin is different than stat m')
        }
        if(current_variant != X$variant){
          stop('Null Table and Statistic not of same variant')
        }
        if(current_size != X$size){
          stop('Null Table and Statistic not of same size')
        }
        if(current_aggregation.type != X$aggregation.type){
          stop('Null Table and Statistic not of same aggregation type')
        }
        if(current_score.type != X$score.type){
          stop('Null Table and Statistic not of same score type')
        }
        if(current_correct.mi.bias != X$additional[1]){
          stop('Null Table and Statistic not of same bias correction')
        }
        if(current_w.sum != X$additional[2]){
          stop('Null Table and Statistic not of same w.sum')
        }
        if(current_w.max != X$additional[3]){
          stop('Null Table and Statistic not of same w.max')
        }
      }else{
      }  
        
    }#end of null table found
    if(XY_is_stat){
      m.stats[1:(current_mmax - current_mmin +1)]  = X$statistic[1:(current_mmax - current_mmin +1)]
    }else{
      if(length(X)!=current_null_table$univariate.object$size){
        stop(paste0('Sample size is ',length(X),' while null table was constructed on sample size ',NullTable$univariate.object$size))
      }        
      m.stats=hhg.univariate.ind.stat(X,Y,current_variant,current_aggregation.type,current_score.type,mmin = current_mmin , mmax = current_mmax,w.sum = current_w.sum , w.max = current_w.max)$statistic        
    }
    #compute statistics and pvalues.
    computed = .hhg.univariate.compute.from.mstats(m.stats,current_null_table$univariate.object,current_mmin,current_mmax)
  
    #pack object
    ret=computed
    ret$m.stats=m.stats
    names(ret$m.stats)=paste0('m.',(current_mmin:current_mmax))
    if(is.null(NullTable)){
      ret$generated_null_table=current_null_table
    }
    if(combining.type == 'MinP'){
      ind.to.remove = which(names(ret) %in% c('Fisher','Fisher.pvalue'))
      ret=ret[-ind.to.remove]
    }else if(combining.type =='Fisher'){
      ind.to.remove = which(names(ret) %in% c('MinP','MinP.pvalue','MinP.m.chosen'))
      ret=ret[-ind.to.remove]
    }else if(combining.type =='Both'){
      #do nothing, no need to remove
    }
    ret$stat.type = 'Independence-Combined'
    ret$mmin = current_mmin
    ret$mmax = current_mmax
    ret$score.type = current_score.type
    ret$aggregation.type = current_aggregation.type
    ret$variant = current_variant
    ret$w.sum = current_w.sum
    ret$w.max = current_w.max
    ret$size = current_size
    class(ret)='UnivariateStatistic'
    return(ret)
}

#function for computing the pvalue of an independence statistic
hhg.univariate.ind.pvalue=function(statistic, NullTable, m=min(statistic$mmax,4)){
  #input checks
  if(!is.list(statistic) | class(statistic)!='UnivariateStatistic'){
    stop(' statistic is not result of hhg.univariate.ind.stat')
  }
  if(statistic$stat.type != 'Independence-Stat'){
    stop(' statistic is not result of hhg.univariate.ind.stat')
  }
  if(!is.list(NullTable) | !('univariate.object' %in% names(NullTable))){
    stop('NullTable supplied is not valid. construct null table using hhg.univariate.ind.nulltable')
  }
  
  uvo=NullTable$univariate.object
  if(class(uvo)!= 'UnivariateObject'){
    stop('Null table supplied does not contain valid univariate object')
  }
  if(length(statistic$size) != length(uvo$size)){
    stop(paste0('statistic group sizes are ',statistic$size,' while null table group sizes are:',uvo$size ))
  }
  if(any(statistic$size != uvo$size)){
    stop(paste0('statistic size is ',statistic$size,' while null table size is ', uvo$size))
  }
  if(statistic$variant != uvo$variant){
    stop(paste0('statistic variant is ',statistic$variant,' while null table variant is ', uvo$variant))
  }
  if(statistic$aggregation.type != uvo$aggregation.type){
    stop(paste0('statistic aggregation.type is ',statistic$aggregation.type,' while null table aggregation.type is ', uvo$aggregation.type))
  }
  if(statistic$score.type != uvo$score.type){
    stop(paste0('statistic score.type is ',statistic$score.type,' while null table score.type is ', uvo$score.type))
  }
  if(statistic$type != uvo$type){
    stop(paste0('statistic type is ',statistic$type,' while null table type is ', uvo$type))
  }
  
  
  if(statistic$additional[1] != uvo$additional[1] ){
    stop(paste0('MI bias correction in test sample is : ',statistic$additional[1] , ' while in null table is : ',uvo$additional[1]))
  }
  if(TRUE){#if(statistic$additional[1]){
    if(statistic$aggregation.type=='max' & statistic$additional[3] != uvo$additional[3]){
      stop(paste0('w.max in test sample: ', statistic$additional[3] , ' while in null table: ', uvo$additional[3]))
    }
    if(statistic$aggregation.type=='sum' & statistic$additional[2] != uvo$additional[2]){
      stop(paste0('w.sum in test sample: ', statistic$additional[2] , ' while in null table: ', uvo$additional[2]))
    }
  }
  
  
  ret=list()
  if(statistic$stat.type == 'Independence-Stat'){
    if(uvo$minm >m | max(uvo$maxm) < m){
      stop(paste0('null table does not include number of partitions m: ',m))
    }
    
    if(statistic$mmin >m | max(statistic$mmax) < m){
      stop(paste0('statistic  does not include number of partitions m: ',m))
    }
    #pvalue can be computed
    m.ind = m - uvo$minm +1
    m.stat.ind = m - statistic$mmin +1
    stat=statistic$statistic[m.stat.ind] #the single m statistic value under the statistic object
    return(uvo$m.stat.negative.ecdf[[m.ind]](-1 * stat))

  }else{
    stop('incorrect statistic type: statistic should be the result of hhg.univariate.ind.stat')
  }  
}

#function for computing K sample statsitic
hhg.univariate.ks.stat=function(x, y,aggregation.type='sum',score.type='LikelihoodRatio', mmax = max(4,round(min(table(y))/3)),mmin=2){
  variant = 'KSample-Variant' #the only type of variant
  type='KSample'
  #input checks
  if(is.null(x) |is.null(y)){
    stop("x & y should be vectors of doubles." )
  }
  if(length(x)!= length(y)){
    stop("X & Y not the same length")
  }
  if(!is.numeric(mmin) | !is.numeric(mmax)){
    stop(" parameters should be 2<= mmin<= mmax, integers")  
  }
  if(mmin !=round(mmin) | mmax !=round(mmax)){
    stop(" parameters should be 2<= mmin<= mmax, integers")  
  }
  if(mmin<2 | mmax<mmin){
    stop(" parameters should be 2<= mmin<= mmax, integers")  
  }
  
  .hhg.univariate.check.inputs(type,variant,(table(y)),score.type,aggregation.type)
  
  #check which types should be computed
  flag_SL = (score.type=='LikelihoodRatio' || score.type=='both') & ( aggregation.type=='sum' || aggregation.type=='both')
  flag_SC = (score.type=='Pearson' || score.type=='both') & ( aggregation.type=='sum' || aggregation.type=='both')
  flag_ML = (score.type=='LikelihoodRatio' || score.type=='both') & ( aggregation.type=='max' || aggregation.type=='both')
  flag_MC = (score.type=='Pearson' || score.type=='both') & ( aggregation.type=='max' || aggregation.type=='both')
  results_SL=NULL
  results_SC=NULL
  results_ML=NULL
  results_MC=NULL
  if(mmax==mmin){
    if(flag_ML| flag_MC){
      test=.hhg.univariate.ks.stat.max(x, y, variant = 'Mm',Max.m =mmax ,m.stats.wanted=T)
      if(length(test$m.stats)< 2*(mmax-1)){ #this is due to use of clumps in the original max computation algorithem
        stop(paste0('Too many clumps in data, cannot compute ', mmax ,' partitions'))
      }
      results_ML=test$m.stats[(mmin-1)] 
      results_MC = test$m.stats[(mmax-1)+(mmin-1)]
    }
    if(flag_SC | flag_SL){
      test=.xdp.test.k.sample(x,y,mmax)
      results_SL = test$sum.lr
      results_SC = test$sum.chisq
    }
  }else{
    if(flag_ML | flag_MC){
      test=.hhg.univariate.ks.stat.max(x, y, variant = 'Mm',Max.m =mmax ,m.stats.wanted=T) #computing the Mm statistic.
      results_ML=test$m.stats[(mmin-1):(mmax-1)] 
      names(results_ML)=paste('M.m_',(mmin:(mmax)),sep = '')  
      results_ML=results_ML[!is.na(results_ML)]
      results_MC=test$m.stats[(mmax-1) + (mmin-1):(mmax-1)] 
      names(results_MC)=paste('M.m_',(mmin:(mmax)),sep = '')  
      results_MC=results_MC[!is.na(results_MC)]
    }
    if(flag_SC | flag_SL){
      test=.xdp.test.k.sample.mk(x,y,mmax) #computing the Sm statistic
      results_SC = test[(mmin-1):(mmax-1),1]
      names(results_SC)=paste('S.m_',(mmin:mmax),sep = '')  
      results_SL = test[(mmin-1):(mmax-1),2]
      names(results_SL)=paste('S.m_',(mmin:mmax),sep = '')  
    }
  }
  #writing to output object
  ret=list()
  if(flag_SL + flag_SC + flag_ML + flag_MC==1){
    if(flag_SL){
      ret$statistic = results_SL
    }else if(flag_SC){
      ret$statistic = results_SC
    }else if(flag_ML){
      ret$statistic = results_ML
    }else if(flag_MC){
      ret$statistic = results_MC
    }
  }else{
    if(flag_SL){
      ret$sum.lr = results_SL
    }
    if(flag_SC){
      ret$sum.chisq = results_SC
    }
    if(flag_ML){
      ret$max.lr = results_ML
    }
    if(flag_MC){
      ret$max.chisq = results_MC
    }
  }
  ret$type = type
  ret$variant=variant
  ret$stat.type = 'KSample-Stat'
  ret$size = table(y)
  ret$score.type = score.type
  ret$aggregation.type = aggregation.type
  ret$mmin = mmin
  ret$mmax = mmax
  class(ret)='UnivariateStatistic'
  return (ret)
}

#function for creating K sample null table object
hhg.univariate.ks.nulltable=function(group.sizes,mmin=2,mmax=max(4,round(min(group.sizes)/3)),aggregation.type='sum',score.type='LikelihoodRatio',nr.replicates=1000,keep.simulation.data=F){
  type='KSample'
  variant = 'KSample-Variant'
  #check inputs
  if(any(is.null(group.sizes)) |any(is.null(mmax)) | is.null(nr.replicates) | is.null(variant) ){
    stop("Size should be a vector of valid group sizes (more than two groups). mmax,mmin , nr.replicates should be integers. mmax and mmin should be less than sum(group.sizes) and >=2")
  }
  
  if(any(is.na(group.sizes)) | is.na(nr.replicates) |is.na(variant) |any(is.na(mmax))){
    stop("Size should be a vector of valid group sizes (more than two groups). mmax,mmin, nr.replicates should be integers. mmax and mmin should be less than sum(group.sizes) and >=2")
  }
  
  if(any(mmax!=as.integer(mmax)) & any(is.na(mmax))){
    stop("mmax must be vector of integers")
  }
  
  .hhg.univariate.check.inputs(type,variant,group.sizes,score.type,aggregation.type)
  
  if(any(group.sizes!= as.integer(group.sizes)) |  nr.replicates!=as.integer(nr.replicates) | any(group.sizes<1) | max(mmax) > sum(group.sizes) | mmin<2 |min(mmax) <2 | nr.replicates < 1){
    stop("Size should be a vector of valid group sizes (more than two groups). mmax,mmin nr.replicates should be integers. mmax and mmin should be less than sum(group.sizes) and >=2")
  }
  
  if(score.type=='both' | aggregation.type=='both'){
    stop('null table can only be constructed of a single statistic, \'both\' not allowed in score.type or aggregation.type')
  }
  #perform permutations on ranks
  m.stats=matrix(rep(NA,(max(mmax)-mmin+1)*nr.replicates),nrow = nr.replicates,ncol = (max(mmax)-mmin+1))
  x=(1:sum(group.sizes))
  y_samp=NULL
  
  for(i in 1:length(group.sizes)){
    y_samp=c(y_samp,rep(i-1,group.sizes[i]))
  }
  
  for(r in 1:nr.replicates){
    y=sample(y_samp)
    current_test=hhg.univariate.ks.stat(x,y,aggregation.type,score.type,mmax = max(mmax),mmin=mmin)
    m.stats[r,]=current_test$statistic
  }
  
  #compile statistics matrix to object
  colnames(m.stats)=paste0('m.',(mmin:max(mmax)))
  univariate.object = .hhg.univariate.object.constructor.from.mstats(m.stats,mmin,mmax,type,variant,group.sizes,score.type,aggregation.type,additional = NULL)
  ret=list()
  if(keep.simulation.data){
    ret$m.stats = m.stats  
  }
  ret$univariate.object = univariate.object
  return(ret)
}

#function for performing the combined K sample test, combining different partition sizes.
hhg.univariate.ks.combined.test=function(X,Y=NULL,NullTable=NULL,mmin=2,mmax=ifelse(is.null(Y),4,max(4,round(min(table(Y))/3))) ,aggregation.type='sum',score.type='LikelihoodRatio' ,combining.type='MinP',nr.perm=1000){
  # function gets either a vector at X and Y, or the result of hhg.univariate.ks.stat at X and null Y. function
  # may also receive null table object or create null table by parameters. since exiting data and null table need to work toghether,
  # or null table and hhg.univariate.ind.stat result, the function has the following priorities:
  
#  if no null table
#     if statisticc is found
#         take statistic parameters
#     else
#         take parameters from function inputs
#  else (there is a null table)
#     if there is a statistic, 
#       take its parameters.
#       try to find an mmax in the null table, that fits statistic.
#       if not found, take the min mmax (it could be the statistic will still not have enough terms, or it wont use all 			terms)
#     else
#       take null table parameters
  
  #these hold the final parameters for the test, according to the above decision logic.
  current_null_table=NULL
  current_mmax=NULL
  current_mmin=NULL
  current_variant = 'KSample-Variant'
  current_aggregation.type=NULL
  current_score.type=NULL
  current_size = NULL
  XY_is_stat=F # check if X is an actual UnivariateStatistic object of the required type.
  
  if(!is.numeric(mmax) | is.infinite(mmax)){ #this is to handle a case were one uses Y as null, and X as statistic
    mmax=mmin
  }
  
  m.stats=rep(NA,mmax-mmin+1)
  if(!is.element(combining.type,c('MinP','Fisher','Both'))){
    stop('combining.type should be \'MinP\' , \'Fisher\'  or \'Both\' ')  
  }
  if(is.null(Y)){ if(is.list(X)){ if( class(X) == 'UnivariateStatistic'){ if(X$stat.type=='KSample-Stat'){
    XY_is_stat=T #we already acctually have the statistics
  }}}}
  if(XY_is_stat == F & is.null(Y)){
    stop('Y not supplied or X is not valid result of hhg.univariate.ks.stat.')
  }
  if((is.null(NullTable))){ #no null table is given, we have to generate by user specific parameters.
    current_mmax=mmax
    current_mmin=mmin
    current_aggregation.type = aggregation.type
    current_score.type = score.type
    if(XY_is_stat){ #handle the case no null table is given, but we have an object in X
      current_mmax = X$mmax
      current_mmin = X$mmin
      current_variant = X$variant
      current_size = X$size
      current_aggregation.type = X$aggregation.type
      current_score.type = X$score.type
    }else{
      current_size = table(Y)
      if(length(X) != length(Y)){
          stop('X and Y not of the same length')
      }
      
    }
    .hhg.univariate.check.inputs('KSample',current_variant,current_size,current_score.type,current_aggregation.type)
    #generate null table
    current_null_table = hhg.univariate.ks.nulltable(current_size,current_mmin,current_mmax,current_aggregation.type,current_score.type,nr.replicates = nr.perm,keep.simulation.data = T)
  }else{ #null table is given, we have to check if its of the correct type
    current_null_table = NullTable
    if(!is.list(current_null_table) | !('univariate.object' %in% names(current_null_table))){
      stop('NullTable supplied is not valid. construct null table using hhg.univariate.ks.nulltable')
    }
    if(current_null_table$univariate.object$type != 'KSample'){
      stop('null table type is not \'KSample\' ')
    }
    current_variant = NullTable$univariate.object$variant
    current_aggregation.type = NullTable$univariate.object$aggregation.type
    current_score.type = NullTable$univariate.object$score.type
    current_mmin = NullTable$univariate.object$minm
    for(m in current_null_table$univariate.object$maxm){
      if(m == mmax){current_mmax=m
      }
    }
    if(is.null(current_mmax)){ current_mmax = min(current_null_table$univariate.object$maxm)}
    current_size = NullTable$univariate.object$size

    if(XY_is_stat){ #handle a case where we have both a null table and object
      current_mmax = NULL
      for(m in current_null_table$univariate.object$maxm){
        if(m == X$mmax){current_mmax=m
        }
      }
      if(is.null(current_mmax)){ current_mmax = min(current_null_table$univariate.object$maxm)}
      if(current_mmax > X$mmax){
        stop('Null Table mmax bigger than computed stat m')
      }
      if(current_mmin != X$mmin){
        stop('Null Table mmin is different than stat m')
      }
      if(current_variant != X$variant){
        stop('Null Table and Statistic not of same variant')
      }
      if(length(current_size) != length(X$size)){
        stop('Null Table and Statistic not of same number of groups')
      }
      if(any(sort(current_size) != sort(X$size))){
        stop('Null Table and Statistic not of same group sizes')
      }
      if(current_aggregation.type != X$aggregation.type){
        stop('Null Table and Statistic not of same aggregation type')
      }
      if(current_score.type != X$score.type){
        stop('Null Table and Statistic not of same score type')
      }
    }else{
    }  
    
  }#end of null table found
  if(XY_is_stat){
    m.stats[1:(current_mmax - current_mmin +1)]  = X$statistic[1:(current_mmax - current_mmin +1)]
  }else{
    tab_y=table(Y)
    if(length(current_size)!=length(tab_y)){
      stop('Test sample group sizes and null table group sizes do not match - Not same number of groups')
    }
    if(any(sort(current_size)!=sort(tab_y))){
      stop('Test sample group sizes and null table group sizes do not match!')
    }
    
    m.stats = hhg.univariate.ks.stat(X,Y,current_aggregation.type,current_score.type,current_mmax,current_mmin)$statistic
    
    actual_mmax = current_mmin + length(m.stats)-1 #might we shorter than wanted for dynamic slicing
    if(actual_mmax != current_mmax){# Mm had too few clumps
      ret=list()
      ret$m.stats=m.stats
      return(ret)
      warning('Too few clumps in data, lower mmax or use \'sum\' in aggregation.type')
    }
  }
  
  # compute pvalues and MinP,Fisher statistics
  computed = .hhg.univariate.compute.from.mstats(m.stats,current_null_table$univariate.object,current_mmin,current_mmax)
  
  #pack the results
  ret=computed
  ret$m.stats=m.stats
  if(aggregation.type=='sum'){
    names(ret$m.stats)=paste0('Sm.',current_mmin:current_mmax)
  }else if(aggregation.type=='max'){
    names(ret$m.stats)=paste0('Mm.',current_mmin:current_mmax)
  }
  if(is.null(NullTable)){
    ret$generated_null_table=current_null_table
  }
  if(combining.type == 'MinP'){
    ind.to.remove = which(names(ret) %in% c('Fisher','Fisher.pvalue'))
    ret=ret[-ind.to.remove]
  }
  else if(combining.type =='Fisher'){
    ind.to.remove = which(names(ret) %in% c('MinP','MinP.pvalue','MinP.m.chosen'))
    ret=ret[-ind.to.remove]
  }else if(combining.type =='Both'){
    #do nothing, no need to remove
  }
  ret$stat.type = 'KSample-Combined'
  ret$mmin = current_mmin
  ret$mmax = current_mmax
  ret$score.type = current_score.type
  ret$aggregation.type = current_aggregation.type
  ret$variant = current_variant
  ret$size = current_size
  class(ret)='UnivariateStatistic'
  return(ret)
}

#function for computing p-values of the K sample statistic
hhg.univariate.ks.pvalue=function(statistic, NullTable,m = min(statistic$mmax,4)){
  #input checks
  if(!is.list(statistic) | class(statistic)!='UnivariateStatistic'){
    stop(' statistic is not result of hhg.univariate.ks.stat')
  }
  if(statistic$stat.type != 'KSample-Stat'){
    stop(' statistic is not result of hhg.univariate.ks.stat')
  }
  if(!is.list(NullTable) | !('univariate.object' %in% names(NullTable))){
    stop('NullTable supplied is not valid. construct null table using hhg.univariate.ks.nulltable')
  }
  uvo=NullTable$univariate.object
  if(class(uvo)!= 'UnivariateObject'){
    stop('Null table supplied does not contain valid univariate object')
  }
  if(length(statistic$size)!=length(uvo$size)){
    stop('Test sample group sizes and null table group sizes do not match - not same number of groups')
  }
  if(any(sort(statistic$size)!=sort(uvo$size))){
    stop('Test sample group sizes and null table group sizes do not match!')
  }
  if(statistic$variant != uvo$variant){
    stop(paste0('statistic variant is ',statistic$variant,' while null table variant is ', uvo$variant))
  }
  if(statistic$aggregation.type != uvo$aggregation.type){
    stop(paste0('statistic aggregation.type is ',statistic$aggregation.type,' while null table aggregation.type is ', uvo$aggregation.type))
  }
  if(statistic$score.type != uvo$score.type){
    stop(paste0('statistic score.type is ',statistic$score.type,' while null table score.type is ', uvo$score.type))
  }
  if(statistic$type != uvo$type){
    stop(paste0('statistic type is ',statistic$type,' while null table type is ', uvo$type))
  }

  ret=list()
  if(statistic$stat.type == 'KSample-Stat'){
    if(uvo$minm >m | max(uvo$maxm) < m){
      stop(paste0('null table does not include number of partitions m: ',m))
    }
    if(statistic$mmin >m | statistic$mmax < m){
      stop(paste0('statistic does not include number of partitions m: ',m))
    }
    #compute pvalue
    m.ind = m - uvo$minm +1
    m.stat.ind = m - statistic$mmin +1
    stat=statistic$statistic[m.stat.ind] #the single m statistic value under the statistic object
    return(uvo$m.stat.negative.ecdf[[m.ind]](-1 * stat))
  }else{
    stop('incorrect statistic type: statistic should be the result of hhg.univariate.ks.stat')
  }  
}

#'@title k-sample DS and Mm test. 
#'@description function for performin k sample DS and Mm test. this function is called by the M.m Test Function
#'@param x vector of sample values
#'@param y vector of sample numbers
#'@param variant should be equal 'ds' or 'mds'. if 'ds' then lambda is valid. on 'mds' prior is valid 
#'@param lambda for 'ds'
#'@param Max.m The maximum m partition to be used.
#'@param m.stats.wanted Boolean value indicating whether to return the list of Mm statistics.
#'@param prior vector of prior/penalty function to be used
.hhg.univariate.ks.stat.max = function(x, y, variant = 'Mm', lambda = 1,Max.m=2,m.stats.wanted=F,prior=NA,DS_type=1){
  w.max = 0
  w.sum = 2
  tables.wanted = F
  perm.stats.wanted = F
  nr.threads = 1
  
  if (variant == 'ds') {
    test_type = .UV_KS_DS
  } else if (variant == 'Mm') {
    test_type = .UV_KS_MDS
  } else {
    stop('Unexpected variant specified')
  }

  nr.perm = 0
  is.sequential = F
  total.nr.tests = 1
  alpha.hyp = NULL
  alpha0 = NULL
  beta0 = NULL
  eps = NULL
  
  wald = .configure.wald.sequential(is.sequential, total.nr.tests, alpha.hyp, alpha0, beta0, eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  is_sequential = as.integer(is.sequential)  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  # y is passed as numbers in 0:(K - 1), sorted according to the order of x
  if (is.factor(y)) {
    y = as.numeric(levels(y))[y]
  }
  y = as.matrix(as.double(y), nrow = length(y), ncol = 1)
  y = as.matrix(y[order(x)])
  
  # Dx and Dy are not used
  Dx = 0
  Dy = 0
  
  w_sum = as.double(w.sum)
  w_max = as.double(w.max)
  
  #fixme:  this assumes that prior length is at least maxk in length
  prior_param=F
  prior_length_param=1
  if(!is.na(prior[1])){
    prior_param=  prior
    prior_length_param=length(prior)
  }
  else{
    prior_param=rep(0,Max.m)  
    prior_length_param = Max.m
  }
  if(prior_length_param <Max.m-2){
    stop('prior supplied is not of sufficient length')
  }
  extra_params = c(as.double(DS_type), as.double(lambda),as.double(Max.m+0.1),as.double(prior_length_param+0.1),as.double(prior_param)) #0.1 is in order to make sure after int casting we get the right number
  
  res = .Call('HHG_R_C', test_type, Dx, Dy, y, w_sum, w_max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(y), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0, extra.stats.wanted = F,m.stats.wanted =m.stats.wanted ,m.stats.size = 2*(Max.m-1)) #changed size of n from nrow(Dx) to y
  
  return (ret)
}

# The k-sample version of the XDP test
.xdp.test.k.sample.mk = function(x, y, ddp.K = 3, w.sum = 0, w.max = 2) 
{
  # The interface may need more work..
  
  if (!is.numeric(x) && !is.ordered(x)) {
    stop('x is expected to be a numeric or ordered vector')
  }
  if (!is.numeric(y) && !is.factor(y) && !is.logical(y)) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  if (is.factor(y)) {
    y = as.numeric(levels(y))[y]
  } else if (is.logical(y)) {
    y = as.numeric(y)
  }
  if (any(y != round(y))) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  if (!is.vector(y) || length(x) != length(y)) {
    stop('x is expected to be a vector, and must have the same length as the vector y')
  }
  if (ddp.K != as.integer(ddp.K) || ddp.K < 2 || ddp.K > length(x)) {
    stop('K must be an integer between 2 and length(x)')
  } 
  
  test_type = .UV_KS_XDP_MK
  
  
  # Dx is used to store ranks of x (a permutation of 1:n)
  Dx = as.matrix(as.double(rank(x, ties.method = 'random')), nrow = length(x), ncol = 1)
  
  # y is passed as numbers in 0:(K - 1)
  y = as.matrix(as.double(y), nrow = length(y), ncol = 1)
  
  # Dy is not used
  Dy = 0
  
  w_sum = as.double(w.sum)
  w_max = as.double(w.max)
  
  # Can make these parameter at some point
  nr.perm = 0
  total.nr.tests = 1
  is.sequential = F
  alpha.hyp = NULL
  alpha0 = NULL
  beta0 = NULL
  eps = NULL
  nr.threads = 1
  tables.wanted = F
  perm.stats.wanted = F  
  
  extra_params = as.double(c(ddp.K))
  is_sequential = as.integer(is.sequential)
  
  wald = .configure.wald.sequential(is.sequential, total.nr.tests, alpha.hyp, alpha0, beta0, eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  res = .Call('HHG_R_C', test_type, Dx, Dy, y, w_sum, w_max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0,m.stats.wanted = T,m.stats.size=ddp.K*2)
  
  if (ddp.K > 3) {
    ret$max.chisq = NA
    ret$max.lr    = NA
    
    if (!is.null(ret$perm.pval.hhg.mc)) {
      ret$perm.pval.hhg.mc = NA
      ret$perm.pval.hhg.ml = NA
    }
  }
  ret_mat=matrix(NA,nrow = ddp.K-1,ncol = 2)
  ret_mat[,1]=ret$m.stats[1:(ddp.K-1)]
  ret_mat[,2]=ret$m.stats[(ddp.K-1) +1:(ddp.K-1)]
  colnames(ret_mat)=c('SC','SL')
  rownames(ret_mat)=paste('m=',as.character(2:(ddp.K)),sep = '')
  return (ret_mat)
}

# The k-sample version of the XDP test
.xdp.test.k.sample = function(x, y, ddp.K = 3, w.sum = 0, w.max = 2) 
{
  # The interface may need more work..
  
  if (!is.numeric(x) && !is.ordered(x)) {
    stop('x is expected to be a numeric or ordered vector')
  }
  if (!is.numeric(y) && !is.factor(y) && !is.logical(y)) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  if (is.factor(y)) {
    y = as.numeric(levels(y))[y]
  } else if (is.logical(y)) {
    y = as.numeric(y)
  }
  if (any(y != round(y))) {
    stop('y is expected to be a numeric, logical, or factor vector with values in {0, 1, ... K-1}')
  }
  if (!is.vector(y) || length(x) != length(y)) {
    stop('x is expected to be a vector, and must have the same length as the vector y')
  }
  if (ddp.K != as.integer(ddp.K) || ddp.K < 2 || ddp.K > length(x)) {
    stop('K must be an integer between 2 and length(x)')
  } 
  
  if (ddp.K == 2) {
    test_type = .UV_KS_XDP2
  } else if (ddp.K == 3) {
    test_type = .UV_KS_XDP3
  } else {
    test_type = .UV_KS_XDP
  }
  
  # Dx is used to store ranks of x (a permutation of 1:n)
  Dx = as.matrix(as.double(rank(x, ties.method = 'random')), nrow = length(x), ncol = 1)
  
  # y is passed as numbers in 0:(K - 1)
  y = as.matrix(as.double(y), nrow = length(y), ncol = 1)
  
  # Dy is not used
  Dy = 0
  
  w_sum = as.double(w.sum)
  w_max = as.double(w.max)
  
  # Can make these parameter at some point
  nr.perm = 0
  total.nr.tests = 1
  is.sequential = F
  alpha.hyp = NULL
  alpha0 = NULL
  beta0 = NULL
  eps = NULL
  nr.threads = 1
  tables.wanted = F
  perm.stats.wanted = F  
  
  extra_params = as.double(c(ddp.K))
  is_sequential = as.integer(is.sequential)
  
  wald = .configure.wald.sequential(is.sequential, total.nr.tests, alpha.hyp, alpha0, beta0, eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  res = .Call('HHG_R_C', test_type, Dx, Dy, y, w_sum, w_max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0)

  
  if (ddp.K > 3) {
    ret$max.chisq = NA
    ret$max.lr    = NA
    
    if (!is.null(ret$perm.pval.hhg.mc)) {
      ret$perm.pval.hhg.mc = NA
      ret$perm.pval.hhg.ml = NA
    }
  }
  
  return (ret)
}

#function for computing ADP for multiple partitions at once
.hhg.test.adp_mk = function(x, y, w.sum = 0, w.max = 2,
                            nr.perm = 0, M = 3,L=M, correct.mi.bias = F, total.nr.tests = 1, 
                            is.sequential = F, alpha.hyp = NULL, alpha0 = NULL, beta0 = NULL, eps = NULL, 
                            nr.threads = 1, tables.wanted = F, perm.stats.wanted = F)
{
  if (!is.vector(y)) {
    stop('y is expected to be a vector')
  }
  if (!is.vector(y) || length(x) != length(y)) {
    stop('x is expected to be a vector, and must have the same length as the vector y')
  }
  if (nr.perm < 0) {
    stop('nr.perm should not be negative')
  }
  if (!is.numeric(y) && !is.ordered(y)) {
    stop('y is expected to be a numeric or ordered vector')
  }
  if (!is.numeric(x) && !is.ordered(x)) {
    stop('x is expected to be a numeric or ordered vector')
  }
  
  test_type = .UV_IND_ADP_MK
  
  
  # Dx is used to store x
  Dx = as.matrix(as.double(rank(x, ties.method = 'random')), nrow = length(x), ncol = 1)
  y  = as.matrix(as.double(rank(y, ties.method = 'random')), nrow = length(y), ncol = 1)
  
  # Dy is not used
  Dy = 0
  
  # For historical reasons, the high-k DDP and ADP variants work on 1-based ranks, while the small-k
  # variants work on 0-based ranks.
  #if (!(variant == 'ddp.obs') && !(variant == 'ddp.all')) {
  #  Dx = Dx - 1
  #  y = y - 1
  #}
  
  w_sum = as.double(w.sum)
  w_max = as.double(w.max)
  
  extra_params = as.double(c(correct.mi.bias,length(M),M,L ))
  is_sequential = as.integer(is.sequential)
  
  wald = .configure.wald.sequential(is.sequential, total.nr.tests, alpha.hyp, alpha0, beta0, eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  res = .Call('HHG_R_C', test_type, Dx, Dy, y, w_sum, w_max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret.raw = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0,m.stats.wanted = T,m.stats.size = length(M)*2)
  ret=list()
  ret$m=M
  ret$l=L
  ret$sum.lr = ret.raw$m.stats[1:(length(ret.raw$m.stats)/2)+ (length(ret.raw$m.stats)/2)]
  ret$sum.chisq = ret.raw$m.stats[1:(length(ret.raw$m.stats)/2)]
  return (ret)
}


#general xdp independence wrapper
.hhg.test.udfree = function(x, y, variant = 'ppr.33.obs', w.sum = 0, w.max = 2,
                            nr.perm = 0, K = 3, correct.mi.bias = F, total.nr.tests = 1, 
                            is.sequential = F, alpha.hyp = NULL, alpha0 = NULL, beta0 = NULL, eps = NULL, 
                            nr.threads = 1, tables.wanted = F, perm.stats.wanted = F)
{
  if (!is.vector(y)) {
    stop('y is expected to be a vector')
  }
  if (!is.vector(y) || length(x) != length(y)) {
    stop('x is expected to be a vector, and must have the same length as the vector y')
  }
  if (nr.perm < 0) {
    stop('nr.perm should not be negative')
  }
  if (!is.numeric(y) && !is.ordered(y)) {
    stop('y is expected to be a numeric or ordered vector')
  }
  if (!is.numeric(x) && !is.ordered(x)) {
    stop('x is expected to be a numeric or ordered vector')
  }
  if (K < 2 || K > length(x)) {
    stop('m, number of partitions,  must be an integer greater than 1')
  }
  
  if (variant == 'spr.obs') {
    test_type = .UV_IND_DDP2
  } else if (variant == 'spr.all') {
    test_type = .UV_IND_ADP2
  } else if (variant == 'ppr.22.obs') {
    test_type = .UV_IND_DDP3_C
  } else if (variant == 'ppr.22.all') {
    test_type = .UV_IND_ADP3_C
  } else if (variant == 'ppr.33.obs') {
    test_type = .UV_IND_DDP3
  } else if (variant == 'ppr.33.all') {
    test_type = .UV_IND_ADP3
  } else if (variant == 'tpr.obs') {
    test_type = .UV_IND_DDP4
  } else if (variant == 'tpr.all') {
    test_type = .UV_IND_ADP4
  } else if (variant == 'ddp.obs') {
    test_type = .UV_IND_DDP
  } else if (variant == 'ddp.all') {
    test_type = .UV_IND_ADP
  } else {
    stop('Unexpected variant specified.')
  }
  
  # Dx is used to store x
  Dx = as.matrix(as.double(rank(x, ties.method = 'random')), nrow = length(x), ncol = 1)
  y  = as.matrix(as.double(rank(y, ties.method = 'random')), nrow = length(y), ncol = 1)
  
  # Dy is not used
  Dy = 0
  
  # For historical reasons, the high-k DDP and ADP variants work on 1-based ranks, while the small-k
  # variants work on 0-based ranks.
  if (!(variant == 'ddp.obs') && !(variant == 'ddp.all')) {
    Dx = Dx - 1
    y = y - 1
  }
  
  w_sum = as.double(w.sum)
  w_max = as.double(w.max)
  
  extra_params = as.double(c(K, correct.mi.bias))
  is_sequential = as.integer(is.sequential)
  
  wald = .configure.wald.sequential(is.sequential, total.nr.tests, alpha.hyp, alpha0, beta0, eps)
  alpha_hyp = as.double(wald$alpha.hyp)
  alpha0 = as.double(wald$alpha0)
  beta0 = as.double(wald$beta0)
  eps = as.double(wald$eps)
  
  nr_perm = as.integer(nr.perm)
  nr_threads = as.integer(nr.threads)
  tables_wanted = as.integer(tables.wanted)
  perm_stats_wanted = as.integer(perm.stats.wanted)
  
  res = .Call('HHG_R_C', test_type, Dx, Dy, y, w_sum, w_max, extra_params, is_sequential, alpha_hyp, alpha0, beta0, eps, nr_perm, nr_threads, tables_wanted, perm_stats_wanted)
  ret = .organize.results(res, n = nrow(Dx), nr.perm, tables.wanted, perm.stats.wanted, grid.len = 0)
  return (ret)
}

