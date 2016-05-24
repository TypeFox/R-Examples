##helper function that reduces imputed matrix to expected genotypes
.calcExpectedGeno<-function(genoData){
  apply(genoData,c(1,2),.mod,as.numeric(dimnames(genoData)[[3]]))
}


.allFreq <-
function(y){
  vecm = c(2,1,0) 
  (vecm[1:length(y)] %*% y) / (2*(sum(y)))
}
.calcOffset <-
function(phenv,family,pheno_cov){
  offsets1 = rep(NA,dim(pheno_cov)[2])
  incl = phenv[2,]==1
  dimres = dim(pheno_cov)[2]
  if(dimres>0) 
    offsets1[incl] =  glm(phenv[1,incl]~pheno_cov[incl,],family=family)$lin else offsets1[incl] = 0
  offsets1
}

.findIn<-function(header,str, critical=FALSE, grep=TRUE){
 for(k in 1:length(str)){
  index =  which(header==str[k])
  if(length(index)>0) break
 }
 index = index[!is.na(index)]
 if(grep & length(index)==0){
   for(k in 1:length(str)){
    index = grep(str[k], header) 
    if(length(index)>0) break
   } 
 }
 index = index[!is.na(index)]
 if(critical & length(index)==0){
   print(head(header))  
 stop(paste('could not match',paste(str,collapse=","), 'in header'))
 }
 if(length(index)>1) index = index[1]
 index
}

###MODIFIED 29-11
.readHeader<-function(header,type){
   formatIndex = .findIn(header,c("FORMAT","imputed"))
   if(length(formatIndex)==0) formatIndex = .findIn(header,c("FORMAT","regionId"))

   pos_index = .findIn(header,c("POS","FirstProbe","start"), critical=TRUE)
   chr_index = .findIn(header,c("CHROM","Chr","chr"), critical=FALSE)
   if(length(chr_index)==0) chr_index = -1
   rs_index = .findIn(header,c("snpid","rsid","regionId","ID","id"), critical=TRUE)
   lastIndex = .findIn(header,"LastProbe")
   NbSnpIndex =.findIn(header,"NbSnp")
   if(length(lastIndex)==0) lastIndex = pos_index
   if(length(NbSnpIndex)==0) NbSnpIndex=-1
   if(length(formatIndex)==0){
       firstGenoIndex = rs_index+1
       if(type=='impute')  firstGenoIndex =1+ .findIn(header,c("ALT"))
   }else  firstGenoIndex = formatIndex+1
   list(formatIndex = formatIndex,pos_index = pos_index,chr_index = chr_index, rs_index = rs_index,lastIndex = lastIndex,
	NbSnpIndex = NbSnpIndex, firstGenoIndex=firstGenoIndex,
	posi = c(pos_index,chr_index,rs_index,lastIndex,NbSnpIndex))
}



.centralise <-
function(vec,na.replace = NA){
  vec1 = (vec-mean(vec,na.rm=TRUE))
  vec1[is.na(vec1)] = na.replace
  vec1
}

.standardise <-
function(vec, na.replace=NA){
  vec1 = (vec-mean(vec,na.rm=TRUE))/sd(vec,na.rm=TRUE)
  vec1[is.na(vec1)] = na.replace
  vec1
}

.destandardise <-
function(vec,mean,sd) vec*sd + mean

.cnt <-
function(str) {
  if(length(grep("N",str))>0) 
    NA 
  else 
    nchar(gsub("A","",str))
}
.cntNa<-function(vec) length(which(is.na(vec)))

.cntVec <-
function(str) apply(as.matrix(str,nrow = length(str)),1,.cnt)
.convAvgToUnc <-
function(vec,n,sc=1){
  res = rep(0,n)
  fl = floor(vec)
  if(fl==vec) 
    res[vec+1] = sc
  if(fl!=vec){
    ceil = fl+1
    val = round(sc*(ceil-vec))
    res[fl+1] = val
    res[ceil+1] = sc-val
  }
  res
}
.convAvgToUncM <-
function(vec,sc,maximum){
 matr = apply(as.matrix(vec),1,.convAvgToUnc,maximum)
 as.vector(t(apply(as.matrix(vec),1,.convAvgToUnc,maximum,sc)))
}
.convAvgToUncMat <-
function(mat,sc){
 maximum = max(mat,na.rm=TRUE)+1
 if(maximum%%1>0) maximum = floor(maximum)+1
 apply(mat,2,.convAvgToUncM,sc,maximum)
}
.convertFactorToMatrix <-
function(f,nme){
   res= apply(as.matrix(levels(as.factor(f))),1,.getOnes,f)
   dimnames(res) = list(names(f),paste(nme,levels(as.factor(f)),sep="_"))
   res
}
.convertMatrixToFactor <-
function(mat){
 nme = as.matrix(data.frame(strsplit(dimnames(mat)[[2]],"_")))
 res = rep("",dim(mat)[1])
 for(i in 1:length(res)) 
   res[i] = paste(nme[2,mat[i,]==1],collase="")
 res
}
.estSigma <-
function(x){
  d = ncol(x)
  v = matrix(ncol=d,nrow=d)
  for( i in 1:d ){
    ptrI = !is.na(x[,i])
    v[i,i] = var(x[ptrI,i])
    if( i < d ){
      for( j in (i+1):d ){
        ptrJ = !is.na(x[,j])
        ptr = as.logical(ptrI*ptrJ)
        v[i,j] = cov( x[ptr,i], x[ptr,j] )
        v[j,i] = v[i,j]
      }
    }
  }
  return(v)
}
.expand <-
function(matr, repl){
 res = matr
 for(i in 2:repl) res = rbind(res,matr)
 res
}
.expand2 <-
function(matr, repl){
  res = matr
  for(i in 2:repl)
    res = abind(res,matr,along = 3)
  res
}
.expandCoeff <-
function(coeff,ind){
  diff =  max(ind) - dim(coeff)[1]
  if(diff>0) coeff = rbind(coeff,matrix(NA, nrow = diff,ncol = dim(coeff)[2]))
  coeff
}
.expandDat <-
function(data,repl){
  if(repl==1) {return(data)}

  inds = sapply(data,.nullDim)
  data[inds==2] = lapply(data[inds==2],.expand,repl)
  data[inds==3] = lapply(data[inds==3],.expand2,repl)
  data
}
.expDist <-
function(y){
   p = (y%*%2:0)/(2*sum(y))
   c(p^2, 2*p*(1-p), (1-p)^2)*sum(y)
}

.extractIds <-
function(vec,nme) {paste(nme[vec],collapse=",")[[1]]}
.extractIdsMat <-
function(mat,spl,caseInd) apply(mat,1,.extractIdsVec,spl,caseInd,dimnames(mat)[[3]])
.extractIdsVec <-
function(vec,spl,caseInd,nme) apply(.getin2(as.numeric(vec),spl,caseInd),2,.extractIds,nme)
.extractSig <-
function(res,ind = 2, totweight = 1,log.p=FALSE){
  vec = summary(res)$coeff[ind,]
  ppvs = if(log.p) log(2) + pt(abs(vec[3]),1,lower.tail=FALSE,log.p=TRUE) else 2*pt(abs(vec[3]),1,lower.tail=FALSE,log.p=FALSE)
  if(length(vec)==3) vec = c(vec, ppvs) 

  c(vec[c(1,4)],  sum(res$weights)/totweight)
}
.extractSig1 <-
function(res,ind=2,totweight=1,hasCov=FALSE,family="binomial",log.p=FALSE,useBF=FALSE){

  if(sum(is.na(res)) == 2){
    res1 = res
  }
  else {
    beta = summary(res)$coeff[ind,1]
   if(useBF){
        pv =  .abfSig(res,log.p=log.p)
        print(pv)
        if(hasCov){
           pv2 = .abfSig(update(res, ~covar,family=family) ,log.p=log.p)
           pv = if(log) pv-pv2 else pv/pv2
	  ###SHOULD CHECK THIS
        }
     }else{
    ll1 = logLik(res)

    nullfam = family=="ordinal"
    if(hasCov & !nullfam) 
      ll2 = logLik(update(res, ~covar,family=family)) 
    else 
      if(!nullfam) 
        ll2 = logLik(update(res, ~1,family=family)) 
      else 
        if(hasCov) ll2 = logLik(update(res, ~covar)) 
        else 
          ll2 = logLik(update(res, ~1))
    df = attr(ll1,"df")[1]-attr(ll2,"df")[1]
    pv = pchisq((2*(ll1 - ll2))/(totweight),df,lower.tail=FALSE,log.p=log.p)
     }
    res1 = c(beta,pv)
  }
  res1
}
.extractSig1P <-
function(res,ind=2,totweight=1,hasCov=FALSE,family="binomial",log.p=FALSE, bf = FALSE){
  if(sum(is.na(res)) == 2){
    res1 = res
  }
  else {
    ll1 = logLik(res)
    coeff = summary(res)$coeff
    beta = coeff[ind,1]
    ve = abs(as.numeric(coeff[ind,3]))
    pvs = if(log.p) log(2) + pt(ve,1,lower.tail=FALSE,log.p=TRUE) else 2*pt(ve,1,lower.tail=FALSE,log.p=FALSE)
    nullfam = family=="ordinal"
    if(!nullfam){
     if(hasCov)       ll2 = logLik(update(res, ~covar,family=family)) 
     else   ll2 = logLik(update(res, ~1,family=family))  
    }  
    else{
	 if(hasCov)	        ll2 = logLik(update(res, ~covar)) 
         else       ll2 = logLik(update(res, ~1))
    }
    pv = pchisq((2*(ll1 - ll2))/(totweight),attr(ll1,"df")[1]-attr(ll2,"df")[1],lower.tail=FALSE,log.p=log.p)
    res1 = cbind(c(beta,NA),c(pvs,pv))
  }
  res1
}
.extractSigExact <-
function(t){
  exactMethod = "wald"
  or1 = fisher.test(t)
  or = oddsratio(t,method = exactMethod)
  c(log(or$measure[2,1]),or1$p.value, sum(t))
}
.extractSigP <-
function(res,ind,totweight = 1,hasCov=FALSE,family="binomial",log.p = FALSE){
  #nobs= sum(res$weights)/totweight
  coeff = .expandCoeff(summary(res)$coeff,ind)
  vec  = coeff[ind,,drop=FALSE]
  ll1 = logLik(res)
#  print(ind)
#  print(vec)
#  print(vec[ind,])
  ppvs = if(log.p) log(2) + pt(abs(vec[,3]),1,lower.tail=FALSE,log.p=TRUE) else 2*pt(abs(vec[,3]),1,lower.tail=FALSE,log.p=FALSE)
  if(dim(vec)[2]==3){
    vec = cbind(vec, ppvs)
  }else{
    vec[,4] = ppvs
  }
  if(hasCov) 
    ll2 = logLik(update(res, ~covar,family=family)) 
  else 
    ll2 = logLik(update(res, ~1,family=family))
  pv = pchisq((2*(ll1 - ll2))/(totweight),attr(ll1,"df")[1]-attr(ll2,"df")[1],lower.tail=FALSE,log.p=log.p)
  res1 = cbind(c(vec[,1],NA),c(vec[,4],pv))

   #cbind(res1,rep(nobs,dim(res1)[[1]]))
  res1
}

.pseudoInv<-function(P){
svd = svd(P)
U = svd$u
d = svd$d
V = svd$u
d1 = diag(1/svd$d)
return(V %*% (d1 %*% t(U)))
}


.fillMissing <-
function( x){
  ptr.ind = which( apply(is.na(x),1,sum)>0 )
  if(length(ptr.ind)==0) return (x)

  v = .estSigma(x)
  m = apply(x,2,mean,na.rm=TRUE)
  lambda = .pseudoInv(v)
  xx = x
  for( ii in 1:length(ptr.ind) ){
    i = ptr.ind[ii]
    ptr.miss = which(is.na(x[i,]))
    fill =  m[ptr.miss] - .pseudoInv(lambda[ptr.miss,ptr.miss]) %*% lambda[ptr.miss,-ptr.miss] %*% as.matrix(x[i,-ptr.miss] - m[-ptr.miss])
    xx[i,ptr.miss] = fill
  }
  return(xx)
 
}
.findInd <-
function(vec,header,ind){x = grep(vec[ind],header)[1]; x}
.findmax <-
function(vec) max(0,max(as.numeric(vec),na.rm=TRUE))
.fisherPvalue <-
function(p){ 
	if(getOption("mPhen.log10p",FALSE)) pchisq(-2*log(10)*sum(p),df = 2*length(p),lower.tail=FALSE,log.p=TRUE)/log(10)
	else  pchisq(-2*sum(log(p)),df = 2*length(p),lower.tail=FALSE)
}
.fixGenoCols <-
function(geno, gCols){
 geno1 = matrix(0,nrow=dim(geno)[1],ncol=dim(geno)[2])
 geno1[,gCols] = apply(geno[,gCols,drop=FALSE],2,.cntVec)
 geno1[,-gCols] = apply(geno[,-gCols,drop=FALSE],2,as.numeric)
 geno1
}
.fixNonNumeric <-
function(pheno){
 nonnumvec = .nonNumeric1(pheno)
 pheno1 = matrix(0,nrow=dim(pheno)[1],ncol=dim(pheno)[2])
 pheno1[,!nonnumvec] = apply(as.matrix(pheno[,!nonnumvec]),2,as.numeric)
 pheno1[,nonnumvec] = apply(as.matrix(pheno[,nonnumvec]),2,.makeNumeric)
 dimnames(pheno1) = dimnames(pheno)
 pheno = pheno1
}
.getArray <-
function(arraynames,v=NA) array(v,dim= .getDim(arraynames), dimnames = arraynames)
.getDim <-
function(arraynames) lapply(arraynames,length)
#c(length(arraynames[[1]]), length(arraynames[[2]]), length(arraynames[[3]]))

.getFamily <-
function(phenv,incl){
 if(length(levels(as.factor(phenv[incl])))>2) "gaussian" else "binomial" 
}

#.getGlmNoCovar <-
#function(phenvec,vars, fam1,weights,offset){
#
#}
#.getGlmRev <-
#function(vars,phenvec,covar,weights, family) .glm1(vars ~ phenvec + covar,family=family,weights=weights)

#.getGlmRevNoCov <-
#function( vars, phenvec, weights, family){
#  .glm1(vars ~ phenvec,family=family,weights=weights)
#}

.abfSig<-function(res, w = .21 * .21,log.p=TRUE){
  fit <- summary(res)$coefficients[2,]
  v <- fit[2]
  z <- fit[3]
  labf <- 0.5*(log(v) - log(w+v)) + (0.5*z*z*w/(v+w))
  if(log.p) return(labf) else return(exp(labf))
}

##log bayes factor
.abf <- function( y, x, w=.21*.21 ){
  .abfSig(glm( y ~ x ,family=gaussian()))
}

.length<-function(x){
 if(is.null(x) )return(0)
 if(is.array(x)) return(dim(x)[2]) else return(length(x))
}

#log bayes factor
.mabf <- function(y1,x,delta =1, tau = 1,log.p=TRUE){
  n = length(y1)
  y = matrix(.standardise(y1),ncol=1,nrow =n)
  p = dim(x)[2]
  if(p==0) return(sum(dnorm(y1,log=T)))
  x = apply(x,2,.standardise)
  k = p 
  numerator = lgamma((n+delta+k)/2)
  denominator1 = (n/2)*log(pi) +  ((n-k)/2)*log(tau) +  lgamma((delta+k)/2)
  Xg = x
  Mg = diag(rep(tau,k)) + t(Xg) %*% Xg
  qg = t(y) %*% y -  t(y) %*% Xg%*% solve(Mg) %*% t(Xg) %*% y
  denominator2 =  (1/2)*log(det(Mg)) +  ((n+delta+k)/2)*log(1+qg/tau) 
  labf = numerator - denominator1 - denominator2
  if(log.p) return(labf) else return(exp(log.p))
}

##logP1 comes from model with more variables (model 1)
## returns posterior probability of model 1 versus model 2
## pi is prior odds (of bigger model)
.calcPPA<-function(logP1,logP2,totweight,pi= 0.05, log.p=TRUE){

logBF = (logP1-logP2)/totweight
PO = exp(logBF + log((pi/(1-pi))))
pv = PO/(1+PO)
if(log.p) pv = log(pv)
pv
}

.lrt<-function(ll1,ll2, totweight){
 df1 =  attr(ll1,"df")[1]
 df2 = attr(ll2,"df")[1]
 if(df2>df1) stop(paste("df2 should be bigger than df1 ",df2,df1))
 pchisq((2*(ll1 - ll2))/(totweight),df1-df2,lower.tail=FALSE,log.p=TRUE)
}

.getGlm <-
function(phenvec,vars, fam1,weights,covar=NULL,offset=NULL){

  if(is.null(covar)) .glm1(phenvec~vars,weights = weights,offset = offset, family = fam1)
  else {
    .glm1(phenvec~vars+covar,family = fam1,weights = weights, offset = offset)
  }
}

.glm1<-function(string, family, weights, offset=NULL){
  resu = NULL
  if(family[1]=="ordinal"){
    if(!is.null(offset)) stop('error')

    resu =  polr(string,Hess=TRUE,method="logistic",weights = weights)

  }else{
    resu = glm(string,family=family, weights = weights,offset = offset)
  }
  if(family=="gaussian") attr(resu,"totweight") = 1
  else attr(resu,"totweight") = mean(weights)
  resu
}


.backwardSelection<-function(vars,phenvec,covar,weights, family, totweight, iterationsMax = 100,minVarLength =1){
hascov=!is.null(covar)
 logthresh = log(getOption("mPhen.defaultBackwardThresh",default=0.01))
 frac1 =getOption("mPhen.defaultBackwardFrac", default = 0.1)
 bf = getOption("mPhen.defaultUseBF",default = FALSE)

if(family=="gaussian") totweight=1
else totweight = mean(weights)
 res = .getGlm(vars,phenvec,family,weights, covar)
 res2 = .getGlm(vars,rep(1, dim(phenvec)[1]),family,weights, covar)
 ll1 = logLik(res)
 ll2 = logLik(res2)

 if(bf){
	abf1 = .mabf(vars,phenvec)  ## should add covar
 }
 pv = .lrt(ll1,ll2,totweight)


 noP = dim(phenvec)[2]
 torem = c()
 pvl = rep(1.0,noP)

 inds = 1:noP
 toremnew = c()

 endloop<-FALSE
 cnt =1

 while(!endloop ){
   frac = frac1[1]
   pvl = rep(0,length(inds))
  # pvlg = rep(0,length(inds))
  # abf = rep(0,length(inds))
  
   for(i in 1:length(inds)){
      if(bf){
	    abfi = .mabf(vars,phenvec[,-c(torem,inds[i]),drop=F])
           pvl[i] =.calcPPA(abfi,abf1,totweight,pi= 0.5)
      }
      else{
          
            
          if(length(inds[-i])==0) lli =ll2
          else{
          
             res2 = .getGlm(vars, phenvec[,-c(torem,inds[i])], family, weights, covar)
             lli = logLik(res2)
	  
          }
	  pvl[i] = .lrt(ll1,lli,totweight)
      }
   }
   ordp = rev(order(pvl))
   names(pvl) = dimnames(phenvec)[[2]][inds]
   len = length(which(pvl[ordp]>logthresh))
   step = max(1,floor(frac*length(pvl)))
   if(len==0) toremnew = c() 
   else toremnew = inds[ordp[1:min(step,len,length(ordp)-minVarLength)]]


  if(cnt <=iterationsMax 
	& length(toremnew)>0
        & length(c(torem,toremnew)) <= noP - minVarLength 
	##& length((1:noP)[-c(torem,toremnew)])>=minVarLength
       ){
          torem = c(torem,toremnew)
          
           res = .getGlm(vars,phenvec[,-torem,drop=F],family,weights, covar)
	   #res = .glm1(vars ~ covar + phenvec[,-torem,drop=F],family=family,weights=weights)

          ll1 = logLik(res)

          if(bf){
           abf1 = .mabf(vars,phenvec[,-torem,drop=F])
	  }

	  pv_new = .lrt(ll1,ll2,totweight)
 	        inds = (1:noP)[-torem]
		cnt = cnt+1
	        pv = pv_new
   }else{
	endloop=TRUE
   }
   
 }
 summ = summary(res)
if(length(torem)==0) torem = c(noP+1)
#without intercept

lastind = 1
if(hascov)  lastind = (lastind+dim(covar)[[2]])
 diff =length(summ$alias)- dim(summ$coeff)[1] 

 coeff = summ$coeff[-(1:lastind),,drop=F]
 alias = summ$alias[-(1:(lastind+diff))]
 pv = .lrt(ll1,ll2,totweight)

 indSig = which(pvl<=logthresh)

 ind = (1:noP)
 beta = rep(0,noP)
 pvs = rep(1,noP) 
if(is.null(alias)){
 beta[-torem][indSig] = coeff[,1][indSig]
}else{
 beta[-torem][which(!alias)][indSig] = coeff[,1][indSig]
}
 pvs[-torem][indSig] = exp(pvl[indSig])
 res1 = cbind(c(beta,NA),c(pvs,exp(pv)))
  res1
}

#.getGlmRevOrd <-
#function( vars, phenvec, covar,weights) {
# tryCatch(polr(vars ~ phenvec + covar,Hess=TRUE,method="logistic",weights = weights), error = function(e) print(NULL))
#}
#.getGlmRevOrdNoCov <-
#function( vars, phenvec, weights){
#  tryCatch(polr(vars ~ phenvec,Hess=TRUE,method="logistic",weights = weights), error = function(e) print(NULL))
#}
.getHist <-
function(vec, spl){ 
  le = vec < spl[1]
  gt = vec > spl[length(spl)]
  c(length(which(le)),hist(as.numeric(vec[!le & !gt]),br=spl,plot=FALSE)$count,length(which(gt)))
}
.getIMatrix <-
function(d,phi,cats,num1,x){
  II = matrix(ncol=d,nrow=d)
  for(j in 1:d){
    for(k in j:d){
      num = 0
      for( i in 1:cats ) num = num + (i-1)*(i-1)*x[,j]*x[,k]*exp(phi[i])
      II[j,k] = sum(num1[[j]]*num1[[k]])
      II[j,k] = II[j,k] - sum(num)
      II[k,j] = II[j,k]
    }
  }
  II
}
.getin <-
function(stend,vec) vec>=stend[1] & vec<stend[2]
.getin1 <-
function(genvec,spl) apply(cbind(spl[1:(length(spl)-1)],spl[2:length(spl)]),1,.getin, genvec)
.getin2 <-
function(genvec,spl,cc) apply(cbind(spl[1:(length(spl)-1)],spl[2:length(spl)]),1,.getin, genvec) & cc
.getIncl <-
function(phenv, na_cov){incl = !((is.na(phenv) | na_cov)); incl}
.getMaf <-
function(matr,base,noallele){
index = !is.na(matr[1,1,]) & matr[1,2,]>0 
vec = matr[1,1,index]
weights = matr[1,2,index]
res=  .minOneMinus(sum(weights*abs(vec -base))/(sum(weights)*noallele))
res
}
.getHWE<-function(matr,totweight){
 index = !is.na(matr[1,1,]) & matr[1,2,]>0 
 vec = matr[1,1,index]
 weights = matr[1,2,index]
 vec1 = .countIn(vec,weights,-2.5:2.5)/totweight
# vec1 = round(vec1) 
min(.hwep(vec1[3:5]),.hwep(vec1[1:3]))
}

.countIn<-function(vec,weights,br){
 
 res = rep(NA, length(br)-2)
 for(i in 1:(length(br-1))){
  min = br[i]
  max = br[i+1]
  ind = vec>=min & vec < max
   res[i] = sum(weights[ind])
 }
res
}
.getOnes <-function(vec,f) f==vec[1]

.getPvLrr <-
function(phend, gvar, family,inclu,covar,totweight){
  log.p = getOption("mPhen.log10p",FALSE)
  phenvec = phend[1,]
  vars = gvar[1,]
  weights = gvar[2,]
  offset = phend[4,]
  include = phend[2,] & inclu & gvar[3,]>0 & weights > 0 & sd(vars,na.rm=T)>0
  if(dim(covar)[2]==0) covar = NULL   else covar = covar[include,]

  if(length(which(include))>0 & var(phenvec[include],na.rm=TRUE)>0){
      glm.res = .getGlm(phenvec[include] ,as.numeric(vars[include]),family,weights[include],covar = covar, offset = offset[include])
    totweight1 = attr(glm.res,"totweight")
    res =  .extractSig1(glm.res,ind=2,totweight = totweight1,hasCov = !is.null(covar),family=family,log.p=log.p)
  }
  else res = rep(NA,2) 
  c(res,sum(weights[include])/totweight)
}

.getPvLrrAllGen <-
function(vars,phen,families,include,covar,totweight,functiAll,functi,jointModel){
 dimv = dimnames(vars)[[1]]
 resn = dimnames(phen)[[1]]
 if(jointModel) resn = c(resn,"all")
  x = .getArray(list(dimv,resn,c("beta","pvalue","Nobs")))
  for(i in 1:length(dimv)){
        ress = functiAll(vars[i,,],phen,families,include,covar,totweight[i],functi)
	  x[i,,] = ress
  }
  x
}
.getPvLrrMult <-
function(vars,phen,fams,include,covar,totweight, functi){
  phenN = dimnames(phen)[[1]]
  res = .getArray(list(phenN,c("beta","pvalue","Nobs")))
  for(k in 1:length(phenN)){   
    res[k,] = functi(phen[k,,],vars, fams[k],include,covar,totweight)
  }
  res
}


.getPvLrrMultiGenVS<-
function(gvar,phend, family,inclu,covar,totweight,functiAll,functi,jointModel){
 .getPvLrrMultiGen(gvar,phend, family,inclu,covar,totweight,functiAll,functi,jointModel, var_sel=TRUE)
}
.getPvLrrMultiGenNoVS<-
function(gvar,phend, family,inclu,covar,totweight,functiAll,functi,jointModel){
.getPvLrrMultiGen(gvar,phend, family,inclu,covar,totweight,functiAll,functi,jointModel,var_sel=FALSE)
}

.getPvLrrMultiGen <-
function(gvar,phend, family,inclu,covar,totweight,functiAll,functi,jointModel, var_sel = FALSE){
  phenN = dimnames(phend)[[1]]
  res = .getArray(list(c(dimnames(gvar)[[1]],"combined"),phenN,c("beta","pvalue","Nobs")))
 for(k in 1:length(phenN)){
  res[,k,] = .getPvLrrMultiGen1(phend[k,,],gvar,family[k],inclu,covar,mean(totweight),functiAll,functi, var_sel = var_sel)
  }
  res
}

.ordTestUnivariate <-
function(phend,gvar, family, inclu,covar,totweight){
 if(dim(covar)[2]>0) stop('ord test not coded up for covariates, use residuals')
  weights = as.numeric(gvar[2,])
  include = phend[2,]>0 & inclu & weights > 0 & gvar[3,]>0
  y = as.numeric(gvar[1,include])
  x = (as.matrix(phend[1,include] - phend[4,include]))
  #print(y)
  #print(head(x))
  resu  = .ordTest2(y,x)
  resu = c(resu[1,], dim(covar)[1])
  resu
}




.getPvRev <-
function(phend, gvar, family, inclu,covar,totweight){ 
 log.p = getOption("mPhen.log10p",FALSE)
  vars = gvar[1,]
  weights = as.numeric(gvar[2,])
  include = phend[2,]>0 & inclu & weights > 0 & gvar[3,]>0

  phenvec = phend[1,include] - phend[4,include]
  emptyCov = dim(covar)[2]==0
  if(emptyCov) covar = NULL else covar = covar[include,]
  var1 = vars[include]
  cont = (family=="gaussian")
  binom = (family=="binomial")
  ord = (family=="ordinal")
  if(ord){
     var1 = as.factor(var1)
   if(length(levels(var1))<=2){
    binom = TRUE
    family = "binomial"
    ord = FALSE
  }
  }
  todo = TRUE
  if(length(which(include))==0 | var(var1,na.rm=TRUE)==0) 
    todo=FALSE
  if(todo){


    glm.res = tryCatch(.getGlm(var1,phenvec,family, weights[include], covar = covar,offset = NULL),error = function(e) print(NULL))
      if(is.null(glm.res)){
	family = "gaussian"
        print("polr failed")
	glm.res = .getGlm(as.numeric(as.character(var1)),phenvec,family, weights[include],covar = covar, offset = NULL)
      }

       totweight1 = attr(glm.res,"totweight")
     c(.extractSig1(glm.res,ind=1, totweight = totweight1,hasCov=!emptyCov,family=family,log.p=log.p),sum(weights[include]/totweight))
 }
  else if(!todo) rep(NA,3) 
 }

.getPvLrrMultiGen1 <-
function(phend, gvar,family,inclu,covar,totweight,functiAll,functi, var_sel = FALSE){
 log.p = getOption("mPhen.log10p",FALSE)
  phenvec = phend[1,]
  vars = gvar[,1,]
  weights1 = gvar[1,2,] 
  nainds =  apply(as.matrix(gvar[,3,]),2,min)
  offset = phend[4,]
  include = phend[2,] & inclu & nainds >0 & weights1 > 0
  # & nainds # & weights > 0
  emptyCov = dim(covar)[2]==0
  if(emptyCov) covar = NULL else covar = covar[include,,drop=F]
  ind = 1:(dim(vars)[1]) + 1
  tobind = sum(weights1[include])/totweight   #,length(families)+1)
  v = var(phenvec[include])
  if(length(which(include))>0  & !is.na(v) & v>0){
    if(var_sel){
      res = .backwardSelection(phenvec[include]-offset[include],t(vars[,include]),covar=covar,weights=weights1[include], family = family, totweight=totweight)
      res = cbind(res,tobind)
     return(res)
    }
    glm.res = .getGlm(phenvec[include],t(vars[,include]),family,weights1[include],covar = covar, offset = offset[include]) 
     totweight1 = attr(glm.res,"totweight")
    cbind(.extractSigP(glm.res,ind=ind,totweight=totweight1,family=family,hasCov=!emptyCov,log.p=log.p),sum(weights1[include])/totweight)
  }
  else rep(NA,3)
}



###uses the first family
.getPvRevPleio_ <-
function(gvar,phen,families,inclu,covar,totweight,functi, variable_selection){
 #print(families)
  log.p = getOption("mPhen.log10p",FALSE)
  phenvec = phen[,1,]-phen[,4,]
  family = families[1]
   vars = gvar[1,]
  weights = gvar[2,]
  nainds  =   gvar[3,] 
  phenincl = apply(as.matrix(phen[,2,]),2,min)
  include =  phenincl>0 & inclu & weights>0 & nainds >0
  var1 = vars[include]
  cont = (family=="gaussian")
  binom = (family=="binomial")
  ord = (family=="ordinal")
  if(ord){
     var1 = as.factor(var1)
   if(length(levels(var1))<=2){
    binom = TRUE
    family = "binomial"
    ord = FALSE
  }
  }
  emptyCov = dim(covar)[2]==0
  if(emptyCov)covar = NULL else covar = covar[include,,drop=F]
  ind = 1:(dim(phenvec)[1])
  if(!ord) ind = ind+1
  todo = TRUE
  tobind = rep(sum(weights[include])/totweight,length(families)+1)
  if(length(which(include))==0 | var(var1,na.rm=TRUE)==0) 
    todo=FALSE
  if(todo){
    if(variable_selection){
      res = .backwardSelection(var1,t(phenvec[,include]),covar=covar,weights=weights[include], family = family, totweight=totweight)
      res = cbind(res,tobind)
      return(res)
    }else{
      glm.res = NULL
      glm.res = tryCatch(.getGlm(var1,t(phenvec[,include]),family, weights[include],covar = covar, offset = NULL),error = function(e) print(NULL))
      if(is.null(glm.res)){
	family = "gaussian"
        print("polr failed, using Gaussian")
	glm.res = .getGlm(as.numeric(as.character(var1)),t(phenvec[,include]),family, weights[include],covar = covar, offset = NULL)
      }
      totweight1 = attr(glm.res,"totweight")
      resu= .extractSig1P(glm.res, ind = ind, totweight = totweight1, hasCov = !emptyCov, family = family,log.p = log.p)
      return(cbind(resu,tobind))
    }
  }  
  if(!todo) 
    return(rep(NA,3) )
}


.getPvRevPleio<-function(gvar,phen,families,inclu,covar,totweight,functi) .getPvRevPleio_(gvar,phen,families,inclu,covar,totweight,functi, FALSE)
.getPvRevPleioVS<-function(gvar,phen,families,inclu,covar,totweight,functi) .getPvRevPleio_(gvar,phen,families,inclu,covar,totweight,functi, TRUE)


.getPvTable <-
function(phend, gvar, family,inclu,covar,totweight){
  phenvec = phend[3,]
  vars = as.factor(gvar[1,])
  weights = gvar[2,]
  offset = phend[4,]
  include = phend[2,] & inclu & weights > 0 & gvar[3,]>0
  emptyCov = dim(covar)[2]==0
  todo = TRUE
  if(length(which(include))==0 | var(vars[include],na.rm=TRUE)==0 | family!="binomial") 
    todo=FALSE
  if(!todo) 
    rep(NA,3) 
  else 
    .extractSigExact(.tabulateWithWeights(phenvec[include],vars[include],weights[include],totweight))
}
.getUNum1Matrix <-
function(d,y,x,phi,num1,cats){
  U = matrix(ncol=1,nrow=d)
  num1 = list(length=d)
  for( j in 1:d ){
    U[j] = sum(y*x[,j])    
    num1[[j]]=0 
    for( i in 1:cats ) 
      num1[[j]] = num1[[j]] + x[,j] * (i-1) * exp(phi[i]) 
    U[j,1] = U[j,1] - sum(num1[[j]])
  }
  list(U,num1)
}
.hwestat <-
function(obs,exp) sum((obs-exp)^2/exp)
.iscase <-
function(vec) !is.na(vec) & vec>mean(vec,na.rm=TRUE)
.joinRes <-
function(mat,ncol,form){
 if(is.null(dim(mat))) 
   mat = matrix(mat,ncol=ncol)
 res1 =  apply(mat,2,.joinResVec,form)
 res1
}

.hwep<-function(vec){
 if(max(vec)==0){
   p = 1.0
 }
 else if(min(vec)<5){
   p = HWExact(vec)$pv
 }
 else{
 p=pchisq(.hwestat(vec, .expDist(vec)),1,lower.tail=F)
}
p
}
.hwepall<-function(y){
h=hist(y,plot=F,br=-2.5:2.5)
min(.hwep(h$counts[3:5]),.hwep(h$counts[1:3]))
}



.joinResVec <-
function(vec,form) paste(sprintf(form,vec),collapse=";")
.mafCount <-
function(y) max(.minOneMinus(.allFreq(y[1:min(3,length(y))])), .minOneMinus(.allFreq(y[3:min(5,length(y))])),na.rm=TRUE)
.makeFactor <-
function(mat1,nme){
  #if(dim(mat1)[2]>0){
  # for(k in 1:(dim(mat1)[2])){
  #   if(length(levels(as.factor(mat1[,k])))>3){
  #       mat1[,k] = .makeTopTail(mat1[,k],c(0.33,0.66))
  #   }
  # }
  #}
  append = dim(mat1)[2]>0
  fact = factor(apply(mat1,1,paste,collapse=".")) 
  lev = levels(fact) 

  matrix = .convertFactorToMatrix(fact,paste(nme,sep="."))

  nmes = dimnames(matrix)[[2]]
  txt = dimnames(matrix)[[2]]
  index = setdiff(seq(length(txt)),grep("_NA",txt))
  res = matrix[,index,drop=FALSE]
  dimnames(res) = list(NULL,nmes[index])
  if(append)  
    res = cbind(res, rep(TRUE,dim(res)[1]))
  if(append)   
    dimnames(res) = list(NULL,c(nmes[index],"all")) 
  res 
}
.makeNumeric <-
function(vec) as.numeric(as.factor(vec))



.makeQuantile <-
function(vec)  qqnorm(vec,plot=FALSE)$x

.makeThresh <-
function(vec,thresh){
 vec1 = rep(NA,length(vec))
 vec1[vec<=thresh[1]] = 0
 vec1[vec>=thresh[2]] = 1
 vec1
}
.makeTopTail <-
function(vec,perc1){
 perc = as.numeric(perc1)/100
 vec1 = rep(NA,length(vec))
 quants = ecdf(sort(vec))(vec)
 vec1[quants<=perc[1]] = 0
 vec1[quants>=perc[2]] = 1
 vec1
}

 ##HAS BEEN MODIFIED 29-11


.mergeAlleleCols <-
function(snp,cols){
  len = 1+ dim(snp)[2] - length(cols) 
  snp1 = matrix(0,nrow = dim(snp)[1],ncol =len) 
  snp1[,1] = .cntVec(apply(snp[,cols],1,paste,collapse=""))
  matr = snp[,-cols,drop=FALSE]
  nme1 = rep("",len)
  nme1[1] = "geno"
  if(len>1) 
    nme1[2:len] = dimnames(snp)[[2]][-cols]
  if(len>1) 
    for(i in 2:len) 
    snp1[,i] = matr[,i-1]
  dimnames(snp1) = list(dimnames(snp)[[1]],nme1)
  snp1
}
.metaresFisher <-
function(matr){
  inds = !is.na(matr[,2] )
  fisherp = .fisherPvalue(matr[inds,2])
  return(c(mean(matr[inds,1]),fisherp,sum(matr[,3])))
}
.metaresInvVarianceFixed <-
function(matr){
  inds = which(!is.na(matr[,2]))
  teSe = apply(matr[inds,,drop=FALSE],1,.seFromPBeta)
  if(length(inds)<getOption("mPhen.metaMinCohortNum",2)) return(c(NA,NA,NA))
  else if(length(inds)==1) return(matr[inds,])
  te = matr[inds,1]
  if(getOption("mPhen.useFisherIfAllBetaNA",FALSE) & .cntNa(te)==length(te)){
    return(.metaresFisher(matr))
  }
  metres = metagen(te,teSe)
  pv = .pFromSEBeta(c(metres$TE.fixed,metres$seTE.fixed))
  ##c(metres$Q,pchisq(metres$Q,1,lower.tail=FALSE)) 
  c(metres$TE.fixed,pv,sum(matr[,3]))
}

.metaresHetMeasure<-function(matr){

  inds = which(!is.na(matr[,2]))
  if(length(inds)<getOption("mPhen.metaMinCohortNum",2)) return(c(NA,NA,NA))
  teSe = apply(matr[inds,,drop=FALSE],1,.seFromPBeta)
  te = matr[inds,1,drop=FALSE]
  metres = metagen(te,teSe)
  pv = .pFromSEBeta(c(metres$TE.fixed,metres$seTE.fixed))
  c(metres$Q,pchisq(metres$Q,1,lower.tail=FALSE),sum(matr[,3]))
}
.pFromSEBeta <-function(v){
  if(getOption("mPhen.log10p",FALSE)){
       x =  (pnorm(abs(v[1]),sd = v[2],lower.tail=FALSE,log.p=TRUE)+log(2))/log(10);
   }else{
	x =  pnorm(abs(v[1]),sd = v[2],lower.tail=FALSE)*2;
   }
   x
}
.seFromPBeta <-function(v){
  if(getOption("mPhen.log10p",FALSE)){
    x = abs(v[1])/qnorm((v[2]*log(10))/2,lower.tail=FALSE,log.p=TRUE); 
  }else{
   x = abs(v[1])/qnorm(v[2]/2,lower.tail=FALSE); 
  }
x
}


.getSuffix<-function(todo,sub=""){
  todo_ = todo[grep(sub,todo[,1]),,drop=F]
    l=levels(as.factor(todo_[,1]))
  if(length(l)>1){

    res = .getSuffix(todo_,sub=l[1])
    for(k in 2:length(l)){
      res = paste(res,"_",.getSuffix(todo_,sub=l[k]),sep="")
    }
    return(res)
  }else{
   l2=levels(as.factor(todo_[,2]))
   if(length(l2)>2){
     return(paste(sub,"_",length(l2),sep=""))
   }else{
     return(paste(sub,"_",paste(l2,collapse="."),sep=""))
}

}
}

.fixDuplicatesInTodo<-function(todo1_){
ph_i = grep("pheno",todo1_[,1])
todo_ = rbind(todo1_[-ph_i,,drop=F],todo1_[ph_i,,drop=F])
type = todo_[,2]
a = levels(as.factor(type))
toexcl = c()
for(i in 1:length(a)){
  match = which(type==a[i])
  if(length(match)>1){
    toexcl = c(toexcl,match[2:(length(match))])
  }
}
if(length(toexcl)>0){
return(todo_[-toexcl,,drop=F])
}else{
return(todo_)
}
}
.minOneMinus <-
function(x) min(x,1-x)
.mod <-
function(vec,vals){
  #vec[vec==""] = 0; 
  #vec[vec=="-"] = 1000 - sum(as.numeric(vec[vec!="-"]),na.rm=TRUE); 
  (as.numeric(vec)%*% vals)/sum(vec)
}
.mod1 <-
function(vec, rescale){
  #vec[vec==""] = 0
  #vec[vec=="-"] = 1000 - sum(as.numeric(vec[vec!="-"]),na.rm=TRUE);
  #vec[vec<200]  =0
  vec1 = floor(rescale * as.numeric(vec))#floor(rescale*(as.numeric(vec)/1000))
  ind = vec1==max(vec1)
  vec1[ind] = vec1[ind]+ rescale - sum(vec1)
  vec1
}

#this function expands a genotype into a distribution over genotypes (just used for testing)
.expandGenoToDist<-function(geno,vals,scale,perc){
	res=rep(perc*scale,length(vals));
       
	res[which(vals==geno)] =res[which(vals==geno)]+scale*(1-perc)
	res
}
.expandGenoToDistAll<-function(geno,spl,scale, perc=0.0){
vals = floor(min(spl)):(floor(max(spl))+1)
vals[[1]]
res = apply(geno,c(1,2),.expandGenoToDist,vals,scale,perc)
res1=aperm(res,c(2,3,1))
dimnames(res1)[[3]] = vals
res1
}

.reduceWeightsToGeno<-function(geno){
geno=data$geno
inds=which(geno[1,2,]>0)
geno1 = geno[1,,inds]
t(geno1[1,,drop=F])
}

.applyTrans<-function(phenoData,limit){
pheno1 = phenoData[,unique(match(limit[,2],dimnames(phenoData)[[2]])),drop=F]


  if(dim(limit)[2]>2)phenoTrans = limit[,3] else phenoTrans = NULL
if(!is.null(phenoTrans)){
    qqfam =  grep("^quantile",phenoTrans)
    if(length(qqfam)>0) pheno1[,qqfam] = apply(pheno1[,qqfam,drop=FALSE],2,.makeQuantile)
    stdfam =  grep("^standard",phenoTrans)
    if(length(stdfam)>0) pheno1[,stdfam] = apply(pheno1[,stdfam,drop=FALSE],2,.standardise)

    toptail = grep("^toptail",phenoTrans)
    if(length(toptail)>0) pheno1[,toptail] = apply(pheno1[,toptail,drop=FALSE],2,.makeTopTail,strsplit(phenoTrans[toptail[1]],"_")[[1]][2])
    thresh = grep("^thresh",phenoTrans)
    if(length(thresh)>0) pheno1[,thresh] = apply(pheno1[,thresh,drop=FALSE],2,.makeThresh,strsplit(phenoTrans[thresh[1]],"_")[[1]][1:2])
    #dimnames(pheno1)[[2]] = dimnames(phenoData)[[2]]
    factor = grep("^factor",phenoTrans)[1]
    while(length(factor)>0 & !is.na(factor[1])){
       k= factor[1]
       newM =  .makeFactor(as.matrix(pheno1[,k]),limit[k,2])
       newM = newM[,-dim(newM)[[2]]]
       pheno1 = cbind(pheno1[,-k,drop=F],newM)
       newLim = cbind(limit[k,1],dimnames(newM)[[2]],"")
        limit = rbind(limit[-k,,drop=F],newLim)
            phenoTrans = limit[,3]
       factor = grep("^factor",phenoTrans)[1]
   }
  }
#    if(!is.null(phenoTrans)) phenN = apply(rbind(phenN,phenoTrans),2,paste,collapse=".")
list(pheno1=pheno1,limit=limit)
}


.applySidakCorrection<-function(pvs,effPhe){
  pv = min(pvs,na.rm=T)
if(getOption("mPhen.log10p",FALSE)){
   if(pv<1e-15) res = pv+log10(effPhe)
   else res = 1-(1-(10^pv))^effPhe
}else{
  if(pv<1e-15) res = pv*effPhe
  else res = 1-(1-pv)^effPhe
}
res
}
.nyholdtSidak<-function(pheno){
M = dim(pheno)[[2]]
if(M==1){
  return(1)
}else{
cor = cor(pheno,use="pairwise.complete.obs")
excl_ind = (apply(apply(cor,c(1,2),is.na),1,max)>0)
v =var(eigen(cor[!excl_ind,!excl_ind])$val)
return((M-1)*(1-v/M)+1)
}
}

.readPhen<-function(phenFile,sep="\t",numHeaderRows = 1){
  if(length(grep("zip",phenFile))>0) phen3 = .readZip(phenFile)
  else{
    phen2 = read.table(phenFile,header=numHeaderRows>0,as.is=T,sep=sep,comment.char="")
    inds = match(unique(phen2[,1]),phen2[,1])
     phen3 = phen2[inds,-1,drop=F]
     rown =     phen2[inds,1] 
    if(numHeaderRows==0){
	coln = paste("phen",1:(dim(phen2)[2]-1),sep="")
       dimnames(phen3) = list(phen2[inds,1],coln)
    }else{
      dimnames(phen3)[[1]] = rown
    }
  }
  if(numHeaderRows>1){
    toincl = (numHeaderRows:(dim(phen3)[1]))
     phen3 = phen3[toincl,,drop=F]
  }
  p4 = as.matrix(apply(phen3,2,.as.numeric1))
  dimnames(p4) = dimnames(phen3)

  return(p4)
}
.as.numeric1<-function(v){
 f = as.factor(v)
 le = levels(f)
 if(length(le)<=3){
     return(as.numeric(f))
 }else{
    return(as.numeric(v))
 }
}
.getLast<-function(vec) vec[length(vec)]

.findInclusion<-function(sampleids, match, any=FALSE){
  if(any){
   rn = c()
   for(k in 1:length(sampleids)){
      rn = c(rn,sampleids[[k]])
   } 
   rn = unique(rn)
   pos = as.numeric(unlist(lapply(strsplit(rn,"_"),.getLast)))
   ord = order(pos)
   rn = rn[ord]
   }else  rn = sampleids[[1]]
  inclrow = list()
  for(i in 1:length(sampleids)){
        if(is.null(match)) inclrow[[i]] = match(rn,sampleids[[i]])
        else if(match) inclrow[[i]] = which(!is.na(match(rn,sampleids[[i]])))
        else{
            if(i==1) inclrow[[i]] = which(!is.na(match(rn,sampleids[[i]])))
            else inclrow[[i]] = which(is.na(match(sampleids[[i]],rn)))
        }
  }
  attr(inclrow,"rn")<-rn
  inclrow
}

.extractNames<-function(phenos,ind){
  sampleids = list()
  for(k in 1:length(phenos)){
    sampleids[[k]] = dimnames(phenos[[k]])[[ind]]
  }
  sampleids
}

.trL<-function(phenos){
res = list()
for(k in 1:length(phenos)){
  res[[k]] = t(phenos[[k]])
  #dimnames(res[[k]]) = rev(dimnames(phenos[[k]]))
}
res
}

###MOD 29-11
.mergeFiles<-function(phenos, 
      markerCol=FALSE,
      rownames = .extractNames(phenos,1),
      colnames = .extractNames(phenos,2),
      inclrow = NULL,
      ind2 = NULL, anyCol=FALSE	      
){
  if(length(phenos)==1) return(phenos[[1]])
  if(is.null(inclrow))       inclrow = .findInclusion(rownames,match=FALSE)
  if(is.null(ind2)) ind2 = .findInclusion(colnames,match=NULL, any=anyCol)
  rn =  attr(ind2,"rn")
  phen1 = array(dim = c(0,length(rn)))
  Cohort_ = c()
  indiv = list()
  for(i in 1:length(phenos)){
     phen2 = phenos[[i]][inclrow[[i]],ind2[[i]],drop=F]
    phen1 = rbind(phen1,phen2)
    indiv[[i]] = rownames[[i]][inclrow[[i]]]
    Cohort_ = c(Cohort_,rep(i,dim(phen2)[[1]]))
 }
 dimnames(phen1)[[2]] =rn
 if(markerCol & length(phenos)>1) phen1 = cbind(phen1,Cohort_)

 if(length(indiv)==1) indiv = indiv[[1]]
 attr(phen1,"indiv")<-indiv
 phen1
}




.rescale<-function(genoData,.rescaleV, integer = FALSE,min = 0.01){
dim= dim(genoData)
for(i in 1:dim[[1]]){
  for(j in 1:dim[[2]]){
    vec = genoData[i,j,]
    vec[vec<min] = 0.0
    vec  = vec * (.rescaleV/sum(vec))
    if(integer){ 
      vec = round(vec)
      maxind = which(vec==max(vec))
      vec[maxind] = vec[maxind] + (.rescaleV - sum(vec))
    }
        genoData[i,j,]=vec
  }
}
genoData
}

.getAvgWeights<- function(x1){
mean((x1[2,x1[3,]>0]))
}


.makeGenoWeights<-function(genoData,rescale = 1, ordinal=FALSE){
   dimg = dim(genoData)
   geno_header = dimnames(genoData)[[2]]
   indiv = dimnames(genoData)[[1]]  
  if(length(dimg)==2){
      genoweights = .getArray(list(geno_header,c("geno","weights","include"),indiv))
      genoweights[,1,] = t(genoData)
      genoweights[,2,] =  t( matrix(1,nrow = length(indiv),ncol = length(geno_header)))
      genoweights[,3,] = !apply(as.matrix(genoweights[,1,]),c(1,2),is.na)
         resc = rep(1,dim(genoData)[2])
  }else{

    genoData = .rescale(genoData,rescale,integer = (rescale>1 | ordinal),min = 0)
    valsg = as.numeric(dimnames(genoData)[[3]])
    len = length(indiv)
    genoweights = .getArray(list(geno_header,c("geno","weights","include"),rep(indiv,length(valsg))))
    for(i in 1:(length(valsg))){
	 genoweights[,1,((i-1)*len +1):(i*len)] = rep(valsg[i],length(indiv))
	 genoweights[,2,((i-1)*len +1):(i*len)] = t(genoData[,,i])
         genoweights[,3,((i-1)*len +1):(i*len)] = t(genoData[,,i]>0)
    }
    # resc = apply(genoweights,1,.getAvgWeights)
     ###else{
          resc = rep(rescale,dim(genoData)[2])
###      }
     #apply(apply(genoData,c(1,2),.entropy),2,sum)/dim(genoData)[1]

     
  }

  attr(genoweights,"dimg")<-dim(genoData)
  attr(genoweights,"rescale")<-resc
#  attr(genoweights,"levels")<-apply(genoweights,1,.getLevels)
#  print(resc)
   genoweights
}

.entropy<-function(v){
v = v/sum(v)
e = 0
for(k in 1:length(v)){
 if(v[k]>0)  e = e+-v[k]*log(v[k])

}
1+e
}

#assumes as baseLevel of 0, so should subtract 2 for CN
#note order of genoData is preserved
.mPhen <-
function(genoweights, pheno_, opts, subinds = (1:dim(genoweights)[3])){
  multiGen = FALSE
  rescale = attr(genoweights,"rescale")[subinds]
  JointModel=opts$mPhen.JointModel
  inverseRegress = opts$mPhen.inverseRegress
  
  if(JointModel){
    if(!inverseRegress) multiGen = TRUE
    else if(length(pheno_$phenN)==1) 	JointModel=FALSE
  }
  if(multiGen) JointModel=FALSE

  geno_header = dimnames(genoweights)[[1]]
  exactTest = !is.null(opts$mPhen.exactMethod)
#  exactMethod = "wald"
  multiPhen = JointModel
  singleOutput=TRUE
  zipOutput=FALSE
  step = 1
  if(length(pheno_$index)==1) multiPhen = FALSE
  if(length(pheno_$index_strat)>0 | length(pheno_$index_cov) > 0) exactTest = FALSE
  if(exactTest){
    multiPhen = FALSE
    multiGen = FALSE
  }
  funct = .getPvLrr
  pvFunct = .getPvLrrMult
  pvFunctMult = .getPvLrrAllGen


  if(inverseRegress) funct = .getPvRev
  if(multiPhen) {
   if(opts$mPhen.variable.selection) pvFunct = .getPvRevPleioVS
   else  pvFunct = .getPvRevPleio
  }
  if(opts$mPhen.scoreTest) {
      if(multiPhen) pvFunct = .ordTest1 
      else funct = .ordTestUnivariate
  }
  if(multiGen){
     if(opts$mPhen.variable.selection) pvFunctMult = .getPvLrrMultiGenVS
     else  pvFunctMult = .getPvLrrMultiGenNoVS
  }
  if(exactTest) funct = .getPvTable
  if(multiGen & dim(genoweights)[1]<=1) stop('need at least two genotypes for MultiGen model')
  geno2 = genoweights[,,subinds,drop=F]
  maxg = apply(geno2[,1,,drop=F],1,max,na.rm=TRUE)
  ming = apply(geno2[,1,,drop=F],1,min,na.rm=TRUE)
  nonint = maxg%%1 > 0
  nonint1 = ming%%1 > 0
  maxg[nonint] = floor(maxg[nonint])+1
  ming[nonint1] = floor(ming[nonint1])
  mini = min(max(maxg,na.rm=T),min(ming,na.rm=T)+1,na.rm=T)
  maxi = (max(maxg,na.rm=T))
  if(mini==-Inf) mini=0
  if(maxi==-Inf) maxi=3
  spl = seq(mini,maxi,step) - 0.5

  if(length(geno_header)==1) multiGen = FALSE
  fams <- pheno_$families
  if(inverseRegress & !is.null(opts$mPhen.link.geno)) fams<-rep(opts$mPhen.link.geno, length(fams))

  stratificNames = dimnames(pheno_$stratMatrix1)[[2]]
  if(multiGen) geno_header=c(geno_header,"comb") 
  if(length(stratificNames)==1) stratificNames = c("all")
  else stratificNames = c(stratificNames,"comb")
  phenN = pheno_$phenN
  if(multiPhen) phenN = c(phenN,"JointModel")

  res = .getArray(list(stratificNames,  geno_header,phenN,c("beta","pvalue","Nobs")))
  
 
  maf=NULL; 
   maf = rep(NA,dim(geno2)[1])
  for(k in 1:length(maf)){
     maf[k] = .getMaf(geno2[k,,,drop=F],0,2)
  }
hwe=NULL;
  if(opts$mPhen.calcHWE){
  hwe = rep(NA,dim(geno2)[1])
for(k in 1:length(hwe)){
     hwe[k] = .getHWE(geno2[k,,,drop=F],rescale[k])
   }
}
  strat = pheno_$stratMatrix1[subinds,,drop=F]
  for(j in 1:(dim(strat)[2])){
    inclu = strat[,j]
    res[j,,,] =  pvFunctMult(geno2,pheno_$phendata[,,subinds,drop=F],fams,inclu,pheno_$pheno_cov[subinds,,drop=F],rescale,pvFunct,funct,JointModel)
    if(getOption("mPhen.log10p",FALSE)){
    ## convert from loge to log10
      res[,,,2] = res[,,,2]/log(10)
     }
 }
  if(length(pheno_$index_strat)>0){ 
    res[dim(strat)[2]+1,,,] =apply(apply(res[1:(dim(strat)[2]-1),,,,drop=FALSE],c(2,3),.metaresInvVarianceFixed),c(3,1),t)
  }
 
  #pheno = pheno_$phendata
  #countsCase = array(NA,dim=dimcounts)
  #countsControl = array(NA,dim=dimcounts)
  #hwe_control = array(NA,dim=dimcounts[1:3])
  #hwe_case = array(NA,dim=dimcounts[1:3])
  #maf_control = array(NA,dim=dimcounts[1:3])
  #maf_case = array(NA,dim=dimcounts[1:3])
  #for(j in 1:(dimcounts[1])){
  #  stra = strat[,j]
  #  for(k in 1:length(families)){
  #     # print(paste(j,k))
  #    inclu = pheno[k,2,] & stra
  #    incluCont = inclu &  pheno[k,3,]
  #    incluCase = inclu &  pheno[k,3,]
  #    if(length(which(incluCont))>0)  countsControl[j,k,,] = t(apply(geno2[,1,incluCont,drop=FALSE],1,.getHist,spl))
  #    if(length(which(incluCont))>0)  countsCase[j,k,,] = t(apply(geno2[,1,incluCase,drop=FALSE],1,.getHist,spl))
  #    hwe_control[j,k,] =  apply(countsControl[j,k,,,drop=FALSE],3,.hwepall)
  #    hwe_case[j,k,] =  apply(countsCase[j,k,,,drop=FALSE],3,.hwepall)
  #    maf_control[j,k,] =  apply(countsControl[j,k,,,drop=FALSE],3,.mafCount)
  #    maf_case[j,k,] =  apply(countsCase[j,k,,,drop=FALSE],3,.mafCount)
  #  }
  #}
    list(Results=res,maf=maf,hwe=hwe)
}


##splits a string into component values
.splitImputed<-function(v1,tot=1000){
v = as.numeric(v1)   ####2014-1-1
ind = which(is.na(v))
v[ind] = tot -sum(v,na.rm=T)
return(v)
}



.getAvg<-function(line,vec) (vec%*%line)/sum(line)
.getAvgAll<-function(gen) apply(gen,c(1,2),.getAvg,as.numeric(dimnames(gen)[[3]]))


.nonNumeric <-
function(vec) sum(is.na(as.numeric(vec)))==length(vec)
.nonNumeric1 <-
function(mat) apply(mat,2,.nonNumeric)
.nullDim <-
function(matr){x = length(dim(matr)); x}
.numLev <-
function(vec) length(levels(as.factor(vec)))
.ordTest <-
function(gvar,phend,families, inclu,covar,totweight,functi){
  include =  apply(as.matrix(phend[,2,]),2,min) & inclu & as.numeric(gvar[2,])>0 &  gvar[3,]
  y = as.numeric(gvar[1,include])
  x = t(as.matrix(phend[,1,include] - phend[,4,include]))
  cats = max(y) + 1
  m = apply(x,2,mean,na.rm=TRUE)
  x = t(t(x) - m)
  d = ncol(x)
  gamma = as.vector(table(y)/length(y))
  phi = log(gamma)
  Unum1 = .getUNum1Matrix(d,y,x,phi,num1,cats)
  U = Unum1[[1]]
  num1 = Unum1[[2]]
  II = .getIMatrix(d,phi,cats,num1,x)
  zstats  = U/sqrt(diag(-II))
  pvs = 2*(1-pnorm(abs(U/sqrt(diag(-II)))))
  overallp = 	1 - pchisq(t(U) %*% solve(-II) %*% U, d)
  results = cbind(c(zstats, NA),c(pvs, overallp))
  results
}
.ordTest1 <-
function(gvar,phend,families, inclu,covar,totweight,functi){
 if(dim(covar)[2]>0) stop('ord test not coded up for covariates, use residuals')
  include =  apply(as.matrix(phend[,2,]),2,min) & inclu & as.numeric(gvar[2,])>0 &  gvar[3,]
  y = as.numeric(gvar[1,include])
  x = t(as.matrix(phend[,1,include] - phend[,4,include]))
  resu  = .ordTest2(y,x)
  resu = cbind(resu, rep(dim(covar)[1],dim(resu)[1]))
  resu
}




.ordTest2 <-
function(y,x,log.p=FALSE){
  cats = max(y) + 1
  m = apply(x,2,mean,na.rm=TRUE)
  x = t(t(x) - m)
  d = ncol(x)
  gamma = as.vector(table(y)/length(y))
  phi = log(gamma)
  U = matrix(ncol=1,nrow=d)
  num1 = list(length=d)
  for( j in 1:d ){
    U[j] = sum(y*x[,j])
    num1[[j]]=0
    for( i in 1:cats ){
      num1[[j]] = num1[[j]] + x[,j] * (i-1) * exp(phi[i])
    }
    U[j,1] = U[j,1] - sum(num1[[j]])
  }

  II = matrix(ncol=d,nrow=d)
  num2 = matrix(ncol=d,nrow=d)
  for( j in 1:d ){
    for( k in j:d ){
      num = 0
      for( i in 1:cats ){
        num = num + (i-1)*(i-1)*x[,j]*x[,k]*exp(phi[i])
      }
      II[j,k] = sum(num1[[j]]*num1[[k]])
      II[j,k] = II[j,k] - sum(num)
      II[k,j] = II[j,k]
    }
  }
  zstats  = U/sqrt(diag(-II))
  pvs = 2*(pnorm(abs(U/sqrt(diag(-II))),lower.tail=F,log.p=log.p))
  chisqv = t(U) %*% solve(-II) %*% U
  overallp = 	pchisq(chisqv, d,lower.tail=F,log.p=log.p)
  results = cbind(c(zstats, NA),c(pvs, overallp))
  results
}

##modified cholesky to keep column order
.chol<-function(cor,pivot=T){
ch = chol(cor,pivot=pivot)
if(pivot){
	 ch = (ch[,match(1:dim(ch)[[2]],attr(ch,"pivot"))])
}
if(!is.matrix(ch)) ch = matrix(ch)
dimnames(ch) = dimnames(cor)
ch
}

.genRandCovar<-function(noPhenos, orthogAll = 0.5, betav = 5, orthog_i = rbeta(noPhenos,orthogAll*betav,(1-orthogAll)*betav), sd = rgamma(noPhenos,shape=10,rate = 10), resample = FALSE,cor1 = NULL,dirichletScale = 5){
if(min(orthog_i)<0) stop(paste("orthog must be <1",min(orthog_i)))
if(noPhenos==1) return(as.matrix(1))

M = array(0,dim = c(noPhenos,noPhenos))

M[1,1] = 1.0
for(i in 2:noPhenos){

     nonorthog = 1-orthog_i[i]
     if(is.null(cor1)){
         v1 = .sampDirichlet(rep(1,i-1),weight=dirichletScale)*nonorthog
	 v = c(v1,1-nonorthog,rep(0,noPhenos - i))
     }else{
	     cor1[i] = min(nonorthog,cor1[i])
             v1 = .sampDirichlet(rep(1,i-2),weight=dirichletScale)*(nonorthog-cor1[i])
	     v = c(cor1[i], v1,1- nonorthog,rep(0,noPhenos - i))
     }
     M[,i] = sqrt(v) * sample(c(-1,1),length(v),replace=T)
}
res = (t(M) %*% M )* outer(sd,sd)
pnames = paste("pheno",(1:noPhenos),sep="")
if(resample){ 
 ord = sample(1:noPhenos)
 res = res[ord,ord]
}
dimnames(res) = list(pnames,pnames)
res
}

.reopenConnection<-function(outC){
  files = outC$files
  for(k in 1:length(files)){
      if(!is.null(files[k])) {
         nme = attr(files[k],"name")
	 files[[k]]<- file(nme, 'a')
  	 attr(files[[k]],"name")<-nme
 	 attr(files[[k]],"first")<-FALSE
 	 #attr(files[[k]],"open")<-TRUE
   }
  }
  outC$files = files
  outC$open=TRUE
  outC
}

##Prepares output files and directories.
##limit is obtained from mPhen.readLimitFile
.openOutputConnection<-function(outFiles="results/",genName,write, plot, limits=NULL, regopts = NULL
  ){
 if(!is.null(regopts)){
    nme = paste(outFiles,"/.opts_",genName,sep='')
    nmes = names(regopts) 
    charc = which(lapply(regopts,class)=="character")
    for(k in charc){
     regopts[[k]] = paste("\'",regopts[[k]],"\'",sep="")
    }
    x =  apply(cbind(paste("opts$mPhen.",nmes,sep=""),regopts),1,paste,collapse="= ")
    for(k in x)     write(k,file = nme)

 }

   nme = paste(paste(outFiles,"/",genName,sep=''),names(write),sep=".")
   files = list()
   for(k in 1:length(write)){
      if(write[[k]]) {
	 files[[k]]<- file(nme[k], 'w')
  	 attr(files[[k]],"name")<-nme[k]
 	 attr(files[[k]],"first")<-TRUE
 	 attr(files[[k]],"open")<-TRUE
      }else {
	files[[k]] = NULL
      }
  }
  
 names(files) = names(write)[1:length(files)]
  nmesplot = names(plot)

  plotFiles = paste(paste(paste(outFiles,"/",genName,sep=''),nmesplot,sep=""),".pdf",sep="")

   
   if(!is.null(limits)){
 	outFileLimit=paste(outFiles,"/limit_",genName,".txt",sep='')
	  write.table(limits,file=outFileLimit,col.names=FALSE,row.names=FALSE,sep="\t",append=F,quote=F)
   }
   list(files = files, plotFiles=plotFiles,pvs=list(), beta = list(), sampleQC=NULL,open=TRUE)
}



#converts multi-dimensional output array into a printable table
.flatten<-function(results, extra=dimnames(results$Res)[[2]],extra_head="pos"){
 res = results$Results
 hwe = results$hwe
 if(is.null(hwe)) hwe = rep(NA,length(results$maf))
 maf = results$maf
 dn = dimnames(res)
 dims=dim(res)
 norows=dims[1]*dims[2]*dims[3]
 matr = .getArray(list(1:norows,c(extra_head,"strata","snpid","phen","maf","hwe","beta","pv","nobs")))
 n=1
 for(i in 1:dims[1]){
     for(j in 1:dims[2]){
      for(k in 1:dims[3]){ 
      line =         c(extra[j],dn[[1]][i], dn[[2]][j],dn[[3]][k],sprintf("%2.3g",c(maf[j],hwe[j],res[i,j,k,])))
          matr[n,] = line
	  n =n+1
      }
    }
 }
matr
}

## method for sub-sampling positions for bootstrap which 
## makes sure it preserves same imputed samples
.subSample<-function(samples, dim3){
  # samples = 1:dim1
  # if(resample) {
  #    samples = sample(samples, dim1,replace=TRUE)
  # }
  dim1 = length(samples)
   samples1 = rep(NA, dim3*length(samples))
   for(k in 1:dim3){
     samples1[(k-1)*(dim1) + (1:(dim1))] = (k-1)*dim1 + samples
   }
     return(samples1)
}





#converts abs(beta*X) to printable form
.flattenBetaX<-function(sumBetaX){
   dims=dim(sumBetaX)
   dn = dimnames(sumBetaX)
   nrows=dims[1]*dims[2]
   matr =.getArray(list(1:nrows,c("PATIENT",dn[[3]])))
   n=1
    for(i in 1:dims[1]){
    for(j in 1:dims[2]){

	matr[n,] = c(paste(dn[[1]][i],dn[[2]][j],sep="_"),round(10*sumBetaX[i,j,])/10)
	n=n+1
    }
   }
   matr
}

##formats table for printing
##table is a table of  numerical values. 
## format is a string in sprintf format (e.g. "%5.3f")
## replaces null values with '-' for ease of reading
.formatTable<-function(table,format,nullv=0,nullrepl = "   -    ",  firstColHeader="SNP"){
  if(!is.matrix(table)) table = t(as.matrix(table))
  nmes=dimnames(table)[[1]]
  dim = dim(table)
  res = array(dim = c(dim[1],dim[2]+1))
  res[,1] = nmes
  for(i in 1:dim[1]){
   for(j in 1:dim[2]){
     if(is.na(table[i,j])){
		res[i,j+1] = sprintf(format,table[i,j])
     }else if(table[i,j]==nullv){
	     res[i,j+1] =  nullrepl

     } else {
	res[i,j+1] = sprintf(format,table[i,j])
     }
   }
  }
  dimnames(res) = list(NULL,c(firstColHeader,dimnames(table)[[2]]))
  res
}


##function for printing out results of simulation with a defined correlation
##matrix cor, effect direction v and betaAll
.writeTable<-function(resultsJ,file, rsid =dimnames(resultsJ$Res)[[2]],  format = "%3.2e",betaAll=NULL, v=NULL,cor=NULL){
  append=!attr(file,"first")
  index = 1
  resu = resultsJ$Res
  strats = dimnames(resultsJ$Res)[[1]]
  len =  dim(resu)[2]
  nmes = dimnames(resu)[[2]]
  numSnps = dim(resu)[[2]]
  numPheno = dim(resu)[[3]]
  label = rep("pval",len)
  label1 = rep("beta",len)
  for(k in 1:length(strats)){
   strata = rep(strats[k],len)
   beta = resu[k,,,1]
   pvalue = resu[k,,,2]
   if(is.null(dim(pvalue))){
    pvalue = array(data = pvalue,dim = c(dim(resu)[2],ncol = dim(resu)[3]),dimnames = dimnames(resu)[2:3])
   }
   if(is.null(dim(beta))){
     beta = array(data = beta,dim = c(dim(resu)[2],ncol = dim(resu)[3]),dimnames = dimnames(resu)[2:3])
   }

   if(!is.null(betaAll)){
     pvalue = cbind(betaAll,pvalue)
     beta = cbind(betaAll,beta)
  }
  write.table(file=file,cbind(strata,label,rsid,.formatTable(pvalue,format,nullv=1.0,nullrepl="   1    ")),row.names=F,quote=F,append=append, col.names=!append)
label = rep("beta",len)
 write.table(file=file,cbind(strata,label,rsid,.formatTable(beta,format,nullv=1.0,nullrepl="   1    ")),row.names=F,quote=F,append=append, col.names=FALSE)

  }
}

##Writes CNV genotype information in a bed format
## input is object returned by mPhen.readGenoConnection
## output is object returned by mPhen.openOutputConnection
##pheno_ is object returend by mPhen.preparePheno
##bed output is split into different groups based on the pheno values in pheno_, according to pre-specified quantiles
.writeBed<-function(pheno_,input,output,index, .quantiles = c(1/3,2/3,1),.base=getOption("mPhen.baseValue",2),
nme=paste(getOption("mPhen.resultsName","result"),index,sep="."),
phenN = pheno_$phenN){

 fil1=file(paste(output$outFileBed,nme,"bed",sep="."),'w')
  g0 = input$genoData
 cols = .makeColsForBed()
 pos = t(as.matrix(data.frame(strsplit(dimnames(g0)[[2]],"_"))))
 if(length(levels(as.factor(pos[,1])))>1) stop('only one chrom at time')
 chr = pos[1,1]
 starts = as.numeric(pos[,2])
 minpos = input$firstPos-1000
 lastP = max(starts[length(starts)]+1000,input$lastPos)
 maxpos = lastP +1000
 linestart=paste("browser position ",paste("chr",chr,sep=""),":",minpos,"-",maxpos,sep="")
    write(linestart,file=fil1,append = FALSE)
 nmes = dimnames(g0)[[1]]
 for(k1 in 1:length(phenN)){
 # phenN = phenN[k1]
    k = which(pheno_$phenN==phenN[k1])
    typenme = phenN[k]
 phend = pheno_$phendata[k,1,]
 quantiles = unique(quantile(phend,.quantiles))
 counts = rep(0,length(quantiles))
 prevj = -Inf
 for(j in 1:length(counts)){
  indsj = which(phend<=quantiles[j] & phend>prevj)
  if(length(indsj)>0){
  
    counts[j] = length(indsj)
    prevj = quantiles[j]
    if(j==1) dirnme = "Controls" 
    else if(j==length(quantiles)) dirnme = "Cases" 
    else dirnme = paste("Q",quantiles[j],sep="")
    line1 = "track name=\"Controls\" description=\"Controls\" itemRgb=\"On\"" 
     trackname = paste(nme,typenme,dirnme,sep="_")
    line1=gsub("Controls",trackname, line1)
 
    write(line1,file=fil1,append=T)

    for(k in indsj){

     arr1 = .getCNProfile(g0[k,],starts,lastP,.base,cols,nmes[k],paste("chr",chr,sep=""))
      write.table(arr1,row.names=F,col.names=F,quote=F,sep="\t",file=fil1,append=T)
    }
    }
 }
}
if(FALSE){
  snps = pos
  bedSnp = cbind(rep(paste("chr",chr,sep=""),dim(snps)[1]),snps[,2],snps[,2])
   typenme = paste("expt","probes",sep="_")
   dirnme = ""
   line1 = "track name=\"Controls\" description=\"Controls\" itemRgb=\"On\"" 
   line1=gsub("Controls",paste(typenme,dirnme,sep="_"), line1)
   write(linestart,file=fil1,append = T)
   write(line1,file=fil1,append=T)
   write.table(bedSnp,row.names=F,col.names=F,quote=F,sep="\t",file=fil1,append=T)
}
stratNames=1:length(output$pvs)
 for(kk in 1:length(output$pvs)){
  topl = output$pvs[[kk]]
  beta = output$beta[[kk]]
  stratN = stratNames[kk]
  .manhattanBed(topl,beta, fil1,paste("pvals",nme,stratN,sep="_"))
  
  }
close(fil1)
}


.summariseBootstrap<-function(res, thresh = 0.05){
resAll = array(0,dim = dim(res[[1]]$Res), dimnames = dimnames(res[[1]]$Res))
for(k in 1:(dim(resAll)[[2]])){
  for(j in 1:length(res)){
    summ = which(res[[j]]$Res[,k,,2]<thresh)
    resAll[,k,summ,2] =  resAll[,k,summ,2]+1
    if(j==1) resAll[,k,summ,1] =  res[[j]]$Res[,k,summ,1]
    else     resAll[,k,summ,1] =  resAll[,k,summ,1]+res[[j]]$Res[,k,summ,1]
  }
  for(m in 1:(dim(resAll)[[3]])){
    resAll[,k,m,1] = resAll[,k,m,1]/resAll[,k,m,2]
  }
   resAll[,k,,2] = resAll[,k,,2] / length(res)
}
resAll[,,,3] = res[[1]]$Res[,,,3]
resAll
}


.parseOpts<-function(){
args = commandArgs()
argsst=which(args=="--args")
if(length(argsst)>0){
  if( argsst <= length(args)){
  args = args[(argsst+1):length(args)]
  print(args)
  eqind = grep("=",args)
  if(length(eqind)>0){
    argsM = as.matrix(data.frame(strsplit(args[eqind],"=")))
    toassign = as.list(argsM[2,])
    names(toassign) = argsM[1,]
    options(toassign)
  }
 }
}
}


.getLimitMatrixFromList<-function(l1){
  limit = cbind(rep('pheno', .length(l1$phenotypes)), l1$phenotypes)
  if(.length(l1$covariates>0)) limit = rbind(limit, cbind(rep('covar', .length(l1$covariates)), l1$covariates))
  if(.length(l1$resids>0)) limit = rbind(limit, cbind(rep('resid', .length(l1$resids)), l1$resids))
  if(.length(l1$strats>0)) limit = rbind(limit, cbind(rep('strat', .length(l1$strats)), l1$strats))
  if(.length(l1$excls>0)) limit = rbind(limit, cbind(rep('excl', .length(l1$excls)), l1$excls))
  limit
}

.readLimitFile<-function(limitfile="./limit.txt", phenNames=NULL){
 todo_ = read.table(limitfile,as.is=T,fill=T,header=F,comment.char='')[,1:3,drop=F]
 hashinds = grep('^#',todo_[,1])
 if(length(hashinds)>0)   todo_ = todo_[-hashinds,,drop=F]
  ind2 = grep(':',todo_[,2])
  while(length(ind2)>0){
        k= ind2[1]
	ln = strsplit(todo_[k,2],':')[[1]]
	#pref = strsplit(ln[1],"[0-9]")[[1]][1]
        #indin = nchar(pref)
        #if(indin>0) ln = substr(ln,start=(indin+1),stop=200)
        if(nchar(ln[[1]])==0)  newt=phenNames[as.numeric(ln[2]):as.numeric(ln[3])]
        else  newt=paste(ln[[1]],as.numeric(ln[2]):as.numeric(ln[3]),sep="")
     
        todo_ = rbind(as.matrix(todo_[-k,]),as.matrix(cbind(todo_[k,1],newt,todo_[k,3])))
        ind2 = grep(':',todo_[,2])
  }
  ind2 = grep('\\*',todo_[,2])
  while(length(ind2)>0){
        k= ind2[1]
	ln = strsplit(todo_[k,2],'\\*')[[1]][1]
        newt = grep(ln,phenNames,value=T)
        todo_ = rbind(as.matrix(todo_[-k,]),as.matrix(cbind(todo_[k,1],newt,todo_[k,3])))
        ind2 = grep(':',todo_[,2])
  }

  todo_ = unique(todo_)##
  
  return(todo_)
#  todo_ = .fixDuplicatesInTodo(todo_)
#  print(todo_)
#outn= .getSuffix(todo_)
#  outn=substr(outn,1,50)
#list(todo = todo_,outn =outn)
}




## updates sample qc based on current set of results
.updateSampleQC<-function(genoData,results,sampleQC = NULL, pvthresh =0.01,maf_thresh = 0.05, hwe_thresh = 1e-4){
   res1 = results$Results
   sumBetaXNew =.getSumBetaX(genoData,res1,results$maf,results$hwe,pvthresh = pvthresh,maf_thresh = maf_thresh,hwe_thresh = hwe_thresh)
   if(is.null(sampleQC)) sampleQC = sumBetaXNew else sampleQC= sampleQC+sumBetaXNew
   sampleQC
}


#calculates sum(abs(beta*x)) for each sample over all snps to identify
#individuals driving the result
## not used if  using inverse regression
.getSumBetaX<-function(genD,res,maf,hwe,pvthresh=0.05,maf_thresh=0.05, hwe_thresh = 1e-4){
 sumBetaX=NULL;
 nosnp = length(dimnames(res)[[2]])
 if(is.null(maf)) maf = rep(0,nosnp)
 if(is.null(hwe)) hwe = rep(1,nosnp)
  geno_header = dimnames(genD)[[1]]
 if(length(dim(genD))>2) genoData=.getAvgAll(genD)
 else genoData = genD
  if(!is.null(maf) & ! is.null(hwe)){
  li = dimnames(res)[c(1,3)]
  li[[3]] = geno_header
  sumBetaX = .getArray(li,v=0)
   dm = dim(sumBetaX)
      for(k in 1:nosnp){
	    nonNAind = !is.na(genoData[,k])
       if( maf[k] < maf_thresh & hwe[k] > hwe_thresh){
       for(i in 1:dm[1]){
        for(j in 1:dm[2]){
           pv = res[i,k,j,2]
           beta = res[i,k,j,1]
	if(!is.na(beta)){
 	   if(pv<pvthresh){
                 sumBetaX[i,j,nonNAind] = sumBetaX[i,j,nonNAind]+abs(beta*(genoData[nonNAind,k]))
              	   }
          }
       }
       }

      }
  }
 }
  sumBetaX
}

  .sample1<-function(probv)  sample(0:1, 1, prob = probv)

##simply samples from a dirichlet with probability distribution p
## and weight w.
.sampDirichlet<-function(p,weight=1){
alpha = p*weight
scale = 1
vec = rep(0,length(alpha))
for(i in 1:length(vec)){
  vec[i] = rgamma(1,shape=alpha[i],scale=scale)
}
vec = vec/sum(vec)
vec
}



###HAS BEEN MODIFIED 29-11
##calculates correlation matrices and effect sizes
#beta0 is expressed as proportion of covar[1,1,]
#v is a rotational vector such that after rotation main effect is in direction of rotation
#returns cholesky decomposition of rotated and non-rotated phenotypes
## phenos is a simulated dataset from the covariance matrix, unless specified (in which case cov has no effect and can be NULL)
.cholesky<-function(cov, .nsamp = 1e5,      
    v =NULL,
    phen = if(is.null(v))NULL else  .sampJoint(.chol(cov2cor(cov),pivot=T),num=.nsamp,sd=sqrt(diag(cov))),
    means= if(is.null(phen)) rep(0,dim(cov)[[2]]) else apply(phen,2,mean))
 { 
   if(!is.null(v)){
       nonzero = which(v!=0)
       if(length(nonzero)==0) v = NULL
       else if(length(nonzero)==1 & nonzero[1]==1) v= NULL
   }
   if(is.null(v)){
      M = diag(rep(1,dim(cov)[2]))
      covT = cov
   }else{
     M = .gs(v) 
  #   Minv = t(M)
 ### note v %*% M is on x axis
     phenosT = phen %*% M  ###t(Minv%*%t(phen))
     covT = cov(phenosT,use="complete.obs")
   }
    corT = cov2cor(covT)
    dimnames(corT) = list(dimnames(cov)[[2]],dimnames(cov)[[2]])
    cholT = .chol(corT,pivot=T)
    list(cholT = cholT,  Minv = t(M),sdT = sqrt(diag(covT)),means=means,corT=corT)
}


.pow<-function(x,base=10,normalise=TRUE,multFactor=1000){
   v = (10^as.numeric(x))
   if(normalise) v = v/sum(v)
   if(multFactor>1){
	 v = round(v*multFactor)
	 maxind = which(v==max(v))
	 v[maxind] = v[maxind] + multFactor - sum(v)
   }
   return(v)
}

.unchanged<-function(v1) as.numeric(v1)
.cnvdistr<-function(v1){
 cn = rep(0,5)
 for(k in seq(1,length(v1),2)){
    p = as.numeric(v1[k+1])
    CN =nchar(sub(",","",v1[k]))+1 
    cn[CN] = round(1000*p)
 }
 cn
}
  
#  if(length(v)==3){
#     v = v[!(v==min(v))]
#   }
#   v
#}



.bafratio<-function(v) v[2]/(v[2]+v[1])
.bafratio<-function(v) v[2]/(v[2]+v[1])

.plotBaf<-function(geno,  usesum=FALSE, thresh = 50){
x = as.numeric(as.matrix(data.frame(strsplit(dimnames(geno)[[2]],'_')))[2,])
sumy = apply(geno,c(1,2),sum)
if(usesum) { plot(x,sumy)
}else{
  y = apply(geno,c(1,2),.bafratio)
  
 plot(x[sumy>thresh],y[sumy>thresh])
}
}

##FOR SIMULATION ######
##Aim is to project out all variations corresponding to vars variables (this is just indices)
## from the matrix Dall.
## removes all variation in P from Dall
## returns the vector v such that Dall - P%*%v is orthogonal to all columns
## in P 
.projOut<-function(Dall, vars = NULL, P = Dall[,vars,drop=F]){
  ncol = dim(Dall)[[2]]
  nrow = dim(Dall)[[1]]
  svd = svd(P)
  U = svd$u
  V = t(svd$v)
  if(dim(P)[2]==1) {
      D = as.matrix(svd$d)
  }else {
      D = diag(svd$d)
  }
  Vinv = solve(V)
  Dinv = solve(D)
  colU = dim(U)[[2]]
  proj = matrix(0,nrow = dim(U)[[2]],ncol = ncol)
  for(j in 1:ncol){ 	
	 for(k in 1:colU){
           proj[k,j] = U[,k] %*% Dall[,j]
	}	
  }  
  Vinv %*% Dinv %*% proj
}

## returns an orthogonal matrix M, such that 
## M %*% v  = c(v%*%v,0,0,0)
.gs<-function(v, thresh = 1e-8){
 ###res = svd(cbind(v,basis))$u  ### not in same direction
 v1 = v/sqrt(v%*%v)
 basis=cbind(v1,diag(rep(1,length(v))))
 W = .projOut(basis,c(1))
 basis1 = basis - basis[,1,drop=F] %*% W
 svd = svd(basis1[,-1,drop=F])
 res = cbind(v1,svd$u[,svd$d>thresh,drop=F])
 res
 ###t(res) ## this is same as inverse as orthonormal
}

#basic function for simulating effect of genotype of phenotype in a linear direction
.sampJoint<-function(ch,x=NULL,num=length(x),varexp=0,
                ncol =  dim(ch)[2],
		means=rep(0,ncol),
			sd = rep(1,dim(ch)[2])){

	
        res = matrix(rnorm(num*ncol,0,1),nrow =num, ncol = ncol)
        if(!is.null(x) & varexp>0){
             beta = sqrt(((varexp/(1-varexp)))/var(x))
             res[,1] =res[,1] + x*beta
        }
	res = res %*% ch     
              #t(t(ch)%*%t(z1))

        for(i in 1:(dim(res)[2])){
          res[,i] = .destandardise(res[,i],means[i],sd[i])
	}

        dimnames(res)[[2]] =  dimnames(ch)[[2]] 
res
}

##END FOR SIMULATION ### 


#### FOR QQ PLOTTING MULTIPLE VARIABLES ###


.qb<-function(np,conf,ntot,sigThresh=1.0){
  logp = getOption("mPhen.log10p",FALSE)
  Obs = rep(0,np)
  if(logp) sigThresh = 10^sigThresh;
  for(s in 1:np)
              {
	        prop = (s/(np+1))*sigThresh
                s1 = ceiling(prop*(ntot+1))
                Obs[s]=qbeta(conf, s1, ntot-s1+1)
              }
 if(logp) return(log10(Obs))
else return( Obs)
}


.getPv<-function(pvs1,effPhe, opts, beta=NULL){
logP = getOption("mPhen.log10p",FALSE)
.sigThresh =opts$mPhen.sigThresh
.limitP = opts$mPhen.limitP
if(logP){
	 .limitP = log10(.limitP)
	.sigThresh = log10(.sigThresh)
}
pool=opts$mPhen.pool
pvs = pvs1
dimp = dim(pvs)
nosnp = dimp[[1]]
nophe = dimp[[2]]
np = nosnp*nophe
ntot = effPhe*nosnp
pv = matrix(0,nrow = np,ncol=7)
nmes = rep(dimnames(pvs)[[1]],nophe)
dimnames(pv)= list(nmes,c("pv","index",0.5,0.95,0.05,0.5,"beta"))
for(i in 1:nophe){
   indexes = ((i-1)*nosnp+1):(i*nosnp)
   pv[indexes ,2] = rep(i,nosnp)
   pv[indexes ,1] = pvs[,i]
   if(!is.null(beta)) pv[indexes, 7] = beta[,i]
}
pv = pv[!is.na(pv[,1]),,drop=F]
pv = pv[pv[,1]<=.sigThresh,,drop=F]
pv[pv[,1]<.limitP,1] = .limitP
pv = pv[order(pv[,1]),,drop=F]
np = dim(pv)[[1]]
if(pool){
 pv[,3] = .qb(np,0.5, ntot,.sigThresh)  #(1:np)/(np+1)
 pv[,4] = .qb(np,0.95, ntot,.sigThresh)
 pv[,5] = .qb(np,0.05, ntot,.sigThresh)
}else{
 for(i in 1:nophe){
        indexes = which(pv[,2]==i)
	nosnp1 = length(indexes)
    	 pv[indexes,3] = .qb(nosnp1,0.5, nosnp1,.sigThresh)  #(1:np)/(np+1)
         pv[indexes,4] = .qb(nosnp1,0.95, nosnp1,.sigThresh)
         pv[indexes,5] = .qb(nosnp1,0.05, nosnp1,.sigThresh)
    }
}
pv[,6] = pv[,3]

pv
}
#.getPvAll<-function(pvs1,effPhe,.sigThresh=1.0){
#  pv = .getPv(pvs1[,1,drop=F],.sigThresh)
#  for(k in 2:(dim(pvs1)[[2]])){
#    pv1 = .getPv(pvs1[,k,drop=F],.sigThresh)
#    pv1[,2] = rep(k,dim(pv1)[[1]])
#    pv = rbind(pv,pv1)
#  }
#  pv=pv[order(pv[,1]),]
#  pv
#}


.manhattanBed<-function(pvs1,logfun,beta, title,fil1,opts){
effPhe = 1.0
 order = order(apply(pvs1,2,min,na.rm=T))
pv <-.getPv(pvs1,effPhe, opts,beta)
if(dim(pv)[[1]]>0){
nophe = dim(pvs1)[[2]]
pos_info = as.matrix(((data.frame(strsplit(dimnames(pv)[[1]],"_")))))
pos=as.numeric(pos_info[2,])
end=as.numeric(pos_info[3,])
chroms=pos_info[1,]
chromsl = levels(as.factor(chroms))
chr = paste("chr",chromsl[1],sep="")
if(length(chromsl)>1) stop('only works for one chrom')
 for(i in 1:nophe){
       subind = pv[,2]==order[i]
       data = pv[subind,1]
       posj = pos[subind]
       endj = end[subind]
       betaj = pv[subind,7]
    line1 = "track type=bedGraph name=\"Controls\" description=\"Controls\" itemRgb=\"On\"  graphType=bar" 
    line1=gsub("Controls",paste(title,dimnames(pvs1)[[2]][order[i]],sep="_"), line1)
 
    write(line1,file=fil1,append=T)
    for(j in 1:length(data)){
     
      logpv = -logfun(data[j])
      if(betaj[j]>0) logpv = -logpv
      posjj = posj[j]
      write(paste(chr,posjj,endj[j],sprintf("%5.3g",logpv),sep=" "),file=fil1,append=T)
    }
 
 }

}
}



### Assumes snp names have format chr_pos
.manhattan<-function(pvs1,logfun, beta1, title=expression(paste("-log10","Manhattan")),opts){
 order = order(apply(pvs1,2,min,na.rm=T))
 effPhe = 1
step=opts$mPhen.noPhenosPerPlot
chromsToDo = opts$mPhen.chromsToDo
epsilon = opts$mPhen.epsilon
pv =.getPv(pvs1,effPhe,opts, beta1)
if(dim(pv)[[1]]>0){
nophe = dim(pvs1)[[2]]
ymax=max(8,logfun(min(pv[,1])))
norep = ceiling(nophe/step)
pos_info = as.matrix(((data.frame(strsplit(dimnames(pv)[[1]],"_")))))
#print(pos_info)
#print(dimnames(pv))
pos=as.numeric(pos_info[2,])
chroms=pos_info[1,]
toincl = !is.na(pos)
if(length(which(toincl))==0){
  pos = 1:length(pos)
  toincl = rep(TRUE,length(pos))
}
if(!is.null(chromsToDo)){
  toincl = chroms%in%chromsToDo & toincl
}

  pos_info = pos_info[,toincl]
  pos= pos[toincl]
  chroms = chroms[toincl]
  pv = pv[toincl,]

chromsl = levels(as.factor(chroms))

if(length(chromsl)>1){
for(k in chromsl){
    indk= which(chroms==k)
    posk = pos[indk]
    mink = 0  #
    maxk=max(posk)
    pos[indk] = ((posk - mink)/(maxk-mink))*(1-2*epsilon)+epsilon

}

   if(length(grep("_",chroms))>0){
     chroms = as.matrix(((data.frame(strsplit(chroms,"_")))))[1,]
   }
   chroms[chroms=="X"]=23
   chroms[chroms=="Y"]=24
   chroms[chroms=="XY"]=25
   chroms = as.numeric(chroms)
   pos = pos + chroms   
}

xstart = min(pos)
xend = max(pos)
xstart = xstart  - 0.1*(xend-xstart)

for(jk in 1:length(chromsl)){
 for(kk in 1:norep){
 start = ((kk-1)*step) 
 index=(start+1):min(kk*step,nophe)
 for(i in index){
       subind = pv[,2]==order[i]
       data = pv[subind,1]
       chrs = chroms[subind]
       subind1 = which(chrs==chromsl[jk])
         col = jk
         if(!opts$mPhen.colourByChroms) col = i-start
       if(jk==1 & i==(start+1)){
	 plot(pos[subind][subind1],logfun(data[subind1]),   pch=i-start, cex=opts$mPhen.cex,col=col,type="p",
		xlim = c(xstart,xend),ylim=c(0,ymax), ylab=expression(paste("-log"[10]," (Expected p-value)")), 
	   xlab=paste("Position on",pos_info[1,1]), main=title)
       }
       else  lines(pos[subind][subind1],logfun(data[subind1]),   pch=i-start, cex=opts$mPhen.cex,col=col,type="p")
 
 }
txtsize = min(0.5, 25/length(index))
legend("topleft", legend =dimnames(pvs1)[[2]][order[index]],cex=txtsize,col=(1:step),pch=1:step)
}
}
}
}



.qqplot<-function(pvs1,effPhe, logfun, beta1, title=expression(paste("-log"[10],"QQ")),opts){
 order = order(apply(pvs1,2,min,na.rm=T))
pv =.getPv(pvs1,effPhe, beta=NULL, opts)
step = opts$mPhen.noPhenosPerPlot
nophe = dim(pvs1)[[2]]


 ymax=logfun(min(pv[,1]))
 xmax=logfun(min(pv[pv[,3]>0,3]))

xmax = min(50,xmax)
norep = ceiling(nophe/step)
if(length(pv[,1])==0){
  xmax = 1
  ymax = 1
}
xmax = max(xmax,ymax)
ymax = xmax
for(kk in 1:norep){
 ## first plot all invisibly
 plot(logfun(pv[,c(3,1)]),xlim = c(0,min(opts$mPhen.minlogpv,xmax)), ylim=c(0,min(opts$mPhen.minlogpv,ymax)), xlab=expression(paste("-log"," (Expected p-value)")), ylab=expression(paste("-log"," (Observed p-value)")), main=title, pch=20, cex=0.5,col=0)

 lines(logfun(pv[,6]), logfun(pv[,6]), col="grey", lw=2)
 lines(logfun(pv[,4]), logfun(pv[,3]), col="grey", lw=2)
 lines(logfun(pv[,5]), logfun(pv[,3]), col="grey", lw=2)
 start = ((kk-1)*step) 
 index=(start+1):min(kk*step,nophe)
 for(i in index){
       data = pv[pv[,2]==order[i],c(3,1),drop=F]
       lines(logfun(data),   pch=i-start, col=i-start,type="p",cex=opts$mPhen.cex)
 }

txtsize = min(0.5, 25/length(index))
legend("bottomright", legend =dimnames(pvs1)[[2]][order[index]],cex=txtsize,col=(1:step),pch=1:step)
}
}


## Makes a QQ plot from a table of pvalues.  Each column is pvalue from a different phenotype, and each row is a different marker
##pvs can also be a list of pvalue tables (e.g. by strata), in which case multiple qqplots are calculated
# dimnames of table is in format ("position\tchr")
.qq<-function(pvs,effPhe,logfun,
              opts =list( 
                         sigThresh=getOption("mPhen.sigThresh",1.0),
                         step=getOption("mPhen.noPhenosPerPlot",20),
                         limitP = 1e-200,
                         title = "",
                         minlogpv=10,
                         cex=0.5,pool=TRUE))
{
  if(!is.list(pvs)) pvs = list(pvs)
  stratNames=names(pvs)
  if(is.null(stratNames)) stratNames = 1:length(pvs)
  for(kk in 1:length(pvs)){
    topl = pvs[[kk]]
    .qqplot(topl,effPhe,logfun,beta1 = NULL,  paste(opts$mPhen.title,stratNames[kk],"logQQ Pool = ",opts$mPhen.pool, "effPhe ",effPhe),opts = opts)
  }
}

## Makes a manhattan plot from a table of pvalues.  Each column is pvalue from a different phenotype, and each row is a different marker
##pvs can also be a list of pvalue tables (e.g. by strata), in which case multiple qqplots are calculated
##dimnames(pv)[[1]]) must be of  format ("position\tchr")
## only p-values less than .sigThresh are plotted
## only 'step' phenotypes are plotted on each manhattan plot
#.sigThresh=.sigThresh,
.manh<-function(pvs,effPhe,logfun, opts = list(
                                               sigThresh=getOption("mPhen.sigThresh",1.0),
                                               step=getOption("mPhen.noPhenosPerPlot",20),
                                               limitP = 1e-200,
                                               title="", 
                                               minlogpv=10,
                                               cex=0.25,
                                               chromsToDo = NULL,colourByChroms=FALSE))

{
  if(!is.list(pvs)) pvs = list(pvs)
  stratNames=names(pvs)
  if(is.null(stratNames)) stratNames = 1:length(pvs)
  for(kk in 1:length(pvs)){
    topl = pvs[[kk]]
    stratN = stratNames[kk]
    .manhattan(topl,logfun, beta1=NULL, paste(opts$mPhen.title,stratN,"Manhattan"),opts = opts)
  }
}

##pv heatmap fuction
.heatm<-function(pvs,effPhe, logfun, opts =list(title="",indexMatch = NULL, Colv=TRUE),symbreaks = FALSE, file ="", xlab="phenotypes",ylab="snp"){
  if(!is.list(pvs)) pvs = list(pvs)
  title=opts$mPhen.title
  if(symbreaks) col = c(rev(colorRampPalette(brewer.pal(9,"Reds"),bias = 0.5)(250)), colorRampPalette(brewer.pal(9,"Blues"), bias = 0.5)(250))
  else col=colorRampPalette(brewer.pal(9,"Blues"))(250)
  Colv=opts$mPhen.Colv
  indexMatch = opts$mPhen.indexMatch
  dendro = Colv
  for(k in 1:length(pvs)){
    main = paste(title,names(pvs)[[k]], if(symbreaks) "beta" else "pv")
    pvm = pvs[[k]]
    if(is.null(pvm)){
      return(NULL)
    }
    if(!is.matrix(pvm)) pvm = as.matrix(pvm)
    heatMatrix66 = logfun(t(pvm))
    d = dim(heatMatrix66)
    if(d[1]==1) heatMatrix66 = rbind(heatMatrix66,rep(0,d[2]))
    d = dim(heatMatrix66)
    if(d[2] ==1) heatMatrix66 = cbind(heatMatrix66,rep(0,d[1]))  #1e-10*runif(d[1]))
    d = dim(heatMatrix66)
    #if(k==1){
     #if(!is.null(Colv) & (is.numeric(Colv) | is.logical(Colv))){
      # if(!is.na(Colv)){
      #  naind = length(which(!is.na(as.vector(heatMatrix66))))
      #   if(Colv & naind==0) dendro<-NA
      # }
      #}

    #}
    dn = dimnames(heatMatrix66)
    cols66<-rep("white",ncol(heatMatrix66))
    if(!is.null(indexMatch)){
      indexChange<-grep(indexMatch,colnames(heatMatrix66))
      cols66[indexChange]<-"red"
    }
    cexCol = 0.5*min(1,0.2 + 1/log10(d[1]))
    cexRow =  0.5*min(1,0.2 + 1/log10(d[2]))
  #  torem = which(apply(heatMatrix66,2,.cntNa)==dim(heatMatrix66)[1])
    #print(heatMatrix66)
    #if(abs(max(heatMatrix66,na.rm=T) - min(heatMatrix66,na.rm=T))<1e-5) return(NULL)
     #print(heatMatrix66)
     #print(k)
     #print(dendro)
#	print(paste(file,main))
      p3 <- tryCatch(heatmap.2(t(heatMatrix66),Rowv=NA,
                    na.rm=TRUE,
                    na.color=getOption("mPhen.nacolor","black"), 
                    Colv=dendro, col = col,dendrogram = "column",
                    density.info="none",trace="none",xlab=xlab,ylab=ylab,key=T,margins=c(8,8),
                    RowSideColors = c(cols66),scale="none",cexRow = cexRow,cexCol=cexCol,main = main,
                    symbreaks = symbreaks),error = function(e) NULL)
    #legend(x="topleft",legend="-log(p)")
     if(!is.null(p3)){

      #  print(heatMatrix66)
    if(symbreaks) legend(x="topleft",legend="beta")
    else legend(x="topleft",legend="-log(p)")
    #legend(x="bottomleft",legend="1a",bty="n",cex=2.5)
    if(k==1)dendro = p3$colDendrogram
    }else{
#	print(heatMatrix66[1:5,1:5])	
	heatmap.2(t(heatMatrix66),Rowv=NA,
                    na.rm=TRUE,
                    na.color=getOption("mPhen.nacolor","black"), 
                    Colv=NA, col = col,dendrogram = "none",
                    density.info="none",trace="none",xlab=xlab,ylab=ylab,key=FALSE,margins=c(8,8),
                    RowSideColors = c(cols66),scale="none",cexRow = cexRow,cexCol=cexCol,main = main,
                    symbreaks = symbreaks)

#	print(dendro)
	warning(paste('unable to make heatmap',file[1], main[1]))	
   }
  }
  invisible(dendro)
}

##fingerprint plot
.fprint<-function(genD,
                  Res,
                  opts,
                  snpinds = 1:(dim(genD)[[2]]),
                  phenoinds =1:(dim(Res$Results)[[3]]), 
                  plot=TRUE){
  colsc = 2  
  #colsc = opts$mPhen.colsc
  res = Res$Results[,,phenoinds,,drop=F]
  maf = Res$maf
  nosnp = length(snpinds)
  geno_header = dimnames(genD)[[1]]
  if(length(dim(genD))>2) genoData=.getAvgAll(genD)   else genoData = genD
  li = dimnames(res)[c(1,2,3)]
  li[[2]] = dimnames(res)[[2]][snpinds]
  li1 = dimnames(res)[c(1,3)]
  li[[4]] = geno_header
  li1[[3]] = geno_header
  sumBetaX = .getArray(li1,v=0)
  intens = .getArray(li,v=0)
  dm = dim(sumBetaX)
  pvs = as.vector(res[,snpinds,,2])
  inds = which(pvs < 0.05)
  if(length(inds)==0){
    sig = colsc*max(abs(as.vector(res[,snpinds,,1])),na.rm=T)  # 
  }else{
    sig = colsc*max(abs(as.vector(res[,snpinds,,1]))[pvs<0.05],na.rm=T)  # 
  }
  sig = max(0.1,sig) 
  colsall = c()

  for(k1 in snpinds){
    k = snpinds[k1]

    nonNAind = !is.na(genoData[,k])
    gx = genoData[nonNAind,k]
    levs = levels(as.factor(gx))
    for(i in 1:dm[1]){
      for(j in 1:dm[2]){
        beta = res[i,k,j,1]
        if(is.na(beta)) beta =0
        sumBetaX[i,j,nonNAind] = beta*gx + sumBetaX[i,j,nonNAind] 
        intensj  =2*(pnorm(as.numeric(levs)*abs(beta),0,sig)-0.5)
        colj = if(beta<0) rgb(0,intensj,0) else rgb(intensj,0,0)
        colsgx = rep(0,length(gx))
        for(jj in 1:length(levs)){
          indjj = gx==levs[jj]
          colsgx[indjj] = jj + length(colsall)
        }
        intens[i,k1,j,nonNAind] = colsgx 
        colsall = c(colsall,colj)            
      }
    }
  }
  resuall = list()
  m = 1
  dmn = dimnames(sumBetaX)
  for(i in 1:dm[1]){
    for(j in 1:dm[2]){
      risk  = sumBetaX[i,j,]
      betas = res[i,snpinds,j,1]
      pvs = res[i,snpinds,j,2]
      ord = rev(order(risk))
      gt = betas>0
      lt = betas<=0        
      ordy = c(which(lt)[(order(pvs[lt]))], which(gt)[rev(order(pvs[gt]))])
      intensij = intens[i,,j,] 
      if(is.null(dim(intensij))) intensij = matrix(intensij,ncol = dim(genD)[[1]],nrow = length(snpinds))
      resu = list(risk =risk,cols=colsall,intens = intensij,
                  ord = ord, ordy = ordy,
                  ynme = dimnames(genoData)[[2]][snpinds], xnme = dimnames(genoData)[[1]], nme =paste(dmn[[1]][i], dmn[[2]][j],sep="_"))
      resuall[[m]] = resu
      m = m+1
    }

  }
  if(plot) {
    .plotFingerprint(resuall)
  }
  invisible(resuall)
}



.plotFingerprint<-function(resuall){
for(i in 1:length(resuall)){
resu = resuall[[i]]
nmes = resu$ynme
risk = resu$risk
intens = resu$intens
cols = resu$cols
ord = resu$ord
ordy = resu$ordy
nmey = resu$ynme
dm = dim(intens)
    x = t(intens)[ord,ordy,drop=F]
    dmx = dim(x)
 if(length(nmes)==dmx[2]) {
 par(mfrow = c(2,1))
    image(1:dmx[1],1:dmx[2],x,col = cols,main = resu$nme,xlab="indiv", ylab="snps",col.axis='white')
   axis(4,1:dmx[2], nmes,cex.axis = 0.5)
    plot(1:dmx[1], exp(risk[ord])/(1+exp(risk[ord])))
 }
}
}




##USED IN MAKING BED



.makeColsForBed<-function(){
red = c(255,0,0)
blue = c(0,0,255)
maxCN=4
minCN=0
cols = list()
for(i in minCN:maxCN){
  p = (maxCN-i)/(maxCN-minCN)
  colv = red * p +  blue * (1-p) 
  cols[[i+1]] = paste(colv,collapse=",")
}
cols[[1]] =  "250,218,221"
cols[[2]] = "255,0,0"
cols[[3]] = "0,0,0"
cols[[4]] = "0,255,0"
cols[[5]] = "0,0,255"

names(cols) =minCN:maxCN
cols
}

.getCNProfile<-function(row,starts, end,.base,cols,nme,chr){
 start = starts[1]
 row0 = row[2:length(row)]
 row1 = row[1:(length(row)-1)]
 br = which(row0!=row1)
 res = array(0,dim = c(1+length(br),9))
 st=1
   col = cols[[round(1+row[st] +.base)]]
 cn = rep(0,1+length(br))
 if(length(br)>0){
  for(k in 1:length(br)){
   en = br[k]+1 
   end1 = starts[en]
   cn[k] = round(row[st] )
   res[k,] = (c(chr,start,end1,nme,1000,'+',  start,end1,col))
   st = en
   start = starts[st]
   col = cols[[1+round(row[st] +.base)]]
 }
 }

  res[length(br)+1,] = (c(chr,start,end,nme,1000,'+',  start,end,col))
  cn[length(br)+1] = round(row[st])
 res[cn!=0,,drop=F]
}



###END PLOTTING ##


####COMMANDS FOR MULTIPHEN-BIGP ###

.orthogonalise<-function(pheno){
W = diag(rep(1,dim(pheno)[2]))
pheno1 = array(dim = dim(pheno))
pheno1[,1] = pheno[,1]
for(k in 2:(dim(pheno)[[2]])){
  ind = apply(pheno[,1:k,drop=F],1,.cntNa)==0
  W[1:(k-1),k] = -.projOut(pheno[ind,k,drop=F],P = pheno[ind,1:(k-1),drop=F])
  pheno1[,k] = pheno[,1:k] %*% W[1:k,k]
 }
dimnames(pheno1) = list(dimnames(pheno1)[[1]],paste(dimnames(pheno)[[2]],"orth",sep="_"))
list(pheno = pheno1, W = W)
}

#.projOutRem<-function(Dall,P){
#   W = .projOut(Dall, P = P)
#  R = P %*% W
#  return(Dall - R)
#}



### END BIGP ###


## genoInputs is a list of connections from different datasets, obtained from mPhen.openGenoConnection.  If there is more than one entry, 
##.batch is number of entries to read i
##.baseValue is the base genotype value, usually 0 for genotypes and 2 for copy number
## the results will be merged
## returns NULL if no lines of input are left
## .thin is used to only select a subset of lines. If value is 10, then only every 10th line included
###Returns a list, which includes $genoData, which is a matrix of genotypes
.readGenoConnection<-function(genoInput,opts=mPhen.options("geno.input")){
  .batch = opts$mPhen.batch
  .baseValue = opts$mPhen.baseValue
  .format = opts$mPhen.format
  .starting = opts$mPhen.starting
  .ending = opts$mPhen.ending
  .thin = opts$mPhen.thin
  rsidToDo = opts$mPhen.rsidToDo
  if(is.null(rsidToDo)) rsidToDo = attr(genoInput,"rsid") 
  if(!is.null(genoInput$genoFiles)) genoInput = list(genoInput)
  res = list()
  for(k in 1:length(genoInput)){
    type = genoInput[[k]]$type
    if(type=="zip"){
      res[[k]] = (.readZipConnection(genoInput[[k]],opts, .batch=.batch,.baseValue=.baseValue,
                                     .format = .format,.starting = .starting,.ending = .ending,
                                     .thin = .thin,rsidToDo = rsidToDo))
    }
    else   
      if(type=="bed"){
        res[[k]] = (.readPlinkConnection(genoInput[[k]],opts, .batch=.batch,.baseValue=.baseValue,
                                       .format = .format,.starting = .starting,.ending = .ending,
                                       .thin = .thin,rsidToDo = rsidToDo))
    }
    else {
      res[[k]] = (.readVCFConnection(genoInput[[k]],opts, .batch=.batch,.baseValue=.baseValue,
                                     .format = .format,.starting = .starting,.ending = .ending,
                                     .thin = .thin,rsidToDo = rsidToDo,type = type))
    }
    if(length(res)<k) break
  }
  if(length(res)!=length(genoInput)){
    #	print(length(names(res)))
    #        print(length(names(genoInput)))
    #     warning(paste("names length mismatch", paste(names(res),collapse=".")
    #        ,paste(names(genoInput),collapse=".")))
    return(NULL)
  }
  names(res) = names(genoInput)
  if(length(res)>1){
    allGenoData = list()
    for(k in 1:length(res)){
      allGenoData[[k]] = res[[k]]$genoData
      ###if(res[[k]]$chrom!=res[[1]]$chrom) stop('chroms should match')
    }
    res[[1]]$genoData = .mergeFiles(allGenoData, markerCol=FALSE, inclrow = attr(genoInput,"inclrow"), 
                                    anyCol = !opts$mPhen.onlyCommon)
  }
  res[[1]]
}


#Closes genotype connection.  genoInput is obtained from mPhen.openVCFConnection
.closeGenoConnection<-function(genoInput){
 if(!is.null(genoInput$type))  genoInput = list(genoInput)

 for(k in length(genoInput)){
  if(genoInput[[k]]$type!="zip"){
   conn = genoInput[[k]]$con
   close(conn)
  }
 }
}

# Will merge across multiple .genoFiles (i.e. if you have two cohorts measured on overlapping markers, you wish to merge in-silico)
## Can currently read 'GT' and 'GL' fields
#Outputs the header and re-ordering to be consistent with ordering of 'indiv' (if provided')
.openGenoConnection<-function(genoFiles,opts, indiv = NULL){
  
   if(is.null(names(genoFiles))) names(genoFiles) = paste("cohort",1:length(genoFiles),sep="")
  usePosition = opts$mPhen.matchPosition

  info_ind = if(usePosition) 1 else 3
  if(!is.list(genoFiles)) genoFiles = list(genoFiles)
  if(!is.list(indiv) & !is.null(indiv)) indiv = list(indiv)
  res = list()
  snpids = NULL
  rsmatch = NULL
  sampleids = list()
  for(k in 1:length(genoFiles)){
   indivk = NULL
   if(!is.null(indiv)) indivk = indiv[[k]]
   bedf = paste(genoFiles[[k]],"bed",sep=".")
   if(length(grep('.zip$',genoFiles[[k]]))>0){
            res[[k]] = (.openZipConnection(genoFiles[[k]],indiv=indivk))
            snpids1 = (res[[k]]$snps[,res[[k]]$posi[info_ind]])
        if(k==1){
              		 snpids = snpids1
                 rsmatch = array(dim = c(length(snpids), length(genoFiles)), dimnames = list(snpids, genoFiles))

        }
                 rsmatch[,k] = match(snpids, snpids1)
   }
   else if(file.exists(bedf)) res[[k]] = (.openPlinkConnection(genoFiles[[k]],indiv=indivk))
   else res[[k]] = (.openVCFConnection(genoFiles[[k]],indiv=indivk))
   sampleids[[k]] = res[[k]]$sampleids
  }
  if(length(res)==1){
     attr(res[[1]],"sampleids")<- res[[1]]$sampleids
     attr(res[[1]],"inclrow")<- 1:length(res[[1]]$sampleids)
     return(res[[1]])
  }
  else{
    names(res) = names(genoFiles)
    onlyCommon = opts$mPhen.onlyCommon
     
      if(onlyCommon & !is.null(rsmatch)){
            print("matching rsids")
            attr(res,"rsid") = dimnames(rsmatch)[[1]][apply(rsmatch,1,.cntNa)==0]
            print(attr(res,"rsid"))
      }
      inclrow = .findInclusion(sampleids,match=FALSE)
      sampleidsall = sampleids[[1]]
      for(i in 2:length(sampleids)){
	   sampleidsall = c(sampleidsall,sampleids[[i]][inclrow[[i]]])
      }
   attr(res,"inclrow")<-inclrow
   attr(res,"sampleids") = sampleidsall
      return(res)
   }
}



### READING ZIP AND VCF FILES
.openVCFConnection<-function(.genoFile, indiv = NULL){
 sampids = c()
 type = 'txt'
 if(length(grep('vcf',.genoFile))>0) type = 'vcf'
 if(length(grep('impute',.genoFile))>0) type = 'impute'
  con<- file(.genoFile, 'r')
if(type=='impute'){
  if(is.null(indiv)) stop('must specify sample order to use impute format via the indiv variable')
  header = c("IMP-STATUS","ID","POS","REF","ALT",indiv)
}else{
  header="##"
  while(length(grep('^##', header))>0){
   header = readLines(con,n=1)[1]
 }
   header = strsplit(gsub("#","",header),'\t')[[1]]
}
   vcf  = .readHeader(header,type)
   if(!is.null(indiv) & type!='impute'){
	   index_sample_inGenoFile<-match(indiv,header)
	 #  index_sample_inGenoFile <- index_sample_inGenoFile[!is.na(index_sample_inGenoFile)]
   }else{
	   index_sample_inGenoFile<-(vcf$firstGenoIndex:length(header))
   }
   sampids = c(sampids,header[index_sample_inGenoFile])
   res =    list(con=con,genoFiles = .genoFile, header=header,index=index_sample_inGenoFile, sampleids = header[index_sample_inGenoFile],
 	posi = vcf$posi, formatIndex=vcf$formatIndex,sampleids = sampids,type = type, firstGenoIndex = vcf$firstGenoIndex)
   res
}


.openZipConnection<-function(genoFiles, indiv= NULL,zipsep="\t"){
  res = list() 
  geno_header = NULL
  SNPs_all = NULL
  marker = c()
  sampleids= NULL
  posi=NULL
  formatIndex = NULL
  header= NULL
  start=c()
  end = c()
  for(i in 1:length(genoFiles)){
   genoFile = genoFiles[[i]]
   SNPs <- read.table(unz(genoFile, "SNPS"),as.is=T,sep=zipsep,header=FALSE)
   Samples<- read.table(unz(genoFile, "Samples"),as.is=T,sep=zipsep,header=FALSE)
   Name<- read.table(unz(genoFile, "Name"),as.is=T,fill=T,sep=zipsep,header=FALSE)
   names(Samples) = grep('^$',Name[3,],value=TRUE,invert=TRUE)[1:dim(Samples)[2]]
   names(SNPs) = grep('^$',Name[2,] ,value=TRUE,invert=TRUE)[1:dim(SNPs)[2]]
   header1 = Samples[,1]
   vcf = .readHeader(Name[2,],'zip')
   SNPs[,vcf$posi[2]] = sub("chr","", SNPs[,vcf$posi[2]])
   start = c(start, SNPs[1,vcf$posi[1]])
   end = c(end, SNPs[dim(SNPs)[1],vcf$posi[1]])
    geno_header1 = grep('^$',Name[1,],invert=T,value=TRUE)
    marker = c(marker,rep(i,dim(SNPs)[[1]]))
    posi1 = vcf$posi
    formatIndex = vcf$formatIndex
    if(i==1){ 
         geno_header = geno_header1
         SNPs_all = SNPs
	 header = header1
         posi=posi1
    }
    else if(length(which(geno_header!=geno_header1))>0) stop('zip files need to match')
    else if(length(which(header!=header1))>0) stop('zip files need to match')
    else if(length(which(posi!=posi1))>0) stop('zip files need to match')
    else{
         SNPs_all = rbind(SNPs_all,SNPs)
    }
   }
   if(!is.null(indiv)){
 	   index_sample_inGenoFile<-match(indiv,header)
 	  # index_sample_inGenoFile <- index_sample_inGenoFile[!is.na(index_sample_inGenoFile)]
    }else{
 	   index_sample_inGenoFile<-(1:length(header))
    }
      ord = order(SNPs_all[,posi[1]])
      list(con = NULL, genoFiles = genoFiles, header =geno_header,snps = SNPs_all[ord,], sampleids =header[index_sample_inGenoFile],
  	posi = posi, formatIndex=formatIndex, index = index_sample_inGenoFile, type="zip", marker=marker[ord], startend =cbind(start,end))
}

#genoConnection <- setClass("genoConnection",  
#       slots = c(genoFiles="character", header="character",
#			snps="data.frame", sampleids = "character", posi="numeric",
#			formatIndex="integer", index ="integer",type="character",
#			marker="integer",startend="matrix"  ))
            


.openPlinkConnection<-function(root, indiv = NULL){
 sampids = c()
  bed.file = paste(root, '.bed', sep = '')   
  bed.file.size = file.info(bed.file)$size # the bed file size in bytes
  famtable = read.table(paste(root, '.fam', sep = ''),as.is=T) # number of individuals
  header = famtable[,1]
  sample.size = dim(famtable)[1]
  pheni = 5:(dim(famtable)[2])
  pheno = famtable[,pheni,drop=F]
  phenN = paste("Phenotype",(pheni-5),sep="")
  phenN[1] = "Sex"
  dimnames(pheno) = list(header,phenN)
  snp.size = ceiling(sample.size/4) # no. bytes storing the data for 1 SNP for all individuals; each byte stores 4 people
  n.snps = round((bed.file.size-3)/snp.size) # the total .bed file size is 3 plus size of ids x snps combo hence removing 3
  bim = read.table(paste(root,'.bim', sep = ''), header = FALSE, colClasses = c('character', 'character', 'character', 'character', 'NULL', 'NULL'))
#  snp.names = as.factor(bim[,2]) # the bim file has the snps names in the 2nd column
#  snp.names = as.vector(snp.names)
  bin.connection = file(bed.file, 'rb') # opens a connection with .bed file
  test.bytes = readBin(bin.connection, what = "raw", n = 3) # bytes 1 and 2 store infor about file type, byte 3 about array type
  if(!identical(as.character(test.bytes), c('6c', '1b', '01'))) {
    stop('BED file not a v0.99 SNP-major BED file, please re-encode the data as v0.99 SNP-major file')
  }
   if(!is.null(indiv)){
	   index_sample_inGenoFile<-match(indiv,header)
#	   #index_sample_inGenoFile <- index_sample_inGenoFile[!is.na(index_sample_inGenoFile)]
   }else{
	   index_sample_inGenoFile<-(1:length(header))
   }
   sampids = c(sampids,header[index_sample_inGenoFile])
   res =    list(con=bin.connection,genoFiles = root, header=header,index=index_sample_inGenoFile, sampleids = header[index_sample_inGenoFile],
 	posi = bim[,c(4,1,2)], formatIndex=NULL,snp.size = snp.size,n.snps = n.snps, sampleids=sampids, type="bed",pheno=as.matrix(pheno))
  res
}

.readLinesPlink<-function(bin.connection,n,snp.size,sample.size){
 genotype = array(NA, dim = c(sample.size,n))
 for(k in 1:n){
    r.bin.snp = readBin(bin.connection, what = 'raw', n = snp.size)
    bin.snp = matrix(as.numeric(rawToBits(r.bin.snp)), ncol = 2, byrow = TRUE)[1:sample.size,]
    genotype[,k] = bin.snp[,1] + bin.snp[,2] - 10 * ((bin.snp[,1] == 1) & (bin.snp[,2] == 0))
  }
  genotype[genotype == -9] = NA 
  genotype
}


###modified 29-11
.readPlinkConnection<-function(genoInput,opts,.batch=20,.baseValue=0,.format="GT",.starting = 0, .ending = 900000000,.thin = 1, rsidToDo = NULL){
  usePosition = opts$mPhen.matchPosition
  info_ind = if(usePosition) 1 else 3 
  subinds = genoInput$index
  sample.size = length(genoInput$sampleids)
  posi = NULL
  if(is.null(rsidToDo)){
    input = matrix(1)
    lastPos = -1
    kj =0
    while(dim(input)[2]>0 & lastPos < .starting & kj <genoInput$snp.size){
      toRead = min(.batch,genoInput$n.snps - kj)
      input<-.readLinesPlink(genoInput$con, n=toRead,snp.size = genoInput$snp.size,sample.size = sample.size)[subinds,,drop=F]
      if(.thin>1) input<-input[,seq(1,length(input),by=.thin),drop=F ]
      len<-dim(input)[2]     
      firstPos = genoInput$posi[kj+1,1]
      lastPos = genoInput$posi[kj+toRead,1]
      posi = genoInput$posi[(kj+1):(kj+toRead),]
      if(firstPos > .ending){
        input = NULL
        break;
      }
      kj = kj+toRead
    }
  }
  else {
    if(length(rsidToDo)==0) return(NULL)
    indsToDo = match(rsidToDo,genoInput$pos[,info_ind])
    indsToDo = sort(indsToDo[!is.na(indsToDo)])
    len<-length(indsToDo)
    posi = genoInput$posi[indsToDo,]
    input =array(NA, dim = c(sample.size,length(indsToDo)))
    for(kj in 1:(genoInput$n.snps)){
      indsI = which(indsToDo==kj)
      if(length(indsI)>0){
        input[,indsI[1]]<-.readLinesPlink(genoInput$con, n=1,snp.size = genoInput$snp.size,sample.size = sample.size)[subinds,,drop=F]
      }
      else {
        tmp<- readBin(genoInput$con, what = 'raw', n = genoInput$snp.size)
      }
    }
    firstPos = genoInput$posi[1,indsToDo[1]]
    lastPos = genoInput$posi[1,indsToDo[length(indsToDo)]]
  }
  nsamp = sample.size
  genoData = input
  dimnames(genoData)[[1]] = genoInput$sampleids 
  chrom=posi[,2]
  rsid = posi[,3]
  #pos_info = apply(posi,1,paste,collapse="\t")
  nbsnp = rep(NA,len)
  snpids=posi[,3]
  indsToKeep = 1:len
  dimnames(genoData)[[2]] = apply(posi[,2:1],1,paste,collapse="_")
  res1 = list(len=length(indsToKeep),nbsnp = nbsnp[indsToKeep], firstPos=firstPos,lastPos=lastPos, genoData=genoData[,indsToKeep,drop=F],
              chrom=chrom, rsids=posi[,info_ind])
  res1
}



###HAS BEEN MODIFIED 29-11
.readVCFConnection<-function(genoInput,opts, .batch=20,.baseValue=0,.format="GT",.starting = 0, .ending = 900000000,.thin = 1, rsidToDo = NULL, type = 'vcf'){
  sep_ = '\t'
  if (type=='vcf') split1 =  '[\t:]' 
  else if(type=='impute'){
    split1 = ' '
    .format = "IMP"
    sep_=' '
  }
  else split1 = '\t'
  vcf = type=='vcf'
  if(is.null(rsidToDo)){
    input = c(1)
    lastPos = -1
    while(length(input)>0 & (lastPos < .starting)){
      input<-readLines(genoInput$con, n=.batch)
      len<-length(input)     
      if( len==0) break
      if( is.null(input[1])) break
      if(.thin>1) input<-input[seq(1,length(input),by=.thin) ]
      splits=strsplit(input[1],split1)[[1]]
      firstPos = as.numeric(splits[genoInput$posi[1]])
      splits=strsplit(input[len],split1)[[1]]
      lastPos = as.numeric(splits[genoInput$posi[4]])
      if(firstPos > .ending){
        input = c()
        break;
      }
    }
  }
  else {
    if(length(rsidToDo)==0) return(NULL)
    for(k in 1:length(rsidToDo)) rsidToDo[k] = paste(rsidToDo[k],sep_,sep="")
    input = rep("",length(rsidToDo))
    done = rep(FALSE,length(rsidToDo))
    while(length(which(!done))>0){
      input1<-readLines(genoInput$con, n=.batch)
      if(length(input1)==0){
        break;
      }
      notd = which(!done)
      for(k in notd){
        subl = grep(rsidToDo[k], input1,value=TRUE)
        if(length(subl)>0){
          input[k] = subl[1]
          done[k] = TRUE
        }
      }
    }
    input = input[which(done)]
    len = length(input)
    splits=strsplit(input[1],'\t')[[1]]
    firstPos = as.numeric(splits[genoInput$posi[1]])
    splits=strsplit(input[len],'\t')[[1]]
    lastPos = as.numeric(splits[genoInput$posi[4]])
    #rsidToDo = c()   #rsidToDo[!done]  
  }
  nsamp = length(genoInput$sampleids)
  if(length(input)==0) return(NULL)
  if(is.null(input[1])) return(NULL)
  chrom=splits[genoInput$posi[2]]
  rsids = rep("",len)
  nbsnp = rep(NA,len)
  snpids=rep("",len)
  posi = genoInput$posi
  formatInd = 1
  .split = .getVCFSplit(.format[1])
  .process = .getVCFProcess(.format[1])
  dn = .getVCFDimNames(genoInput$sampleids,snpids,.format[1])
  genoData = .getArray(dn[[1]])
  formatIndex = genoInput$formatIndex
  firstIndex = if(vcf) formatIndex else formatIndex +1
  if(type=='impute') firstIndex = genoInput$firstGenoIndex
  nsamp1 =  length(genoInput$header)-(firstIndex-1)
  chrpos = rep(0,len)
  subinds = genoInput$index
  dimg = dim(genoData)
  #if(length(formatIndex)==1)
  for (i in 1:len){
    line = strsplit(input[i],split1)[[1]]
    formatLen = (length(line) - (firstIndex-1))/(nsamp1)
    rsids[i] = line[posi[3]]  #paste(line[posi[1:3]],collapse='\t')
    if(opts$mPhen.matchPosition) rsids[i] = line[posi[1]]
    chrpos[i] = as.numeric(line[posi[1]])
    if(posi[2]<0) {
      fileline =  strsplit(genoInput$genoFiles,"/")[[1]]

      chrom_ = strsplit(fileline[length(fileline)],"\\.")[[1]][1]
      snpi = paste(chrom_,line[unique(posi[1])],sep='_')
    }
    else snpi = paste(line[unique(posi[c(2,1)])],collapse='_')
    snpids[i] = snpi 
    if(posi[5]>=0) nbsnp[i]  =  as.numeric(line[posi[5]])
    if(type=='impute'){
      line1 = 1000*t(apply(matrix(line[firstIndex:length(line)],ncol = 3,byrow = TRUE),c(1,2),as.numeric))
    }
    else {
      if(formatLen>1){
        .formatV = line[seq(formatIndex,formatIndex+formatLen-1)]
        formatInd = which(.formatV==.format)
        if(length(formatInd)==0){
          warning('did not find format id')
          next
        }
        subinds = (genoInput$index - firstIndex)*formatLen+firstIndex+(formatInd-1)
      }
      if(.split==""){
        line1 = as.numeric(line[subinds])
      } 
      else {
        m = as.matrix(data.frame(strsplit(line[subinds],.split)))
        line1 = apply(m,2,.process)
      }
    }
    if(length(dimg)>2) genoData[,i,] = t(line1) else genoData[,i] = line1 - .baseValue
  }
  indsToKeep = which(chrpos>=.starting & chrpos<=.ending)
  firstPos = chrpos[indsToKeep[1]]
  lastPos =chrpos[indsToKeep[length(indsToKeep)]]
  print(paste("first","last",firstPos,lastPos))
  dimnames(genoData)[[2]] = snpids
  list(len=length(indsToKeep),nbsnp = nbsnp[indsToKeep], firstPos=firstPos,lastPos=lastPos, 
       genoData=if(length(dimg)>2) genoData[,indsToKeep,,drop=F] else genoData[,indsToKeep,drop=F] 
       ,rsids=rsids[indsToKeep],chrom=chrom)
}


##reads full zip file
.readZip<-function(genoFile,zipsep="\t"){
  SNPs <- read.table(unz(genoFile, "SNPS"),as.is=T,sep=zipsep,header=FALSE)
  Samples<- read.table(unz(genoFile, "Samples"),as.is=T,sep=zipsep,header=FALSE)
  Name<- read.table(unz(genoFile, "Name"),as.is=T,fill=T,sep=zipsep,header=FALSE)
  types = grep('^$',Name[1,],invert=T,value=T)
   names(Samples) = grep('^$',Name[3,],value=TRUE,invert=TRUE)[1:dim(Samples)[2]]
   names(SNPs) = grep('^$',Name[2,] ,value=TRUE,invert=TRUE)[1:dim(SNPs)[2]]
  res = array(dim = c(dim(Samples)[[1]],dim(SNPs)[[1]],length(types) ),dimnames = list(Samples[,1],SNPs[,4],types))
  res = res[,,]
  len = length(dim(res))
  for(k in 1:(dim(SNPs)[[1]])){
      t = read.table(unz(genoFile, SNPs[k,4]),as.is=T,fill=T,sep=zipsep,header=FALSE)
      if(len==2) res[,k] = t[,]
      else res[,k,]  = t 
  }
  res
}


.getVCFSplit<-function(.format){
 .split = rep("",length(format))
 for(k in 1:length(.format)){
    if (.format[k] == "GT"){
     .split[k] =   "[\\|\\/]"
   }else if (.format[k] == "GL"){
     .split[k] =  ","
   }else if(.format[k] =="IMP"){
     .split[k] =  ","
   }else if(.format[k] =="AD"){
     .split[k] =  ","
   }else if(.format[k] =="distr"){
     .split[k] =  "[;=]"
   }
  }
  .split
}

.sum1<-function(vec) sum(as.numeric(vec))

.getVCFProcess<-function(.format){
 # .process = rep(NA, length(format))
 k = 1
# for(k in 1:length(.format)){
    if (.format[k] == "GT"){
      .process =.sum1
   }else if (.format[k] == "GL"){
      .process = .pow
   }else if(.format[k] =="IMP"){
      .process = .splitImputed
   }else if(.format[k] =="AD"){
      .process = .unchanged
   }else if(.format[k]=="distr"){
     .process = .cnvdistr
   }else{
	.process = .unchanged
   }
 #}
.process
}

 .getVCFDimNames<-function(sampleids, snpids,.formats){
   dn = vector("list",length(.formats))
   for(k in 1:length(.formats)){
     .format = .formats[k]
    if(length(which(c("GL","IMP")==.format))>0){
   ### NOTE, ASSUMES 3 GENOTYPES HERE
     dn[[k]] = list(sampleids,snpids,c(0,1,2))
   }else if("distr"==.format){
   ### NOTE, ASSUMES 5 CNV GENOTYPES HERE
     dn[[k]] = list(sampleids,snpids,c(0,1,2,3,4))
   }else if(length(which(c("AD")==.format))>0){
   ### NOTE, ASSUMES 3 GENOTYPES HERE
     dn[[k]] = list(sampleids,snpids,c(0,1))
   }else{
        dn[[k]] = list(sampleids,snpids)
   }
   }
dn
}

##modified 1/2014
.readZipConnection<-function(genoInput,opts, .batch=10,.baseValue=0,.format="geno",.starting = 0, .ending = Inf,.thin = 1, rsidToDo = NULL, zipsep="\t"){
  genoFile = genoInput$genoFiles
  posi = genoInput$posi
  if(length(genoInput$header)==1) genoInd = 1 
  else{
    genoInd = .findIn(genoInput$header,.format,critical=TRUE)
  }
  .format = genoInput$header[genoInd]
  snps = genoInput$snps
  pos =  snps[,posi[1]]
  rsids = snps[,posi[3]]
  if(is.null(rsidToDo) ){
    todo =  which(pos>=.starting & pos<=.ending)
    if(length(todo)==0){
      return(NULL)
    }
    if(.thin>1) todo<-todo[seq(1,length(input),by=.thin) ]
    todo = todo[1:min(.batch,length(todo))]
  }
  else {
    if(opts$mPhen.matchPosition) tomatch = pos else tomatch = rsids
    todo = match(rsidToDo,tomatch)
    todo = todo[!is.na(todo)]
    todo = todo[pos[todo]>=.starting & pos[todo]<=.ending]
    if(length(todo)==0){
      return(NULL)
    }
  }
  marker = genoInput$marker[todo]
  #startend = genoInput$startend
  pos = pos[todo]
  snps = snps[todo,,drop=F]
  rsids = rsids[todo]   
  len = length(todo)
  if(len==0){
    next     
  }
  firstPos = pos[1]
  lastPos =  pos[len]
  nsamp = length(genoInput$sampleids)
  input = array(NA, dim = c(len,nsamp))
  todrop = rep(FALSE, length(todo))
  for(j in 1:length(todo)){
    rsidj = rsids[j]
    if(length(grep(':',rsidj))>0){
      t1 = system(paste("unzip -p",genoFile[marker[j]],rsidj),intern=TRUE)
      t1 = t(as.matrix(data.frame(strsplit(t1,"\t"))))
    }
    else t1 = read.table(unz(genoFile[marker[j]], rsidj),as.is=T,fill=T,sep=zipsep,header=FALSE)
    if(dim(t1)[2]>0)  input[j,]<- t1[genoInput$index,min(genoInd, dim(t1)[2])]	
    else todrop[j] = TRUE
  }
  input = input[,!todrop,drop=FALSE]
  pos = pos[!todrop]
  snps = snps[!todrop,,drop=F]
  rsids = rsids[!todrop]   
  chrom=snps[1,2]
  pos_info = rep("",len)
  nbsnp = rep(NA,len)
  snpids=rep("",len)
  formatInd = 1
  .split = .getVCFSplit(.format[1])
  .process = .getVCFProcess(.format[1])
  dn = .getVCFDimNames(genoInput$sampleids,snpids,.format[1])
  formatIndex = genoInput$formatIndex
  chrpos = rep(0,len)
  #subinds = genoInput$index
  #if(length(formatIndex)==1)
  genoData = .getArray(dn[[1]])
  dimg = dim(genoData)
  for (i in 1:len){
    line = input[i,]
    line_ = snps[i,]
    formatLen = 1
    pos_info[i] = line_[[posi[3]]] #paste(line_[posi[1:3]],collapse='\t')
    if(opts$mPhen.matchPosition) pos_info[i] = line_[[posi[1]]]
    chrpos[i] = as.numeric(line_[posi[1]])
    snpids[i] = paste(line_[unique(posi[c(2,1)])],collapse='_')
    if(posi[5]>=0) nbsnp[i]  =  as.numeric(line_[posi[5]])
    if(.split==""){
      line1 = as.numeric(line)
    }
    else {
      line2 = sub(',,',',0,',sub("^,","0,",sub(",$",",0",line)))
      line1 = as.matrix(data.frame(lapply(strsplit(line2,.split),.process)))
      #	   line1 = apply(m,2,.process)
    }
    if(length(dim(genoData))>2) genoData[,i,] = t(line1) else genoData[,i] = line1 - .baseValue
  }
  print(paste("first","last",firstPos,lastPos))
  #   if(length(unique(snpids)) < length(snpids)) snpids = pos_info
  dimnames(genoData)[[2]] = snpids
  list(len=length(todo),nbsnp = nbsnp, firstPos=firstPos,lastPos=lastPos, genoData=genoData,rsids=pos_info,chrom=chrom)
}

## END ZIP FILES

##Carries out a meta-analysis on a list of pvalue tables
##pvs is a list of tables. 
##pvs can be obtained from output$pvs where output is obtained from mPhen.writeToOutputConnection
##NEED TO UPDATE THIS WITH BETAS
.metaAnalysis<-function(pvs, betas){
#   print(names(pvs))
  arr = .getArray(list(names(pvs), pvs[[1]][,3],dimnames(pvs[[1]])[[2]][9:11]))
  for(j in 1:length(pvs)){
    arr[j,,] = as.matrix(pvs[[j]][,9:11])
  }
  results = pvs[[1]]
  for(j in 1:(dim(results)[[1]])){
  results[j,9:11] = .metaresInvVarianceFixed(arr[,j,])
 }
 results
}


.updateRelatedness<-function(geno, K,
na.replace=0,standardise = TRUE){
 if(standardise){
	 geno1 = apply(geno,2,.standardise,na.replace)
 }else{
	 geno1 = apply(geno,2,.centralise,na.replace)
 }
 K = K + geno1 %*% t(geno1)
 K
 }



.proc <-
function(vec) (vec - 0.5)*vec*(1-vec)
.proc1 <-
function(vec) vec





.Random.seed <-
c(403L, 10L, 1828790116L, -106126681L, 1781369853L, 1907359220L, 
1618369242L, 575036453L, -2100244161L, 551208646L, -1305195696L, 
1620063091L, 1467360657L, 1338402752L, -407563874L, 1163320169L, 
-255471717L, -158550358L, -1812528948L, 870089103L, -2106141243L, 
211777916L, 1557350866L, 765670733L, 772941831L, -102791186L, 
1511796424L, -922664501L, 923776041L, -1944870280L, 3379110L, 
748652993L, -993182125L, 896554242L, -1217702860L, -1306090057L, 
2099857837L, 407210692L, -1737355542L, 1201774933L, -538159377L, 
-553786122L, -1913441248L, -846060061L, -13789567L, -1530732304L, 
2066358158L, -1000476423L, -1517618549L, -808068422L, -434104900L, 
366411519L, -176328107L, 1197720940L, 541494722L, -1351586275L, 
-1199638889L, 1595621374L, -2062865480L, 883558299L, 48508409L, 
1971645768L, 1890658422L, -625486415L, 391830499L, 1401219250L, 
1693943044L, 627410503L, 826265885L, 836231060L, -1567440070L, 
1936527237L, 799186015L, 2111821542L, -798433936L, 214460051L, 
-619888847L, -1104476704L, 1268081278L, 975065353L, -227675333L, 
-1163343094L, 1183580908L, 742096751L, 752066725L, -1886844324L, 
-1424381390L, -1714749139L, -2044360729L, 1780256526L, 696557800L, 
-1060023573L, -1304265847L, 421308184L, 1772530246L, 1173047969L, 
436612403L, 907048098L, 510622612L, -99073513L, 588695053L, -21142748L, 
-186884342L, 1021375093L, 935653839L, -943932714L, -1855201664L, 
987873475L, -148496799L, -268937520L, -463310802L, -2094795303L, 
216467947L, 1178714L, -1310595556L, -220788961L, 256547701L, 
1313866252L, 58938338L, 1657077949L, 562966327L, 1956354014L, 
-2133269224L, -850507653L, -354292199L, -873658712L, 1861354582L, 
629260497L, -350462845L, -1155268846L, -1692808540L, -404089241L, 
-509179971L, 180214068L, -80487782L, -1605232411L, -1124676993L, 
548295686L, -604797040L, -1442723277L, -1314774575L, 271139712L, 
1914543710L, -1870977623L, -438590373L, -1825637270L, 1091986700L, 
1618650959L, -1953541115L, -2011927876L, 2141445010L, -718038899L, 
352493255L, -1450838098L, -98934264L, -1350973941L, -247853335L, 
-1708947016L, 1573512806L, -808721663L, 164767123L, -1373274558L, 
-979254796L, -579232777L, -236740627L, 873048196L, -135508310L, 
-608639979L, -478874449L, -392725834L, -1995065120L, -1945648093L, 
-1487013055L, -1015814608L, -2053000114L, -133243719L, 2117245131L, 
-1190533766L, 233822332L, -392090433L, -1893832299L, -859732L, 
731091842L, 1594415453L, 2098780247L, -309293250L, 213235576L, 
-390931365L, -2105439303L, 210172424L, -493133642L, 1119786993L, 
-2109586909L, 1614735986L, -2142539580L, -1736042873L, -1653567907L, 
-40243116L, -757099910L, 1442687429L, 962288415L, 220288678L, 
-2095527888L, 1185323731L, -1800736527L, -1510212576L, -881038658L, 
723377865L, 1452063867L, -1720675766L, 899458732L, 975348143L, 
-799651483L, -371927268L, -1924834446L, -1738934291L, -261392345L, 
-183861170L, -1935317848L, -1367930197L, 259866569L, -2025516072L, 
-946851194L, 702964577L, -115480845L, 915970402L, -1678104620L, 
-641232169L, 1432828109L, 1588173988L, 1355174584L, 378907482L, 
-244912752L, 1971642748L, -1715249900L, -1916810942L, -1587313920L, 
-1220977732L, -1079637168L, 1634935074L, -784724664L, -818554868L, 
-558778276L, -1417685806L, 1832355840L, -980579500L, -1375424024L, 
-263313414L, -1682839856L, 1755518252L, -82126604L, 1581143954L, 
1969654880L, 993758476L, 601548512L, 1731833186L, 117244584L, 
1848158796L, -956102788L, -1030862030L, -903672208L, -1188405948L, 
1808870360L, 2087575674L, -1421879760L, -1740449092L, -1466861740L, 
-257000894L, 1729813408L, -732276548L, 122778064L, -1493219998L, 
1270047720L, 22874060L, 1754171132L, -176476174L, 941980800L, 
-1432786700L, -1779190648L, -817253062L, -1883757744L, 2026843244L, 
862097556L, -1038128814L, -1826783520L, -1779121972L, -1426401120L, 
226531714L, -1983538968L, 271255788L, 1276282940L, 133980274L, 
-2112370576L, 1073144356L, -1525611976L, 452794906L, -762478384L, 
-123566532L, -615830380L, 1010746498L, -1723594816L, -57639620L, 
-1773605232L, -1395129310L, 602992648L, 394058892L, -1480273892L, 
-1999743854L, 1335471232L, 1998210068L, -1462066392L, 400266298L, 
844926160L, -601648020L, 770639540L, -280604654L, -1108058144L, 
1195748300L, 278093280L, -1909482462L, -530242264L, 1150146828L, 
-1894519236L, 646722802L, -1195464464L, 238741316L, 692441176L, 
1483227706L, 847804656L, -1144958276L, -2075887084L, 424504514L, 
9266784L, -1685228164L, -1098451760L, -134339550L, 347072680L, 
-857289332L, 2007257532L, 415012914L, -594705600L, -409137484L, 
1772183880L, 1236640442L, -582516976L, 1160601836L, -1628886956L, 
330368786L, -196009568L, -2121275444L, 1424031584L, -771990334L, 
425586216L, 965286316L, -794836164L, -2090478734L, 1202966320L, 
-1045405660L, 1558490680L, -194686246L, -1497729136L, -619645700L, 
1253571988L, -991743166L, -1176811008L, -290715588L, 1416933968L, 
-100274014L, -1520868536L, 2045348364L, 1119746780L, -1166135214L, 
503532672L, 63899476L, -1927794712L, -571833606L, 1668007760L, 
1823646764L, 1678170484L, 452791954L, 1791140448L, -4738932L, 
68628832L, -461760286L, 1481641128L, 1694287564L, 874420732L, 
1892891570L, 631243760L, -1714925116L, 1769721176L, -197320326L, 
-1398648400L, -1820137028L, 1889386068L, -1925604926L, 1583604512L, 
-312806468L, 815820880L, -229770270L, 872456040L, -834504884L, 
2018127612L, 440350450L, 520568064L, 972582004L, 1854803464L, 
-2003409478L, 1514501584L, 1426307052L, -910906476L, 2125056210L, 
2008925280L, 570740044L, 272438048L, -1205822846L, -1482980888L, 
-526509076L, -88637892L, -621902222L, -174847120L, -297124572L, 
1951472184L, -529813094L, 877096400L, -1678449220L, -2007175532L, 
1818272386L, -1837548096L, -2034690244L, -967765360L, -1418315998L, 
1563948552L, -1987861748L, -1613634660L, 2039700498L, 992816000L, 
675248916L, -820447960L, -823685062L, 84902864L, -1944880532L, 
1402227764L, -1540042606L, 954366304L, -1444373172L, 2021266912L, 
242326434L, 27900840L, -267903860L, -1647070276L, -2129304462L, 
-530438160L, -2044279740L, 1171236952L, 1201823802L, -517806928L, 
21230286L, -1107269637L, -2020904867L, -2019846262L, -461367384L, 
282926353L, -172383229L, 60722948L, -565042382L, -964116089L, 
-984471599L, 1511098710L, -1235720396L, 830425669L, 1796781711L, 
1142185992L, -2020883418L, 804019411L, -727269835L, -273582702L, 
-210400928L, -1091948343L, 287851515L, -2123291700L, 1447288410L, 
-1837317297L, 190667289L, -1850140626L, -143960516L, 2070916013L, 
-1911946153L, 1420489120L, 2105019134L, -132044341L, -934178867L, 
312902394L, -846268104L, 1263157985L, 634307987L, 717565364L, 
1748872898L, 1600356055L, 1723985697L, -1694742810L, 500908004L, 
626061077L, 912869183L, 1044533976L, -1902405898L, 149094659L, 
-1990458811L, -1699343774L, 1583550928L, 1142894329L, -1792086997L, 
1685564316L, -951172150L, -1268319041L, -540680695L, -1584359650L, 
740105932L, -1864180067L, -1232757977L, -62818992L, -1563362322L, 
-902871909L, -1572360451L, 1707505130L, -697727352L, -726151119L, 
-2016236893L, 927011492L, -1638071470L, 833915559L, 977008369L, 
-207858122L, -1080575276L, 868634213L, 1657633903L, -2107341144L, 
-1177123834L, -1523216077L, 1720618389L, 313238002L, 807374400L, 
264826793L, -747458021L, -732935188L, 923077626L, 2096850991L, 
-1278121927L, -1579550130L, -768779620L, -259637235L, -1779196809L, 
-264848256L, 12188382L, -1761859157L, 1233478317L, 1828206362L, 
-8591400L, -362781631L, 576762611L, 1649702164L, -1608671070L, 
1424555447L, 710046081L, 1055770502L, 1793762756L, -290458379L, 
-484868769L, 130265016L, 1929831062L, -1928737629L, -894728603L, 
-1228512382L, -1027522320L, -924211047L, -1085947765L, 1221929340L, 
1409980330L, 52873183L, 1444250985L, -1255106178L, -1627494036L, 
2118521277L, 437674631L, 227391088L, 1147608334L, 1772855867L, 
-2097583843L, 1388436426L, 1643812072L, 1497319377L, -1988175421L, 
125797060L, 1876795762L, 208817223L, -1378934767L, -436647530L, 
118074612L, 747079301L, 425430991L, -1115038008L, -1820453018L, 
95047059L, 1445646069L, -913502510L, 982776096L, -208234231L, 
90951227L, 1463313292L, 692838938L, 553525903L, 1327017433L, 
-1959169426L, 1006297596L, -1512598675L, 1061666071L, 1340897504L, 
771710654L, 10914955L, -1193645299L, 863370042L, 2144383864L, 
599515425L, 1357515219L, 603994356L, 1176567778L)

.tabulateWithWeights <- function(phen,vars,weights,totweight){
  tab = table(list(vars,phen))
  nme = dimnames(tab)
  nme1 = as.numeric(nme[[1]])
  nme2 = as.numeric(nme[[2]]) 
  tab1 = matrix(0,nrow = length(nme1), ncol  = length(nme2))
  for(i in 1:length(nme1)) 
    for(j in 1:length(nme2)) 
      tab1[i,j] = sum(weights[vars==nme1[i] & phen==nme2[j]])
  tab2 = tab1/totweight
  tab2
}


