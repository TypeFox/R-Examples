###samples with mean p value of meanp, with mu controlling the variance
###via a beta distribution.  Mu>0, with higher values indicating tighter
## distribution.  Imputed controls whether to generate imputed data, and dirichlet controls how much noise in this data.
mPhen.sampleGeno<-function(n = 100, sampSize = 100, chr="0",pos = 1:n, snpids = paste(chr,pos,sep="_"),meanAlleleFreq=0.2, mu = 10,samples =paste("id",1:sampSize,sep="_"),imputed = FALSE, dirichlet = 1){
  genos = 0:2

  sampSize = length(samples)
  mean = meanAlleleFreq
  if(imputed){
    geno = array(dim = c(sampSize,length(snpids),length(genos)),dimnames = list(samples,snpids,genos))
  }else{
    geno = array(dim = c(sampSize,length(snpids)),dimnames = list(samples,snpids))
  }
  p = rbeta(length(snpids),mu*mean, mu*(1-mean))
 for(i in 1:length(snpids)){
    probs = c(p[i]^2,2*p[i]*(1-p[i]),(1-p[i])^2)
    g = sample(genos,sampSize,prob = probs,replace=TRUE)

    if(imputed){
     for(k in 1:sampSize){
       v = rep(1,length(genos))
       v[g[k]+1] = v[g[k]+1] + 1000
       geno[k,i,] = 1000*.sampDirichlet(v,dirichlet)
      }
    } else{
       geno[,i] = g
   } 
 }
geno
}
## generates correlation which is in a block format. Order is between, within(on-diagonal), within (off-diagonal)
mPhen.sampleCovar<-function(noPhenos,blockSize, orthogAll = c(0.9,0.5),dirichletScale = 50,  resample = FALSE, sd = rgamma(noPhenos,shape=10,rate = 10)){
      blockSize = blockSize
      noBlocks = ceiling(noPhenos/blockSize)
      if(noPhenos %% blockSize>0) stop('noPhenos needs to be multiple of block size for now')
        blocks = array(FALSE,dim = c(noPhenos,noBlocks))
        chols = list()
	M2 = diag(rep(1,noPhenos))
	sgnM = array(0,dim = c(noPhenos,noPhenos))
        covarB = .genRandCovar(noBlocks,orthogAll = orthogAll[1],dirichletScale = dirichletScale,resample = FALSE)
        cholB = chol(cov2cor(covarB))
        sgnB = sign(cholB)       
	cholB2 = cholB^2   
        Minner2 = array(0,dim = c(blockSize,blockSize))
	sgnInner = array(0,dim = c(blockSize,blockSize))

        for(k in 1:noBlocks){
            indsk = (k-1)*blockSize+(1:blockSize)
            blocks[indsk,k] = TRUE
            covark = .genRandCovar(blockSize,orthogAll = orthogAll[2],dirichletScale = dirichletScale,resample = FALSE)
            cholsk = chol(cov2cor(covark))
            sgnk = sign(cholsk)
            cholsk2 = cholsk^2
            M2[indsk,indsk] = cholsk2 * cholB2[k,k]
            sgnM[indsk,indsk] = sgnk
            if(k>1){
	      for(j in 1:(k-1)){
                  indsj = (j-1)*blockSize+(1:blockSize)
                  for(i in 1:blockSize){
			 Minner2[,i] = .sampDirichlet(rep(1,blockSize),weight=dirichletScale)*cholB2[j,k]
	                 sgnInner[,i] =  sample(c(-1,1),blockSize,replace=T)
                  }
                  M2[indsj,indsk] = Minner2
		  sgnM[indsj,indsk] = sgnInner
              }
  

            }
        }
        M = sqrt(M2) * sgnM ##squares every element    
       res = (t(M) %*% M) * outer(sd,sd)
       pnames = paste("pheno",(1:noPhenos),sep="")
       if(resample){ 
           ord = sample(1:noPhenos)
           res = res[ord,ord]
       }
       dimnames(res) = list(pnames,pnames)
        res
}




##Simulates multi-dimensional phenotype data
##x is a continuous variable (usually genotype data)
## sample_names are vector of  sample ids, and has same length as x
##corrSt is an object with cholesky decomposition, obtained from mPhen.cholesky
##varexp is a number between 0 and 1 representing the target variance explained in the major direction of effect
mPhen.simulate<-function(x,sample_names,covar,effDir, varexp, inverse = FALSE, geno.link="gaussian", 
 effDirInReverseEigenspace = FALSE, freq = 0.1){
 if( effDirInReverseEigenspace){
      ev = eigen(covar)
      beta = effDir/ sqrt(effDir %*% effDir)
      effDir<- (ev$vectors %*% rev(beta))[,1]
 }
  if(inverse){
     x = rep(0,length(sample_names))
     v = NULL
     varexp0=0
  }else{
     v = effDir
     y = x
     varexp0 = varexp
  }
   corrSt = .cholesky(covar,v)  ## gets a cholesky decomposition required for simulation
   cholT=corrSt$cholT;  
   Minv = corrSt$Minv;
   sdT = corrSt$sdT
   means = corrSt$means
   outpT = .sampJoint(cholT,x=as.vector(x),varexp = varexp0,sd=sdT)
   outp =means +  outpT %*% Minv  ##t(M %*% t(outpT))
   dimnames(outp) = list(sample_names,dimnames(cholT)[[1]]) 
 
  if(inverse){
    y = .sampJoint(matrix(1), x = (outp %*% effDir)[,1], varexp = varexp, sd = 1)
    if(geno.link=="binomial"){
	  prob = pnorm(y)
	    caseStatus = as.matrix(apply(cbind(1-prob,prob),1,.sample1))
            y = caseStatus[,1]
    }else if(geno.link=="ordinal"){
	thresh = c((1-freq)^2,2*freq*(1-freq),freq^2)        
        prob1 = pnorm(y-qnorm(thresh[1]))
        mat1 = cbind(1-prob1,prob1)
        prob2 = pnorm(y-qnorm(sum(thresh[1:2])))
        mat2 =  cbind(1-prob2,prob2)
        caseStatus = as.matrix(apply(mat1,1,.sample1))
        ind = caseStatus==1
	caseStatus[ind,] = (as.matrix(apply(mat2,1,.sample1))+1)[ind,]
        y = caseStatus
   }
 }
 limit = list(phenotypes = dimnames(outp)[[2]],covariates = NULL, resids = NULL, strats = NULL, excls = NULL)
 list(pheno = outp, geno = y,effDir=effDir, limit = limit)
}



## Plots phenotype data, coloured by corresponding genotype
##pheno is a matrix of phenotypes
mPhen.plotCorrelation<-function(pheno_to_plot,geno,title="",cex=0.25,cols = c(1,2,3)){
if(dim(geno)[2]>10) stop('too many genotypes to plot')
if(dim(pheno_to_plot)[2]<=1) return(NULL)
dims = rev(c(dim(geno)[2],ceiling(dim(pheno_to_plot)[2]/2)))
par(mfrow=dims)
for(i in seq(1,dim(pheno_to_plot)[2],2)){
 i1 = i
 if(i1==dim(pheno_to_plot)[2]) i1 = i-1
 for(j in 1:dim(geno)[2]){
  x = geno[,j]
  zero = which(x<0.5)
  one = which(x>=0.5 & x<1.5) 
  two = which(x>=1.5)
  
  nmes = dimnames(pheno_to_plot)[[2]][i1:(i1+1)]
  print(nmes)
  title = paste("Correlation at",dimnames(geno)[[2]][j],sep=" ")
  plot(pheno_to_plot[,i1], pheno_to_plot[,i1+1], type='p', col="white",xlab = nmes[1],ylab = nmes[2],main=title,cex=cex)
      lines(pheno_to_plot[zero,i1], pheno_to_plot[zero,i1+1], type='p', col=cols[1],main=title,xlab=nmes[1],ylab=nmes[2],cex=cex)
      lines(pheno_to_plot[one,i1], pheno_to_plot[one,i1+1], type='p', col=cols[2],main = title,xlab=nmes[1],ylab=nmes[2],cex=cex)
      lines(pheno_to_plot[two,i1], pheno_to_plot[two,i1+1], type='p', col=cols[3],main=title,xlab=nmes[1],ylab=nmes[2],cex=cex)
  legend("bottomright", legend =0:2, cex=cex,col=1:3,pch=1:3)
 }
 }
}


