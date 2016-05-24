Fst = function(gen,fam){ 
  anyNA = function(x) any(is.na(x))
  # Weir and Cockerham (1984): Unbiased Fst from weighted AOV 
  if(anyNA(fam)) stop("No missing fam is allowed")
  if(anyNA(gen)) gen[is.na(gen)]=1 # replacing missing
  # number of markers and populations
  m = ncol(gen)
  r = length(unique(fam))
  n = table(fam)
  nbar = mean(n)
  # function to calculate binomial variance
  Var = function(snp){p=mean(snp,na.rm=T)/2;return(p*(1-p))}
  # function to estimate p and heterozigousity by family
  p_i = function(snp,fam) tapply(snp,fam,mean)/2 
  h_i= function(snp,fam) tapply(snp,fam,function(x)mean(x==1,na.rm=T))
  # Whole-Population parameters
  C2 = var(n)/(mean(n)^2) # squared CV of sample sizes
  nc = nbar*(1-C2/r)
  # Parameters by population
  M = apply(gen,2,p_i,fam=fam) # observed P by family by allele
  N = matrix(n,r,m)
  Pbar = colSums(M*N)/(r*nbar) # expected P
  eP = t(matrix(Pbar,m,r))
  s2 = colSums(N*(M-eP)^2)/((r-1)*nbar) # var over pop
  H = apply(gen,2,h_i,fam=fam) # observed heterozigousity
  Hbar = colSums(H*N)/(r*nbar)
  V = Pbar*(1-Pbar)
  # Three variance components: A, B and C
  Va = nbar/nc*(s2-(1/(nbar-1))*(V-(r-1)/r*s2-Hbar/4))
  Vb = (nbar/(nbar-1))*(V-(r-1)/r*s2-(2*nbar-1)/(4*nbar)*Hbar)
  Vc = Hbar/2
  # F - statistics!
  fst = Va/(Va+Vb+Vc); fst[is.nan(fst)|fst<=0]=0
  fit = 1-Vc/(Va+Vb+Vc); fit[is.nan(fit)|fit<=0]=0
  fis = 1-Vc/(Vb+Vc); fis[is.nan(fis)|fis<=0]=0
  Fstat = data.frame(fst,fit,fis)
  class(Fstat) <- "fst"
  return(Fstat)
}

plot.fst = function(x,..., p=NULL,chr=NULL){
  anyNA = function(x) any(is.na(x))
  class(x)=NULL; FST=x$fst
  C = 1 - 0.1 * (4 + 6 * abs(FST)/max(abs(FST)))
  
  if(is.null(chr)!=T){
    plot(FST,...,xaxt = "n",main="Fixation Index",pch=20,xlab="Genome",col=(rgb(C,C,C)))
  }else{
    plot(FST,...,main="Fixation Index",pch=20,xlab="Genome",col=(rgb(C,C,C)))
  }
  K = smooth(smooth(FST))
  lines(K,type="l",lwd=1,col=rgb(0,0,1,0.6))
  if(is.null(chr)!=T){
    # Adding Chromosome in X axis
    medians=rep(NA,length(chr))
    CHR=c(); for(i in 1:length(chr)) CHR=c(CHR,rep(i,chr[i]))
    for(i in 1:length(chr)) medians[i] = median(which(CHR==i))
    axis(1, at=round(medians), labels=1:length(medians))
    abline(v=cumsum(chr[-length(chr)])+0.5,lty=3)
  } 

  
  # P-value  
  if(is.null(p)!=TRUE){
    par(ask=TRUE); CDF = ecdf(FST);
    Pval = -log(1-CDF(FST)) ; THR = -log(p/length(FST))
    plot(Pval,ylab="Empirical -log probatility",
         xlab="Genome",main="P-values",pch=20,col=1+1*(Pval>THR))
    if(is.null(chr)!=T) abline(v=cumsum(chr[-length(chr)])+0.5,lty=3)
    abline(h=THR,col=rgb(1,0,0,0.4),lty=2);par(ask=FALSE)
    KK = ksmooth(1:length(FST),Pval,"normal",bandwidth = 10)
    lines(KK,type="h",lwd=2,col=rgb(1,0.5,0,0.6))}
  }