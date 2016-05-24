gbtest <-
function(dta,design,design0=NULL,graphthresh=0.05,nbsamples=1000) {
       # computes the significant intervals
       # nbsimul controls the precision in the calculation of the
       # significant sequence length

   lseq = function(x) {  # x is a vector of binary values (0 and 1)
      T = length(x)      # lseq returns the lengths of the sequences of 1
      changepoints = c(x[1],diff(x))
      start = (1:T)[changepoints== 1]
      end = (1:T)[changepoints== -1]-1
      if (length(end)<length(start)) end[length(start)] = T
      list(start=start,length=end-start+1)
   }

   signiflength = function(n,T,rho,graphthresh,nbsamples=1000) {    # returns the length of the shortest run 
      sigma = matrix(rep(0:(T-1),T),nrow=T,ncol=T,byrow=TRUE)       # length that occurred in less than 5 % of 
      sigma = abs(sigma-t(sigma))                                   # nbsimul simulations
      sigma = rho^sigma                                             # covariance of an AR(1) process
      ttest = rmt(nbsamples,S=sigma,df=n-1)                         # simulations of nbsimul ttests under the null hypothesis
      over = abs(ttest)>qt(1-graphthresh/2,n-1)                     # turns into binary values (1 if significant, 0 otherwise)
      lgth = apply(over,1,function(s) lseq(s)$length)               # calculates the lengths of significant sequences 
                                                                 # for each simulation
      maxlgth = unlist(lapply(lgth,function(x) ifelse(length(x)>0,max(x),0))) 
                                                                 # keeps the longest run
      qtile = quantile(maxlgth[maxlgth!=0],probs=0.95)              # returns the 95 % quantile of these maximum run lengths
      check = mean(maxlgth>=qtile)<=0.05                            # check if the returned quantile will ensure
      if (!check) qtile = qtile+1                                   # that the risk is lower than the significance level
      return(qtile)
   }

   acfdp = function(dta) {           # estimates autocorrelation 
      acfs = apply(dta,1,function(dp) 
          acf(dp,lag.max=1,plot=FALSE,type="covariance",demean=FALSE)$acf[,1,1])
                                  # estimates autocovariance for each difference potential
      list(rho=mean(acfs[2,]/acfs[1,]),std=sqrt(mean(acfs[1,])))    
                                  # returns the estimates of rho and sigma
   }

   if (is.null(design0)) design0 = matrix(1,nrow=nrow(dta),ncol=1)
   erpdta = as.matrix(dta)
   design = as.matrix(design)
   design0 = as.matrix(design0)
   if (typeof(erpdta)!="double") stop("ERPs should be of type double")
   if (nrow(erpdta)!=nrow(design)) stop("dta and design should have the same number of rows")
   if (nrow(erpdta)!=nrow(design0)) stop("dta and design0 should have the same number of rows")
   if (ncol(design)<=ncol(design0)) stop("design0 should have fewer columns than design")
   idsignal = NULL
   for (j in 1:ncol(design)) {
      cj = apply(design0,2,function(x,y) all(x==y),y=design[,j])
      if (all(!cj)) idsignal = c(idsignal,j)
   }
   if (length(idsignal)!=(ncol(design)-ncol(design0))) stop("the null model design0 should be nested into the non-null model design")
   if (typeof(graphthresh)!="double") stop("graphthresh should be of type double")
   if ((graphthresh<=0)|(graphthresh>=1)) stop("graphthresh should be in ]0,1[, typically 0.05")

   n = nrow(erpdta)                  # sets the value for n
   T = ncol(erpdta)                  # sets the value for T 

   pdesign = solve(t(design)%*%design)%*%t(design)
   Proj = design%*%pdesign
   Proj0 = design0%*%solve(t(design0)%*%design0)%*%t(design0)
   rdf1 = nrow(design)-ncol(design)
   rdf0 = nrow(design0)-ncol(design0)

   beta = (pdesign%*%erpdta)[idsignal,]
   if (length(idsignal)==1) beta = matrix(beta,nrow=1)

   lsamp = lapply(1:nbsamples,function(i,p) sample(1:p),p=n)
   lres1 = lapply(lsamp,function(samp,x,Proj) Proj%*%x[samp,],x=erpdta,Proj=diag(n)-Proj) 
   lres0 = lapply(lsamp,function(samp,x,Proj) Proj%*%x[samp,],x=erpdta,Proj=diag(n)-Proj0) 
   lscer1 = lapply(lres1,function(res) as.vector(t(rep(1,nrow(res)))%*%res^2))
   lscer1 = matrix(unlist(lscer1),nrow=nbsamples,byrow=TRUE)
   lscer0 = lapply(lres0,function(res) as.vector(t(rep(1,nrow(res)))%*%res^2))
   lscer0 = matrix(unlist(lscer0),nrow=nbsamples,byrow=TRUE)
   lfstat = ((lscer0-lscer1)/(rdf0-rdf1))/(lscer1/rdf1)
   accoef = acfdp(lfstat)              # computes autocorrelation
   maxl = signiflength(n=rdf1,T=T,rho=accoef$rho,graphthresh=graphthresh,nbsamples=nbsamples)
       # computes the significant sequence length

   res1 = erpdta-Proj%*%erpdta
   scer1 = as.vector(t(rep(1,n))%*%res1^2)
   res0 = erpdta-Proj0%*%erpdta
   scer0 = as.vector(t(rep(1,n))%*%res0^2)
   fstat = ((scer0-scer1)/(rdf0-rdf1))/(scer1/rdf1)
   pval = pf(fstat,df1=rdf0-rdf1,df2=rdf1,lower.tail=FALSE)

   over = pval<=graphthresh       # over=1 if the t-test is significant
                                  # 0 otherwise
   lgth = lseq(over)              # computes the lengths of sequences
   if (length(lgth$length)==0) {
      intervals = integer(0)                      # computes the significant intervals
      signifseq = logical(0)
   }
   if (length(lgth$length)>0) {
      intervals = lapply(1:length(lgth$length),function(i,l) l$start[i]:(l$start[i]+l$length[i]-1),l=lgth) # computes the significant intervals
      signifseq = lgth$length > maxl    # says which sequences are significant
   }
   list(nbsignifintervals=sum(signifseq),intervals=intervals[signifseq],significant=unlist(intervals[signifseq]),signal=beta,rho=accoef$rho,r2=(1-1/(1+fstat*((rdf0-rdf1)/rdf1))))
}
