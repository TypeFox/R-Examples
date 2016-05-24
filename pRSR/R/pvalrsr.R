pvalrsr<-function(x,t=1:length(x), nf=150, Numpq=11){
  n<-length(x) 
  if(!(n%in%(6:200)))stop("Series size must be between 6 and 200")
  IND<-which((6:200)==n)
  qhat<-TPout$qhatM[,IND]
  probs<-TPout$probs
  pn<-341 # length(probs)
  obs.test<-GetFitHReg(x, t=t, nf=nf)[1] # observed test statistic
  diff<-abs(obs.test-qhat)  # finding Numpq observations closest to observed test statistic
  min.ind<-order(diff)[1]
  ind.after<-pn-min.ind
  NumDist<-(Numpq-1)/2
  # we need to find the Numpq numebr of indices for interpolation
  if(min.ind>=NumDist && ind.after>=NumDist){
    ind<- (min.ind-NumDist):(min.ind+NumDist)} else {
         if(min.ind>=NumDist && (ind.after)<NumDist){
            ind<- (min.ind-(Numpq-(ind.after+1))):pn} else
              ind<-1:Numpq}
  fit<-lm(qnorm(probs)~qhat+I(qhat^2)+I(qhat^3), subset=ind)
  pval<-1-pnorm(predict(fit, data.frame(qhat=obs.test)))
  pval
}

