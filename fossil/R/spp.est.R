'spp.est' <-
function(x, rand = 10, abund = TRUE, counter = FALSE, max.est = 'all')
{
  #to mark off it it is a single vector or a matrix
  if (ncol(as.matrix(x))==1) vec<-TRUE
  else vec<-FALSE
  #first we separate out if it is abundance or occurrence based that we're running
  if (abund==TRUE) {
    if (max(x)==1) warning("cannot use incidence data for abundance-based analyses. If the data is incidence based, please run this function again with the option of 'abund=FALSE'")
    if (vec==TRUE) {
      n<-length(x)
      m<-sum(x)
    }
    else {
      n<-ncol(x)
      m<-n
    }
    #we use the abundance based versions of Chao's species estimators if the matrix is abundance based
    ests<-list(chao1,chao.sd,ACE,jack1)
    cn<-c("N.obs", "S.obs", "S.obs(+95%)", "S.obs(-95%)", "Chao1", "Chao1(upper)", "Chao1(lower)", "ACE", "ACE(upper)", "ACE(lower)", "Jack1", "Jack1(upper)", "Jack1(lower)")
  }
  else {
    if (vec==TRUE) warning('You are using a single sample for an occurrence based analysis, which will give invalid results. Please check your data and ensure it is a matrix with multiple samples and multiple species')
    x[x>1]<-1
    n<-ncol(x)
    m<-n
    #these are the incidence/occurrence based estimators we use
    ests<-list(chao2,chao.sd,ICE,jack1)
    cn<- c("# Samples", "S.obs", "S.obs(+95%)", "S.obs(-95%)", "Chao2", "Chao2(upper)", "Chao2(lower)", "ICE", "ICE(lower)", "ICE(lower)", "Jack1", "Jack1(lupper)", "Jack1(lower)")
  }
  if (max.est != 'all' && max.est < m) m <- max.est
  ##se is final output table with all the values returned, with a dimension of 13 (the various estimators) by m, the number of samples we can resample (ie typically the number of localities)
  se<-matrix(,m,13)
  se[,1] <- 1:m ##number of samples taken
  le<-length(ests)
  ##if x is a abundance vector, this creates a single vector to sample of length==sum(x)
  if (abund==TRUE & vec==TRUE) ss<-rep(1:n,x)
  else ss<-1:n
  for (i in 1:m) {
    avg<-matrix(,rand,le+1)
    for (j in 1:rand) {
      ssc<-sample(ss,i,FALSE)
      if (abund==TRUE & vec==TRUE) {
        b<-numeric(n)
        for(k in 1:i) b[ssc[k]]<-b[ssc[k]]+1
      }
      else {
        if (i>1) b<-rowSums(x[,ssc])
        else b<-x[,ssc]
      }
##run all stats at this point
## avg is the matrix which contains all the randomizations for each loop (ie if there are 100 samples and 10 randomizations, there will be 100 avg tables made each with 10 rows)
      avg[j,1]<-length(b[b>0])
      for (k in 1:le) avg[j,(k+1)]<-ests[[k]](b)
    }
    se[i,2]<-mean(avg[,1]) ##sobs
    se[i,3]<-se[i,2]+1.96*(sd(avg[,1])) ##sobs + stdev
    se[i,4]<-2*se[i,2]-se[i,3] ##sobs - stdev
    se[i,5]<-mean(avg[,2]) ##mean chao
    se[i,6]<-se[i,5]+mean(avg[,3],na.rm=TRUE) ##chao-mean chaosd
    se[i,7]<-se[i,5]-mean(avg[,3],na.rm=TRUE) ##chao-mean chaosd
    ##next line makes sure if the chao.sd is a number, otherwise it is recalculated as just a 95% confidence interval 
    ##chao.sd is a little finicky
    if (is.nan(se[i,6])==TRUE) se[i,6]<-se[i,5]+1.96*sd(avg[,2],na.rm=TRUE)/sqrt(rand)
    if (is.nan(se[i,7])==TRUE) se[i,7]<-se[i,5]-1.96*sd(avg[,2],na.rm=TRUE)/sqrt(rand)
    ## check to see if any numbers are real
    re<-avg[,4][avg[,4]<Inf & is.na(avg[,4])==FALSE]
    se[i,8]<-mean(avg[,4][avg[,4]<Inf],na.rm=TRUE)##mean ACE
    if (length(re)>0) se[i,9] <- se[i,8]+1.96*(sd(avg[,4][avg[,4]<Inf],na.rm=TRUE))  ##ACE+acesd
    se[i,10]<-2*se[i,8]-se[i,9]   ##ACE-acesd
    ## check to see if any numbers are real
    re<-avg[,5][avg[,5]<Inf & is.na(avg[,5])==FALSE]
    se[i,11]<-mean(avg[,5][avg[,5]<Inf],na.rm=TRUE)##mean jack
    if (length(re)>0) se[i,12]<-se[i,11]+1.96*(sd(avg[,5][avg[,5]<Inf],na.rm=TRUE))  ##jack+jacksd
    se[i,13]<-2*se[i,11]-se[i,12]   ##jack-jacksd
    ##etc  more to be added
    if (counter==TRUE && i%%100==0) print(paste(i, "out of", m, "runs"))
  }
  if (abund==TRUE) attr(se,"data.type")<-"abundance"
  else attr(se,"data.type")<-"presence\absence"
  colnames(se)<-cn
  return(se)
}


