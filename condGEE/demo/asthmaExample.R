# data is available from: http://www.blackwellpublishing.com/rss/Volumes/Cv52p3.htm
# it was previously used by Duchateau et al. JRSSC volume 52 (2003), part 3, pages 355-363

n.subjects <- length(unique(asthma$id.w))

med.nar <- 4       #median time not at risk before current gap     
med.rec <- 40      #median length of at risk period

i <- 1 ; k <- 1
asth <- NULL
while(k <= n.subjects) 
{
  n <- asthma$nn[i]         #number of gaps for subject k
  inds <- i:(i+n-1)         #indices of subject k's gaps 
  gaps <- asthma$stop.w[inds] - asthma$start.w[inds]

  if(n > 1)
  {
    nar <- (asthma$start.w[inds[-1]]-asthma$stop.w[inds[-n]]) > med.nar
    rec <- gaps[1:(n-1)] > med.rec
  }
  age1 <- asthma$start.w[inds]>182 & asthma$start.w[inds]<366
  age2 <- asthma$start.w[inds]>=366

  subj.k <- cbind(asthma$id.w[i], log(gaps), asthma$st.w[i:(i+n-1)], asthma$trt.w[i],
                !asthma$fevent[i:(i+n-1)], c(0,nar), c(0,rec), asthma$trt.w[i]*c(0,nar),
                 age1, age2)

  asth <- rbind(asth, subj.k)

  k <- k + 1
  i <- i + n
}  

#the results look at little different than Clement and Strawderman (2009) because our code
#has been updated. Qualitatively nothing has changed though.

start <-  c(4.4, 0.4, -1, 0.3, 0.8, -0.5, -0.3, -0.5, 3)
condGEE(asth,start,k1=K1.t3,k2=K2.t3)