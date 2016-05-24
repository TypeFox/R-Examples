multi.cont <- function (DAT, baseline = NULL) {

  
    P.DAT  <- DAT/sum(DAT)
    p.diag <- diag(P.DAT)
    PC     <- sum(diag(P.DAT))
    p.forc <- apply(P.DAT, 1, sum)
    p.obs  <- apply(P.DAT, 2, sum)
    K      <- nrow(DAT)
    S      <- matrix(NA, ncol = nrow(DAT), nrow = nrow(DAT))
    kap    <- 1/(K - 1)
    p.base <- p.obs
 
if(!is.null(baseline) ){
 p.base <- baseline
  if(nrow(DAT)!= ncol(DAT) | nrow(DAT)!= length(baseline)){ stop("Dimension of contingency table does not
correspond to length of baseline probabilities") }
} ## close if null


## uses notation from Joliffe p. 89-90.

b <- 1/(K-1)
a <- numeric()

for(k in 1:K){

      a[k] <- (1-sum(p.base[1:k]))/sum(p.base[1:k])
    }
### fill diagonal
for(i in 1:K){
if(i ==1){
S[i,i]<- b* sum(a[1:(K-1)] )
    }else if(i==K){
  S[i,i]<- b*(sum(1/a[1:(i-1)] ) )} else{
      S[i,i]<- b*(sum(1/a[1:(i-1)] ) +  sum(a[i:(K-1)] )) } ## close else
} ## close i

### fill off diagonal   

 for(i in 1:(K-1)){
    for(j in (i+1):K){
if(i == 1 & j == K){
S[j,i] <-   S[i,j] <- b*(-1*(j-i))
}else
if(i == 1){    
S[j,i] <-   S[i,j] <- b*(-1* (j-i) +  sum(a[j:(K-1)] ))}else
if(j == K){
S[j,i] <-   S[i,j] <- b*( sum(1/a[1:(i-1)] )- (j-i) )
}else
{S[j,i] <-   S[i,j] <- b*( sum(1/a[1:(i-1)] )- (j-i) +  sum(a[j:(K-1)] ))}


} ## close i
} ## close j
  
    GS <- sum(P.DAT * S)
    bias <- p.forc/p.obs
    pc <- diag(P.DAT)/p.obs
    f <- h <- ts <- far <- bias2<- pc2<- numeric()
for(i in 1:nrow(P.DAT)) {
a <-P.DAT[i,i]
b <- P.DAT[i, -i]  
c <- P.DAT[-i, i]
d <- P.DAT[-i, -i]

h[i] <- sum(a)/ (sum(a) + sum(c))  ## hit rate also pod
f[i] <- sum(b)/ (sum(b) + sum(d) ) ## false alarm rate
far[i]<- sum(b)/(sum(a) + sum(b) ) ## false alarm ratio
pc2[i]<- (sum(a)  + sum(d))/ sum(P.DAT)
ts[i] <- sum(a)/(sum(a) + sum(b) + sum(c) )
# bias2[i] <- (sum(a) + sum(b))/ (sum(a) + sum(c))
}

    
   # far <- (p.forc - p.diag)/p.forc
   # d <- numeric()
   # for (i in 1:nrow(P.DAT)) {
   #     d[i] <- sum(P.DAT[-i, -i])
   #   }
   # ts <- p.diag/(1 - d)

    hss <- (sum(p.diag) - sum(p.base * p.forc))/(1 - sum(p.base * 
        p.forc))
    pss <- (sum(p.diag) - sum(p.base * p.forc))/(1 - sum(p.base * 
        p.base)) 


    return(list(pc = PC, bias = bias,  ts = ts, hss = hss, pss = pss, gs = GS,   pc2 = pc2, h = h, f = f,
           false.alarm.ratio = far))

  } ## close function

############### example

#USFMATemp <- read.table("~/Desktop/verify/USFMA.txt",header=FALSE)
#FMA <- table(USFMATemp)
#
#JJA <- matrix(c(3,8,7,8,13,14,4,18,25), nrow = 3 )
#
#A<- multi.cont(FMA)
#B<- multi.cont(FMA, baseline = c(0.3, 0.4, 0.3) )
#
#D<- verify(FMA, baseline = c(0.3, 0.4, 0.3), frcst.type = "cat", obs.type = "cat"  )
#
#summary(D)
#
#
#A$HSS
#B$HSS
#
#multi.cont.alt(JJA)

