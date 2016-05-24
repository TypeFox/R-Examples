nptl <-
function(n,p=0.95,gam=0.95){
#  function nptl(n,p,gam)  
#     For a random sample of size n calculate largest value
#      of m such that with  confidence level gamma
#      100p percent of population lies  below the
#      mth largest data value in the sample... see Section 3.6
# USAGE: nptl(n,p,gam)
# ARGUMENTS: n: sample size p: defined above
#            gam:  confidence level for one-sided interva
# VALUE: m 
# DETAILS: Requries R function qbeata(p,par1,par2)
# REFERENCES:
#   Sommerville, P.N. (1958) Annals Math Stat pp 599-601
k <- ceiling(n*p)
pv <- qbeta(1-gam,k,n+1-k)
while( pv < p && k < n+1){
k <- k + 1
if( k == n + 1) next
pv<-qbeta(1-gam,k,n+1-k)
}
if( k <= n) m<- n+1-k else m<- NA
m
}

