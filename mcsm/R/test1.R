test1=function(Nsim=10^4,df=6){
# simple derivation of a chi-square sample

midf=df/2
U=matrix(data=runif(midf*Nsim),nrow=midf)
2* apply(-log(U),2,sum)
}

test2=function(Nsim=10^4,df=6){
# simple derivation of a chi-square sample

rchisq(Nsim,df)
}
