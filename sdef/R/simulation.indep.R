simulation.indep <-
function(n,GammaA=2,GammaB=2,epsilonM=0,epsilonSD=1,r1,r2,DEfirst,DEsecond){

FC1=c()
FC2=c()
delta=rgamma(n,shape=GammaA,scale=GammaB)
epsilon1=rnorm(n,epsilonM,epsilonSD)
epsilon2=rnorm(n,epsilonM,epsilonSD)
names=c()
#Group 1 : DE in the first experiment
for(i in (1):(DEfirst)){
    x=rbinom(1,1,0.5)
    if(x==1){FC1[i] <- delta[i] + epsilon1[i]*r1}
    if(x==0){FC1[i] <- -delta[i] - epsilon1[i]*r1}
    FC2[i] <- epsilon2[i]*r2
    names[i] <- "DEfirst"
}

#Group 2 : DE in the second experiment
for(i in (DEfirst+1):(DEfirst+DEsecond)){
    x=rbinom(1,1,0.5)
    FC1[i] <- epsilon1[i]*r1
    if(x==1){FC2[i] <- delta[i] + epsilon2[i]*r2}
    if(x==0){FC2[i] <- -delta[i] - epsilon2[i]*r2}
    names[i] <- "DEsecond"

}

#Group 3 : Not DE in Both experiments
for(i in (DEfirst+DEsecond+1):(n)){
    FC1[i] <- epsilon1[i]*r1
    FC2[i] <- epsilon2[i]*r2
    names[i] <- "Null"
}

##############################################
#Assign the Pvalues

Pval1 = c()
Pval2 = c()

for(i in 1:n){
Pval1[i] <- 2*pnorm(-abs(FC1[i]/r1))
Pval2[i] <- 2*pnorm(-abs(FC2[i]/r2))
}

##############################################
return(list(names=names,FC1=FC1,FC2=FC2,Pval=cbind(Pval1,Pval2)))
}

