nPPV <-
function(propP, se, sp, prev, PPV0, conf.level=0.95, power=0.8)
{
propQ <- (1-propP)
# Eqn 2.9
sigmasq1 <- (1-se)/(se*propP)+(sp)/((1-sp)*propQ)
# Eqn 2.8
sigmasq2 <- (se)/((1-se)*propP)+(1-sp)/((sp)*propQ)
# Eqn 2.5
phi1 <- log(1-sp)-log(se)
# Eqn. 2.6
phi2 <- log(1-se)-log(sp)

TRUEppv<-ppv(p=prev, se=se, sp=sp)
PPVlimit <- 1/(1+exp(phi1)*(prev/(1-prev)))
z1alpha <- qnorm(p=conf.level)
z1beta <- qnorm(p=power)

if(TRUEppv<=PPV0){n<-rep(NA)}else{
  nest <- ((z1alpha+z1beta)^2 * sigmasq1) / (phi1-log((prev/(1-prev))*((1-PPV0)/PPV0)))^2
  n<-ceiling(nest)
 }

if(!any(is.na(n))){
n1se <- pmin(n*propP*se, n*propP*(1-se))
n0sp <- pmin(n*(1-propP)*sp, n*(1-propP)*(1-sp))

if(any(n1se<5)){
wn1 <- which(n1se<5)
if(length(wn1) <=4 ){wn1show <- paste(n[wn1], collapse=", ")}else{wn1show <- paste(n[wn1[1]],", ", n[wn1[2]],", ..., ", n[wn1[length(wn1)]], sep="")}
warning(paste("Some sample sizes n (", wn1show, ") might be too small to expect validity of asymptotic formulas for anticipated sensitivity (se), and the corresponding proportion of positives (propP).", sep=""))
}

if(any(n0sp<5)){
wn0 <- which(n0sp<5)
if(length(wn0) <=4 ){wn0show <- paste(n[wn0], collapse=", ")}else{wn0show <- paste(n[wn0[1]],", ", n[wn0[2]], ", ..., ", n[wn0[length(wn0)]], sep="")}
warning(paste("Some sample sizes n (", wn0show, ") might be too small to expect validity of asymptotic formulas for anticipated specificity (sp), and the corresponding proportion of negatives (1-propP).", sep=""))
}}

out<-list(n=n, se=se, sp=sp, prev=prev, PPV0=PPV0, TRUEPPV=TRUEppv, propP=propP, power=power, conf.level=conf.level)
return(out)
}