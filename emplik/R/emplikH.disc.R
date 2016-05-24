########################################################
############ emplikH.disc() #########################
########################################################

emplikH.disc <- function(x, d, y= -Inf, K, fun, 
	                     tola=.Machine$double.eps^.25, theta)
{
n <- length(x) 
if(n <= 2) stop("Need more observations")
if(length(d) != n ) stop("length of x and d must agree") 
if(any((d!=0)&(d!=1))) stop("d must be 0/1's for censor/not-censor")
if(!is.numeric(x)) stop("x must be numeric values --- observed times")

#temp<-summary(survfit(Surv(x,d),se.fit=F,type="fleming",conf.type="none"))
#
newdata <- Wdataclean2(x,d)
temp <- DnR(newdata$value, newdata$dd, newdata$weight,y=y)

otime <- temp$times         # only uncensored time?  Yes. 
orisk <- temp$n.risk
odti <- temp$n.event

###if the last jump is of size 1, we need to drop last jump from computation
last <- length(orisk) 
if (orisk[last] == odti[last]) {
                     otime <- otime[-last] 
                     orisk <- orisk[-last]
                     odti  <- odti[-last]
                     }
######## compute the function g(ti, theta) 
gti <- fun(otime,theta) 

Kcenter <- sum(gti * log(1- odti/orisk))

### the constrain function. To be solved in equation later.

constr <- function(x, Konst, gti, rti, dti, n) { 
               rtiLgti <- rti + x*n*gti
               OneminusdH <- (rtiLgti - dti)/rtiLgti
               if( any(OneminusdH <= 0) ) stop("estimator not well defined")
               sum(gti*log(OneminusdH)) -  Konst } 

##############################################################

differ <- constr(0, Konst=K, gti=gti, rti=orisk, dti=odti, n=n)

if( abs(differ) < tola ) { lam <- 0 } else {
    step <- 0.2/sqrt(n) 
    if(abs(differ) > 200*log(n)*step )   #Why 200 ? 
    { print( Kcenter )
      stop("the given K value is too far away from K_0. \n Move K closer to the above value")
    }

    mini<-0
    maxi<-0   
######### assume the constrain function is increasing in lam (=x) 
    if(differ > 0) { 
    mini <- -step 
    while(constr(mini, Konst=K, gti=gti, rti=orisk, dti=odti, n=n) > 0
 	          && mini > -200*log(n)*step )
    mini <- mini - step 
    } 
    else { 
    maxi <- step 
    while(constr(maxi, Konst=K, gti=gti, rti=orisk, dti=odti, n=n) < 0
                  &&  maxi < 200*log(n)*step )
    maxi <- maxi+step 
    }

    if(constr(mini, Konst=K, gti=gti, rti=orisk, dti=odti, n=n)*constr(maxi, 
                Konst=K, gti=gti, rti=orisk, dti=odti, n=n) > 0 )
    stop("given theta/K is/are too far away from theta0/K0")

# Now we solve the equation to get lambda, to satisfy the constraint of Ho

    temp2 <- uniroot(constr,c(mini,maxi), tol = tola, 
                  Konst=K, gti=gti, rti=orisk, dti=odti, n=n) 
    lam <- temp2$root 
}
####################################################################
rPlgti <- orisk + n*lam*gti

loglik <- 2*sum(odti*log(rPlgti/orisk) +
           (orisk-odti)*log(((orisk-odti)*rPlgti)/(orisk*(rPlgti-odti)) ) )

#?is that right? YES the -2log lik ratio. 
# Notice the output time and jumps has less the last point.
list("-2LLR"=loglik, lambda=lam, times=otime,
                jumps=odti/rPlgti)
}

