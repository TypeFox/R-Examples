bayesNI <-
function(x1,x2,n1,n2, dm='OR',rho, m=10,noninform.prior=TRUE,w1,w2,TWE=1,zeta=0.5, plot.prior=FALSE){
  
if(rho>1|rho<0)stop("0<=rho<=1 condition violated")
if(rho==1 & dm=='RD') stop("rho=1 for risk difference, p[H0]=0")
if(rho==0 & (dm=='RR'|dm=='OR')) stop("rho=0 for relarive risk/odds ratio, p[H0]=0")


if(dm=='RD'){gfun<-function(x) x-rho
             cat(c("H0: theta2<=theta1-",rho,"vs. H1: theta2>theta1-",rho),fill=TRUE)
            }
  else{if(dm=='RR'){gfun<-function(x)rho*x
                    cat(c("H0: theta2/theta1<=",rho,"vs. H1: theta2/theta1>",rho),fill=TRUE)
                    }
        else{gfun<-function(x)rho*x/(1-x+rho*x)}
             cat(c("H0: odds(theta2)/odds(theta1)<=",rho," vs. H1: odds(theta2)/odds(theta1)>",rho),fill=TRUE)
            }

cat(c("weight assignment in TWE: ",1-zeta," Type I Error | ",zeta," Type II error"),fill=TRUE)

#Construct Noninformative Prior
  if(noninform.prior==TRUE){
  temp.prior<-noninform_prior(m=m,gfun=gfun)
  w2<-w1<-temp.prior$w
  pH01=temp.prior$pH0
      }
   else{
       if(length(w1)!=m|length(w2)!=m|sum(w1)!=1|sum(w2)!=1)stop("Please check input values for w1 and w2")
       pH01=prior_prob(m=m,gfun=gfun)
       }


#plot prior if plot.prior=TRUE
#if(plot.prior=TRUE) plotprior(m=m,w1=w1,w2=w2)

#Calculate observed Log(Bayes Factor)
logBF<-calcBF(x1,x2,n1,n2,m=m,w1=w1,w2=w2,pH0=pH01,gfun=gfun)


#Calculate the cut-off value (the critical value for the rejection region) LogBF0
e1w=1-zeta
if(TWE==1){L0<-log((1-zeta)/zeta)}
else{L0<-X1X2BFtable(n1,n2,m=m,w1=w1,w2=w2,pH0=pH01,e1w=e1w,gfun=gfun)$L0_BE}


cat(c("logBF(x1,x2)=",logBF,"L0=",L0),fill=TRUE)

list(logBF=logBF, L0=L0, w1=w1,w2=w2)
}

