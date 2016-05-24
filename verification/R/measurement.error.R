measurement.error <- function( obs, frcs = NULL, theta = 0.5, CI = FALSE, t = 1, u = 0, h = NULL, ...){
### if frcs is null, its assumed that obs is a vector with
### assume data entered as c(n11, n10, n01, n00)
if(is.null(frcs) & length(obs) ==4 ){
  print(" Assume data entered as c(n11, n01, n10, n00)")
n11<- obs[1]
n10<- obs[3]
n01<- obs[2]
n00<- obs[4]} else{

### check to see if frcs is [0,1] if not convert
if( prod(unique(obs) %in% c(0,1) ) == 1 ){ ## if obs is not binomial

  if(is.null(h)){ frcs <- as.numeric(frcs > theta)  } else {frcs <- as.numeric(frcs > h) }
}# close if not unique

A<- table(data.frame(obs, frcs) )

n11 <- A[2,2]
n00 <- A[1,1]
n10 <- A[2,1]
n01 <- A[1,2]
}# close is.null else


    # No error checking, but either n_ij can be a vector or theta can be;
    # both can NOT be vectors.
    # briggs wib2004@med.cornell.edu
    
    n<-n11+n01+n10+n00;
    # to determine which of RARE of COMM(ON) is used for each set
    z<-((n10+n11)/n);
    z<-(z <= theta)*1;

        q11<-n11/(n11+n01);
        p11<-(q11-u)/(t-u);
        q00<-n00/(n10+n00);
        p00<-(q00-(1-t))/(t-u);
        px<-(n11+n01)/n
        py<-(n11+n10-n*u)/(n*(t-u))
            #k_RARE <- (n11*(1-u-theta*(t-u))-n01*(u+theta*(t-u)))/((n11+n10-n*u)*(1-theta));
    k_RARE <- (p11-theta)*px/((1-theta)*py);
    G_RARE <- 2*n11*log(p11/theta) + 2*n01*log((1-p11)/(1-theta));
            #k_COMM <- (n00*(t-(1-theta)*(t-u))-n10*(1-t+(1-theta)*(t-u)))/((n00+n01-n*(1-t))*theta);
    k_COMM <- (p00-1+theta)*(1-px)/(theta*(1-py));
    G_COMM <- 2*n00*log(p00/(1-theta))  +2*n10*log((1-p00)/(theta));
    
    k<-k_RARE*z+k_COMM*(1-z);
    G<-G_RARE*z+G_COMM*(1-z);

    # if k<0 then G=0 and p=.5
    z<-(k>0)*1;
    G<-G*z;
    p<-pchisq(G,1,lower.tail=FALSE)/2;
    p<-p*z+0.5*(1-z);
    bigP<-p11*z+p00*(1-z)

    if(CI){
        ciLO<-0;ciHI<-0;
        rootLO<-0;rootHI<-0;
        if(length(theta>=1)){
            f_RARE<-function(p11,n11,n01,theta)
                2*n11*log(p11/theta) + 2*n01*log((1-p11)/(1-theta))-
                  (2*n11*log((n11/(n11+n01))/theta) +
                   2*n01*log((n01/(n11+n01))/(1-theta))-4.25);
                # 3.84 for climate
                # 4.25 for markov
            f_COMM<-function(p00,n00,n10,theta)
                2*n00*log(p00/(1-theta)) + 2*n10*log((1-p00)/theta)-
                  (2*n00*log((n00/(n00+n10))/(1-theta)) +
                   2*n10*log((n10/(n00+n10))/theta)-4.25);
            for (i in 1:length(theta)){
                    zz<-((n10[i]+n11[i])/n[i]);
                    zz<-(zz<=theta[i]);            
                    tol<-0.0001
                    if(zz){
                        hi<-p11[i]-tol
                        lo<-p11[i]+tol;# print(c(0,hi,lo,n11[i],n01[i],theta[i]))
                        rootLO[i]<-uniroot(f_RARE,lower=0.001,upper=hi,n11=(n11[i]+tol),
                                           n01=(n01[i]+tol),theta=theta[i])$root
                        ciLO[i] <- (rootLO[i]-theta[i])*px[i]/((1-theta[i])*py[i]);
                        if(n01[i]==0) {
                            rootHI[i]<-1
                        } else {
                            rootHI[i]<-uniroot(f_RARE,lower=lo,upper=(1-tol),n11=(n11[i]+tol),
                                               n01=(n01[i]+tol),theta=theta[i])$root
                        }
                        ciHI[i] <- (rootHI[i]-theta[i])*px[i]/((1-theta[i])*py[i]);
                    } else {
                        hi<-p00[i]-tol
                        lo<-p00[i]+tol; # print(c(1,hi,lo,n00[i],n10[i],theta[i]))
                        rootLO[i]<-uniroot(f_COMM,lower=0.001,upper=hi,n00=(n00[i]+tol),
                                           n10=(n10[i]+tol),theta=theta[i])$root    
                        ciLO[i] <- (rootLO[i]-1+theta[i])*(1-px[i])/(theta[i]*(1-py[i]));                        
                        rootHI[i]<-uniroot(f_COMM,lower=lo,upper=(1-tol),n00=(n00[i]+tol),
                                           n10=(n10[i]+tol),theta=theta[i])$root    
                        ciHI[i] <- (rootHI[i]-1+theta[i])*(1-px[i])/(theta[i]*(1-py[i]));                        
                    }
            }
        }
    } else {
        ciLO <- CI
        ciHI <- CI
    }
    answer<-list(z=z,n=n,k=k,G=G,p=p,theta=theta,ciLO=ciLO,ciHI=ciHI);
return(answer)
  }

#skill <- function(x,y,theta=0.5,CI=FALSE,t=1,u=0) {
#    # briggs wib2004@med.cornell.edu   
#    n<-length(theta);
#    n11<-0;n01<-0;n10<-0;n00<-0;
#    for (i in 1:n){
#        fit<-table((x>theta[i])*1,y)
#        n11[i]<-fit[2,2];
#        n01[i]<-fit[2,1];
#        n10[i]<-fit[1,2];
#        n00[i]<-fit[1,1];
#    }
#    answer<-skillTABLE(n11,n01,n10,n00,theta,CI,t=1,u=0);
#return(answer)
#  }



