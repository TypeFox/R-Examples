############################################################           
WW <- function(x,a,h,method=c("forward","backward")) UseMethod("WW")
############################################################
WW.default <-function(x,a,h,method=c("forward","backward"))
{ method <- match.arg(method)
  n  <- length(x) #Calculating the output matrix (costs matrix)
                                                    CM<- matrix(NA,nrow=n,ncol=n)
   if (method=="forward"){CM<- matrix(NA,nrow=n,ncol=n)
                      f<- matrix(NA,nrow=1,ncol=n+1)
                      hc<-matrix(NA,nrow=n,ncol=n)

                     for(t in 2:n){  CM[1,1]<- a[1]
                                        f[1]<- 0
                                        f[2]<- CM[1,1]
                                        for(k in 1:t){ if(k==t){  CM[t,k]<- f[k]+a[k]
                                                                }
                                                             else { total = 0;
                                                                    for (r in k:(t-1)){hc<-sum(h[k:r])*x[r+1];
                                                                                      total <- total+hc;
                                                                                     CM[t,k]<- f[k]+a[k]+total
                                                                                      }
                                                                         
                                                                   }  
                                                         
                                                      }
                                             f[t+1]<-min(CM[t,],na.rm=TRUE)
                                             
                                        }
                       v=f[-1];
                       TC<- f[n+1];
                             }
  else {V<- matrix(NA,nrow=1,ncol=n)
                           hc<-matrix(NA,nrow=n,ncol=n)
                     
                              for(i in (n-1):1){ CM[n,n]<- a[n]
                                                  V[n]<- CM[n,n]
                                                   V[n+1]<- 0 
                                                 for(j in i:n){ if(j==i){  CM[i,j]<- V[j+1]+a[i]
                                                                         }
                                                                  else { total = 0;
                                                                         for (k in i:(j-1)){hc<-sum(h[i:k])*x[k+1]
                                                                                            total <- total+hc
                                                                                            CM[i,j]<- V[j+1]+a[i]+total
                                                                                            }
                                                                        }                                                           
                                                                }
                                                    V[i]<-min(CM[i,],na.rm=TRUE)
                                                 }
                                v<-V[-(n+1)];
                                TC<-V[1];
         }
                       s <- apply(CM, 1, function(y) which(y == min(y, na.rm = TRUE)));
                       Cuts<-sapply(s, paste, collapse = ' or ');
                       
                       ww<-list(TVC=TC, Jt = Cuts, Solution=CM, call=sys.call()) 
                       class(ww)<-"WW"  
                       ww
}

############################################################
print.WW <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\nTVC:\n")
print(x$TVC)
cat("\nSolution:\n")
print(x$Solution)
cat("\nJt:\n")
print(x$Jt)
}
##########################################################




