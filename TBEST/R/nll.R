nll<-function(parmhat,z){
                a<-parmhat[1]
                k<-parmhat[2]
                n<-length(z)
                m<-max(z)
                if(k>1)L<- -Inf
                if(k<=1&a<max(0,k*m))L<- -Inf
                if(k<=1&a>=max(0,k*m)){
                        if(k<.Machine$double.eps)L<- -n*log(a)-1/a*sum(z)
                        if(k>=.Machine$double.eps)L<- -n*log(a)+(1/k-1)*sum(log(1-k*z/a))
                }
                return(-L)
        }
