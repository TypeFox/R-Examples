simpson.coeff   <-function(n){
                w=array(1,n)
                w[seq(2,n-1,by =2)]=4
                w[seq(3,n-2,by =2)]=2
                return(as.double(w/3.))
                }
