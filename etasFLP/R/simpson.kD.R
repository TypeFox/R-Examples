simpson.kD  <-function(n,k=2){
            a   =simpson.coeff(n)
            if (k==1){
                a
                }
            else{
                as.vector(outer(a,Recall(n,k-1)))
                }
                }
