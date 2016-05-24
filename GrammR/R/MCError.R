MCError <-
function(True,Est){
        n <- length(True);
        m <- length(Est);
        A <- 0;
        if (m == n){
                    for (i in 1:(n-1)){
                            for (j in (i+1):n){
                            A <- A + ((True[i] == True[j])*(Est[i] != Est[j]) + (True[i] != True[j])*(Est[i] == Est[j]));
                    }
          }
          A <- A/choose(n,2);
          return(A);
         }
          else{ return(NA);}
}
