get.bj <-
function(c, B, p,j){
         
       Bj=NULL
           for(k in 1:c-1){
            if(c>2 & j<=p*(c-1)){

             Bj=cbind(Bj,B[,j])
             j=j+p
             
            }
            if(c==2){
               Bj=B[,j] 
            }

          }
       
     return(Bj)
}
