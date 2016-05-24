#' VN SD Efficient Set
#'
#' Determines the set of stochastic dominant variables for various degrees.
#' @param A data.frame of variables.
#' @param degree Degree of stochastic dominance test
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126. \url{http://www.scirp.org/Journal/PaperInformation.aspx?PaperID=63817}.
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100); z<-rnorm(100)
#' A<-data.frame(x,y,z)
#' \dontrun{VN.SD.Efficient.Set(A,1)}
#' @export



VN.SD.Efficient.Set <- function(A,degree) {
  n <- ncol(A)
  max_target <- max(A)
  LPM_order<- numeric(0)
  Dominated_set<- numeric(0)
  current_base<- numeric(0)


  for (i in 1:n){
    LPM_order[i] <- LPM(1,max(A),A[,i])

    }

  final_ranked <- A[,order(LPM_order)]

  current_base<-1

  if(degree==1){
  for (i in 1:(n-1)) {

      base<- final_ranked[,current_base[length(current_base)]]

      challenger <- final_ranked[,i+1]

      if (FSD(base,challenger)==1){ current_base[i]<- current_base[length(current_base)]}


      if (FSD(base,challenger)==0){

        for (j in current_base){
          base<- final_ranked[,j]
          if (FSD(base,challenger)==0){ next }
          else
            {Dominated_set[i] <- i+1  }
          }

       current_base[i]<- i+1}

      else {Dominated_set[i]<- i+1 }
  }
}


 if(degree==2){
   for (i in 1:(n-1)) {

     base<- final_ranked[,current_base[length(current_base)]]

     challenger <- final_ranked[,i+1]

     if (VN.SSD.uni(base,challenger)==1){ current_base[i]<- current_base[length(current_base)]}


     if (VN.SSD.uni(base,challenger)==0){

       for (j in current_base){
         base<- final_ranked[,j]
         if (VN.SSD.uni(base,challenger)==0){ next }
         else
         {Dominated_set[i] <- i+1  }
       }

       current_base[i]<- i+1}

     else {Dominated_set[i]<- i+1 }
   }
 }


 if(degree==3){
   for (i in 1:(n-1)) {

     base<- final_ranked[,current_base[length(current_base)]]

     challenger <- final_ranked[,i+1]

     if (TSD(base,challenger)==1){ current_base[i]<- current_base[length(current_base)]}


     if (TSD(base,challenger)==0){

       for (j in current_base){
         base<- final_ranked[,j]
         if (TSD(base,challenger)==0){ next }
         else
         {Dominated_set[i] <- i+1  }
       }

       current_base[i]<- i+1}

     else {Dominated_set[i]<- i+1 }
   }
 }


  if(length(Dominated_set)>0){
    SD_A  = (data.frame((final_ranked[-na.omit(Dominated_set)])))
    return(SD_A)
    }

    else {print(colnames(final_ranked))}


}
