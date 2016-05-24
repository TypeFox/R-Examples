
#' Implements Minimum Cost Algorithm to solve transportation problem
#'
#' This function implements Minimum Cost Algorithm to resolve transportation problem and get optimized allocation matrix
#'
#' @param ex_matrix A cost matrix where last column must be the supply and last row must be the demand.
#' Input matrix should not have any missing values (NA), otherwise function will throw an error.
#' @return A List which contrains the allocation matrix and the total optimized cost.
#'
#' @details
#' This function takes a cost matrix (with Supply and Demand) and using North-West Corner approach gives
#' the allocation matrix as well as the calcualted optimized cost.
#' This function checks for degenerated problem but it can't resolve it. User need to resolve by seeing final allocation matrix.
#' If Supply and Demand are not equal Balance Supply/Demand will be stored in Dummy variable.
#'
#' @examples
#' \dontrun{
#'
#' #Input matrix where last row is the Demand and last column is the Supply
#' ex_matrix=data.frame(M1=c(13,10,25,17,210),M2=c(25,19,10,24,240),
#'                      M3=c(8,18,15,18,110),M4=c(13,5,14,13,80),M5=c(20,12,18,19,170),
#'                      Supply=c(430,150,100,130,810),
#'                      row.names = c("W1","W2","W3","W4","Demand"))
#'
#' ex_matrix
#'          M1  M2  M3 M4  M5 Supply
#' W1      13  25   8 13  20    430
#' W2      10  19  18  5  12    150
#' W3      25  10  15 14  18    100
#' W4      17  24  18 13  19    130
#' Demand 210 240 110 80 170    810
#'
#' mincost(ex_matrix)
#'
#' $Alloc_Matrix
#'      M1  M2  M3 M4  M5
#' W1 140 140 110  0  40
#' W2  70   0   0 80   0
#' W3   0 100   0  0   0
#' W4   0   0   0  0 130
#'
#' $Total_Cost
#' [1] 11570
#'
#'}
#'
#' @export
#'
mincost=function(ex_matrix){

        if(sum(is.na(ex_matrix))>0)
          stop("Your matrix has NA values")

        Demand=as.vector(ex_matrix[nrow(ex_matrix),-ncol(ex_matrix)])
        Supply=as.vector(ex_matrix[-nrow(ex_matrix),ncol(ex_matrix)])
        High_Values=max(ex_matrix) + 999999999
        Alloc_Matrix=ex_matrix[-nrow(ex_matrix),-ncol(ex_matrix)]
        ex_matrix=Alloc_Matrix
        Alloc_Matrix[,]=0
        Total_Cost=0
        Total_alloc=0

        while(sum(Supply) != 0 & sum(Demand) != 0)
        {
            tc=which.min(apply(ex_matrix,MARGIN=2,min))  #column of minimum value
            tr=which.min(apply(ex_matrix,MARGIN=1,min))  #row of minimum value
  
            min_curr=min(Demand[tc],Supply[tr])
  
            Demand[tc]=Demand[tc] - min_curr
            Supply[tr]=Supply[tr] - min_curr
            Alloc_Matrix[tr,tc]=min_curr
            Total_Cost=Total_Cost+(min_curr*ex_matrix[tr,tc])
  
            if(Demand[tc]==0)
            {
              ex_matrix[,tc]=rep(High_Values,nrow(ex_matrix))
            }else if(Demand[tc]==Supply[tr])
            {
              ex_matrix[tr,]=rep(High_Values,ncol(ex_matrix))
              ex_matrix[,tc]=rep(High_Values,nrow(ex_matrix))
            }else{
              ex_matrix[tr,]=rep(High_Values,ncol(ex_matrix))
            }
            Total_alloc=Total_alloc+1
          }

          output=list()
          output$Alloc_Matrix=Alloc_Matrix
          output$Total_Cost=Total_Cost
          
          #If Supply and Demand are not same
          if(sum(Demand) != 0)
            output$Dummy_demand=sum(Demand)
          else if(sum(Supply) != 0)
            output$Dummy_supply=sum(Supply)
          
          if(Total_alloc < (nrow(Alloc_Matrix) + ncol(Alloc_Matrix)-1))
            warning("Degenracy in Transporation Problem Occurred")
          return(output)
}
