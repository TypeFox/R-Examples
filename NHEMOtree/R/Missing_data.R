Missing_data <-
function(data){

  # Missing variables
    mv<- c()
    for (i in 1:nrow(data)) mv[i]<- sum(is.na(data[i,1:ncol(data)]))
    mv2<- which(mv>0)

  # Deletion of rows with missing values
    if (length(mv2)>0) data<- data[-mv2,] 
    return(data)
}
