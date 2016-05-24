
sortmat <-
function (Mat, Sort, decreasing=FALSE){
# useful function for genes ranking.  
#  Sort matrix or dataframe 'Mat', by column(s) 'Sort'.
#  e.g. sortmat(datafr, c(2,4,1))     
# Sort is a number ! 
    m <- do.call("order", c(as.data.frame(Mat[, Sort]), decreasing=decreasing) )
    Mat[m, ]
}
