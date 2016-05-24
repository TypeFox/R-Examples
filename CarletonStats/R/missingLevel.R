missingLevel <-
function(data){
 if (!is.data.frame(data)) stop("Data must be a data frame")
 m <- dim(data)[2]
 for (i in 1:m)
 {
   if(is.factor(data[,i]))
    { if (levels(data[,i])[1]=="")
       is.na(levels(data[,i])[1]) <- TRUE
    }
 }

return(data)
}
