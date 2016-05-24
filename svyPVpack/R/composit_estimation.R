composit <- function(inputlist, labcolumn, statcolumn, weightcolumn, SEcolumn=NULL)
{
# compute the composit estimate from any statistic
# inputlist = is list of tables which stems from the repeated application of a svyPV function. each element of the list must contain an data.frame with the following columns:
  ##1) a column which labels the rows of the data.frame eg.: "Group1" which perhaps indicates "Gender" and contains values like Male/Female etc.
  ##2) a column with the statistic of interest, e.g.: mean or sd
  ##3) a column with weights
  #### ----> typically all these values are provided by the functions of the svyPVpack package

#   inputlist    <- pmlist
#   labcolumns   <- "Group1"
#   statcolumn   <- "PVPSL_mean"
#   weightcolumn <- "Sum.of.weights"  
#   SEcolumn=NULL
#   SEcolumn     <- "PVPSL_mean.SE"
  
# wie wird aber mit mehreren gruppen umgegangen
  
####### controls ################################  

NA11 <- c(labcolumn,statcolumn,weightcolumn)

check1 <- sapply(inputlist, function(x)
            {
              sum(NA11 %in% colnames(x)  )
            })

if(!all(length(NA11) == check1))
{
stop(paste("Check names in data.frames in inputlist - list element: ",which(length(NA11) != check1),"\n"))
}


#################################################  
  
inl <- lapply(inputlist, function(x)
            {
            x[,c(weightcolumn,statcolumn)]  
            })  
  
  
inl_array <- array(unlist(inl), dim=c(dim(inl[[1]]), length(inl)), dimnames = list(rownames(inputlist[[1]]),c(weightcolumn,statcolumn)))
  
# compute the composite statistic - multi group input is possible - an composit estimate for each group is returned.
comp_stat <- apply(inl_array,1, function(x) sum(x[1,] * x[2,])/sum(x[1,]))
#sum(inl_array[,1,]  * inl_array[,2,]) / sum(inl_array[,1,])
  

######## if a standard error of the component is desired
if(!is.null(SEcolumn))
  {
    
  sel <- lapply(inputlist, function(x)
      {
        x[,c(weightcolumn,SEcolumn)]  
      })    
      
  sel_array <- array(unlist(sel), dim=c(dim(sel[[1]]), length(sel)), dimnames = list(rownames(inputlist[[1]]),c(weightcolumn,SEcolumn)))
  
  # compute the standard error of the composit statistic
  comp_se <- apply(sel_array,1,function(x)
                    {
                    sqrt(sum(x[1,]^2*x[2,]^2))/sum(x[1,])
                    })

  NAMES <- paste0("composite_",c(statcolumn,SEcolumn))
  } else { 
         comp_se <- NA 
         NAMES <- paste0("composite_",c(statcolumn,"SE"))
         }




comp_allp <- data.frame(comp_stat,comp_se)
colnames(comp_allp) <- NAMES


# add group information

labcol <- data.frame(inputlist[[1]][,labcolumn])
colnames(labcol) <- labcolumn

comp_all <- data.frame(labcol,comp_allp)
  
return(comp_all) 
  
}




