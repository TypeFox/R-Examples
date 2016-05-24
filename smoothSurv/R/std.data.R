###########################################
#### AUTHOR:    Arnost Komarek         ####
####            (2003)                 ####
####                                   ####
#### FILE:      std.data.R             ####
####                                   ####
#### FUNCTIONS: std.data               ####
###########################################

### ============================================================================
### std.data: Standardize data (subtract mean and divide by standard deviation)
### ============================================================================
## datain ..... input dataframe
## cols ....... which cols of the original data set should be transformed
##
## OUTPUT .... dataout --> data.frame
std.data <- function(datain, cols){
   dataout <- datain
   changecols <- colnames(dataout) %in% cols
   leavecols <- !changecols

   options(warn = -1)
   means <- sapply(dataout, mean, na.rm = TRUE)
   sds <- sapply(dataout, sd, na.rm = TRUE)

   options(warn = 1)
   changed <- 0
   for(i in 1:ncol(dataout)){
     if(changecols[i]){
        if(is.na(means[i]) | is.na(sds[i])){
           str <- paste("Missing mean or sd for variable ", colnames(dataout)[i],
            ", it is not standardized.", sep="")
           warning(str, call.=FALSE)
        }
        else{
           dataout[[i]] <- (dataout[[i]] - means[i])/sds[i]
           changed <- changed + 1
        }
     }
   }

   options(warn = 0)     ## default value
   cat("\nNumber of standardized columns: ", changed, "\n")

   means <- means[changecols]
   sds <- sds[changecols]
   tab <- rbind(means, sds)
   rownames(tab) <- c("mean","sd")

   cat("\nUsed means and sd's: \n")
   print(tab)

   return(dataout)
}

