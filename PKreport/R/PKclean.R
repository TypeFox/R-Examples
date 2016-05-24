#############################################################################################
## File: PKclean.R
## Author: Xiaoyong Sun
## Date: 10/22/2009
## Goal: clean all non-use folder and files
## Notes:
##      -
#############################################################################################

## clean all folders and files
## clean hidden code, and code notes
PKclean <- function()
{
    general.list <- .pkplot$getGlobalConfig()
    sapply(5:length(general.list), function(i)
            {
               if(file.exists(general.list[[i]])) unlink(general.list[[i]], recursive=TRUE)
            })
    .pkplot$cleanPKCode()
    .pkplot$cleanPKCodeNote()
    

}

