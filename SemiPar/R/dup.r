############ S-function: dup  ############

# For determined whether a value
# in a vector is duplicated

# Last changed: 07 JUL 2001

dup <- function(dat)
{
   dat <- match(dat,unique(dat))
   id <- order(dat)
   look <- c(F, ifelse(diff(sort(dat))==0,T,F))
   dat[id] <- look
     
   return(as.logical(dat))
}

########### End of dup #############
