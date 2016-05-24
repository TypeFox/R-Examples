########## R function: intervalGrow ##########

# For growing an interval by adding an
# interval of `radius' r to the ends.

# Last changed: 11 DEC 2008

intervalGrow <- function(sttInterval,r)
{
   grownInterval <- c(sttInterval[1]-2*r,sttInterval[2]+2*r)
  
   return(list(grownInterval=grownInterval,
               volume=abs(grownInterval[2]-grownInterval[1])))
}

############ End of intervalGrow ############
