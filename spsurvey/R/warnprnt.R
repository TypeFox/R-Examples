warnprnt <- function(warn.df=get("warn.df", envir=.GlobalEnv),
   m=1:nrow(warn.df)) {

################################################################################
# Function: warnprnt
# Programmer: Tom Kincaid
# Date: September 19, 2005
# Revised: April 17, 2014
# Description:
#   This function prints the warnings data frame. 
#   Input:
#      warn.df = a data frame that contains warning messages.  The default is
#         "warn.df", which is the name given to the warnings data frame created
#          by functions in the spsurvey package.
#      m = the vector of indices for warning messages that are to be printed.
#         The default is a vector containing the integers from 1 through the
#         number of rows in warn.df, which will print all warning messages in
#         the data frame.
#   Output: None.
################################################################################

   for(i in m) {
      cat(paste("Warning Message", i, "\n"))
      cat(paste("Function:", warn.df$func[i], "\n"))
      if(!is.na(warn.df$subpoptype[i]))
         cat(paste("Population Type:", warn.df$subpoptype[i], "\n"))
      if(!is.na(warn.df$subpop[i]))
         cat(paste("Subpopulation:", warn.df$subpop[i], "\n"))
      if(!is.na(warn.df$indicator[i]))
         cat(paste("Indicator:", warn.df$indicator[i], "\n"))
      if(!is.na(warn.df$stratum[i]))
         cat(paste("Stratum:", warn.df$stratum[i], "\n"))
      cat(paste("Warning:", warn.df$warning[i]))
      cat(paste("Action:", warn.df$action[i], "\n"))
   }

   invisible(NULL)
}

