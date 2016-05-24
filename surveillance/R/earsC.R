#         \|||/
#         (o o)
# ,~~~ooO~~(_)~~~~~~~~~,
# |        EARS        |
# |   surveillance     |
# |       methods      |
# |   C1, C2 and C3    |
# '~~~~~~~~~~~~~~ooO~~~'
#        |__|__|
#         || ||
#        ooO Ooo


######################################################################
# Implementation of the EARS surveillance methods.
######################################################################
# DESCRIPTION
######################################################################
# Given a time series of disease counts per month/week/day
# this function determines whether there was an outbreak at given time points:
# it deduces for each time point an expected value from past values,
# it defines an upperbound based on this value and on the variability
# of past values
# and then it compares the observed value with the upperbound.
# If the observed value is greater than the upperbound 
# then an alert is flagged.
# Three methods are implemented.
# They do not use the same amount of past data 
# and are expected to have different specificity and sensibility
# from C1 to C3
# the amount of past data used increases,
# so does the sensibility
# but the specificity decreases.
######################################################################
# PARAMETERS
######################################################################
#   range : range of timepoints over which the function will look for
# outbreaks.
#   method : which of the three EARS methods C1, C2 and C3 should be used.
#
######################################################################
# INPUT
######################################################################
# A R object of class sts
######################################################################
# OUTPUT
######################################################################
# The same R object of class sts with slot alarm and upperbound filled
# by the function
######################################################################

earsC <- function(sts, control = list(range = NULL, method = "C1",
                                      alpha = 0.001)) { 

######################################################################
  #Handle I/O
  ######################################################################
  
  #If list elements are empty fill them!

  # Method
  if (is.null(control[["method",exact=TRUE]])) {
    control$method <- "C1"
  }
  
  # Extracting the method
  method <- match.arg( control$method, c("C1","C2","C3"),several.ok=FALSE)
  
  # Range
  # By default it will take all possible weeks
  # which is not the same depending on the method
  if (is.null(control[["range",exact=TRUE]])) {
    if (method == "C1"){
      control$range <- c(8:dim(sts@observed)[1])
    }
    if (method == "C2"){
      control$range <- c(10:dim(sts@observed)[1])
    }
    if (method == "C3"){
      control$range <- c(12:dim(sts@observed)[1])
    }    
  }
  
  # zAlpha
  if (is.null(control[["alpha",exact=TRUE]])) {
    # C1 and C2: Risk of 1st type error of 10-3
    # This corresponds to an Z(1-zAlpha) of about 3
    if (method %in% c("C1","C2")) {
      control$alpha = 0.001
  }
    # C3: Risk of 1st type error of 0.025
    # This corresponds to an Z(1-zAlpha) of about 2
    if (method=="C3") {
      control$alpha = 0.025
    }
  }

  # Calculating the threshold zAlpha
  zAlpha <- qnorm((1-control$alpha))
 
  
  #Deduce necessary amount of data from method
  maxLag <- switch(method, C1=7, C2=9, C3=11)
  
  # Order range in case it was not given in the right order
  control$range = sort (control$range)

  ######################################################################
  #Loop over all columns in the sts object
  #Call the right EARS function depending on the method chosen (1, 2 or 3)
  #####################################################################
   for (j in 1:ncol(sts)) {

     # check if the vector observed includes all necessary data: maxLag values.
     if((control$range[1] - maxLag) < 1) { 
       stop("The vector of observed is too short!")
     }      

     ######################################################################
     # Method C1 or C2
     ######################################################################
     if (method %in% c("C1","C2")) {

       # Create a matrix with time-lagged vectors
       refVals <- NULL
       for (lag in maxLag:(maxLag-6)) {
         refVals <- cbind(refVals, observed(sts)[(control$range-lag),j])
       }

       # calculate the upperbound 
       sts@upperbound[control$range,j] <- apply(refVals,1,mean)+
        zAlpha*apply(refVals,1,sd)
     }
          
     if (method=="C3") {
       # Create a matrix with time-lagged vectors
       refVals <- NULL
       rangeC2 = ((min(control$range) - 2) : max(control$range))
       for (lag in 9:3) {
         refVals <- cbind(refVals, observed(sts)[(rangeC2-lag),j])
       }

       # Calculate C2
       C2 <- (observed(sts)[rangeC2,j] - apply(refVals,1,mean)) / apply(refVals,1,sd)
       # Calculate the upperbound
       # first calculate the parts of the formula with the maximum of C2 and 0 for     # two time lags.
       partUpperboundLag2 =  pmax(rep(0,length=length(C2)-2),C2[1:(length(C2)-2)]-1)
       partUpperboundLag1 =  pmax(rep(0,length=length(C2)-2),C2[2:(length(C2)-1)]-1)
     
       sts@upperbound[control$range,j] <- observed(sts)[control$range,j] +
         apply(as.matrix(refVals[3:length(C2),]),1,sd)*(zAlpha - (partUpperboundLag2 + partUpperboundLag1))

       # Upperbound must be superior to 0 which is not always the case
       #with this formula
       sts@upperbound[control$range,j] = pmax(rep(0,length(control$range)),sts@upperbound[control$range,j])     
     }  # end of loop over j
   } #end of loop over cols in sts

  #Make sts return object
  control$name <- paste("EARS_",method,sep="")
  control$data <- paste(deparse(substitute(sts)))
  sts@control <- control
  #Where are the alarms?
  sts@alarm[control$range,] <- matrix(observed(sts)[control$range,]>upperbound(sts)[control$range,] )
  
  #Done
  return(sts[control$range,])
}

