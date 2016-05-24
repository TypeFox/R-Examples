########## R function: transUnitInt ##########

# For transformation to the unit interval, with
# 2% spacing at ends. This means the actual
# interval is (0.02,0.98).

tranUnitInt <- function(x,a,b)
   return((b-49*a+48*x)/(50*(b-a)))

############ End of transUnitInt ############
