####################
#This function calculates the variance of a Fisher Z transformed estimate x, given its original variance 
#########################
FisherZ.Var <-
function(x,var.x)
  {
    var.x/(1-x^2)^2
    
  }

