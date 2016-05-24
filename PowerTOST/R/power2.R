#-----------------------------------------------------------------------------
# Function(s) for power calculation for unbalanced (sequence) groups
# 
# Author: dlabes
#-----------------------------------------------------------------------------
# Jan 2015: power.TOST itself now handles balanced and unbalanced designs
# thus this function is depreciated
power2.TOST <- function(alpha=0.05, logscale=TRUE, theta1, theta2, theta0,
                        CV, n, design="2x2", method="exact", robust=FALSE)
{
  warning("This function is depreciated. User power.TOST() instead.")
  # simply call power.TOST
  pow <- power.TOST(alpha=alpha, logscale=logscale, theta1=theta1, theta2=theta2,
                    theta0=theta0, CV=CV, n=n, design=design, method=method,
                    robust=robust)
  return(pow)
}
