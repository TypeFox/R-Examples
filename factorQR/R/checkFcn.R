checkFcn <- function(x, p=.5){
return((x < 0)*-1*x*(1-p) + (x>0)*p*x)
 		}