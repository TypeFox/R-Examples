"superLogistic" <-
function (sumT, par, lambda) 
{
    ## par is a list, elements theta_1, theta_2, theta_3 
   
    for (i in 1:length(par)) {
        theta_1 <- par[[i]][1]
	theta_2 <- par[[i]][2] 
	theta_3 <- par[[i]][3]
	sumT <- sumT + (lambda + theta_1) / (theta_2 + exp( - theta_3) )
    }
    sumT
}
