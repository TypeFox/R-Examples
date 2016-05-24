logctablepost=function (theta, data) 
{
    theta1 = theta[1]
    theta2 = theta[2]
    s1 = data[1]
    f1 = data[2]
    s2 = data[3]
    f2 = data[4]
    logitp1 = (theta1 + theta2)/2
    logitp2 = (theta2 - theta1)/2
    term1 = s1 * logitp1 - (s1 + f1) * log(1 + exp(logitp1))
    term2 = s2 * logitp2 - (s2 + f2) * log(1 + exp(logitp2))
    return(term1 + term2)
}
