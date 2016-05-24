anglesArc <-
function (v, theta) 
{
    theta.OX <- ifelse(v[2] >= 0, acos(v[1]), 2 * pi - acos(v[1]))
    angs <- c(theta.OX - theta, theta.OX + theta)
    return(angs)
}
