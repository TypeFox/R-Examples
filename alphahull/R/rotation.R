rotation <-
function (v, theta) 
{
    v.rot <- numeric(2)
    v.rot[1] <- cos(theta) * v[1] + sin(theta) * v[2]
    v.rot[2] <- -sin(theta) * v[1] + cos(theta) * v[2]
    return(v.rot)
}
