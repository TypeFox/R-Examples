`cusp.Like` <-
function (p, x) 
{
    print(p)
    -2 * sum(log(dcusp((x - p[3])/p[4], p[1], p[2])/p[4]))
}

