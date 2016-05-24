sigmashrink <-
function (s2, d) 
{
    df = d
    n = length(s2)
    if (length(d) == 1) 
        d = rep(d, n)
    s2t = rep(0, n)
    d0 = 0
    s20 = 0
    bern = c(1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510, 
        43867/798, -174611/330, 854513/138, -236364091/2730, 
        8553103/6, -23749461029/870, 8615841276005/14322, -7709321041217/510, 
        2577687858367/6, -26315271553053478912/1919190, 2929993913841559/6, 
        -2.61082718496449e+20/13530)
    dicoeffs = tricoeffs = quadcoeffs = rep(0, 20)
    for (i in 1:20) {
        dicoeffs[i] = bern[i]/(2 * i)
        tricoeffs[i] = bern[i]
        quadcoeffs[i] = bern[i] * (2 * i + 1)/2
    }
    gamcoeffs = c(dicoeffs, tricoeffs, quadcoeffs)
    temp = .C("sigmashrink", as.double(s2), as.double(d), as.double(s2t), 
        as.double(d0), as.double(s20), as.double(gamcoeffs), 
        as.integer(n))
    sigma2 = temp[[3]]
    if (temp[[4]] > 0) 
        df = df + temp[[4]]
    else df = Inf
    return(list(sigma2 = sigma2, df = df))
}
