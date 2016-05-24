`applytaper` <-
function(f, p=0.05)
  {
    ##################  apply a cosine taper to a time series
    if(missing(p)) { p = 0.1 }
    ###################################   default is a 10 percent taper
    n = length(f)
    
    L=round((n)*p)

    s = seq(from=1, to=n, by=1)

    vwin=rep(1,times=n)
    
    bend = seq(1,  L )

    nx = pi + (pi) * (bend - 1)/(L-1)
    
    vwin[s<=L] = 0.5*cos(nx)+0.5

    bend = seq(n-L+1, n )

    nx = 0 + (pi) * (bend - ( n-L+1  ) )/(L-1)
    
    vwin[s>= n-L+1] = 0.5*cos(nx)+0.5

    newf = vwin*f
    return(newf)
}

