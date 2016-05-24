RESCALE<-function (x, nx1=0, nx2=1, minx=0, maxx=1) 
{
    if(missing(minx)) minx = min(x)
    if(missing(maxx)) maxx = max(x)
    
    if(minx == maxx)
        {
            nx = rep(mean(c(nx1, nx2)) , length(x) )
         }
        else
            {
                nx = nx1 + (nx2 - nx1) * (x - minx)/(maxx - minx)
            }
    return(nx)
}

