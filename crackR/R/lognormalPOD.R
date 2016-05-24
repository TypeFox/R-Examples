lognormalPOD <-
function(a, median, slope, a.min.detectable=0, poi=1, far=0)
{
    ## lognormal CDF pod curve

    ## false alarm rate parameter is simply a minimum return value

    pcd <- plnorm(a, meanlog=log(median), sdlog=slope)

    ## zero detectability region
    pcd[a < a.min.detectable] <- 0

    ## downweighting all detections by probability of inspection
    pcd <- pcd * poi

    ## apply false alarm rate
    pcd[ pcd < far ] <- far
    
    return(pcd)
}
