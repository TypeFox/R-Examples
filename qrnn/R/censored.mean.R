censored.mean <-
function(x, lower, trim = 0)
{
    x.mean <- mean(x, trim = trim)
    x.median <- median(x)
    if (isTRUE(all.equal(x.median, lower))){
        return(x.median)
    } else{
        return(x.mean)
    }
}

