`RCRng` <- 
function(n)
    if (length(n) == 1) c(1,n[1]) else c(n[1],n[length(n)])
