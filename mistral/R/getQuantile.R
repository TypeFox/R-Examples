## -----------------------------------------------------------------------------
## Fonction getQuantile
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

getQuantile = function(data,alpha,floor=0) {

    data_sorted = sort(data);
    N = length(data)
    ind = (N-1)*alpha + 1
    y = data_sorted[floor(ind)] + (ind-floor(ind))*(data_sorted[floor(ind)+1]-data_sorted[floor(ind)])
    cat("y =",y,"\n")
    if(y<floor) {y=0;cat("y<",floor," => last subset & y = ",floor,"\n",sep="")}

    return(y)
}