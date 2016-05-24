.make.diss.names <- function(column.names){

    n <- length(column.names);

    diss.names <- character(n * (n - 1) / 2);

    k <- 1;

    for(i in 1:(n-1)){

        for(j in (i+1):n){

            diss.names[k] <- paste(column.names[c(i,j)],collapse = "v");

            k <- k + 1;
        }
    }

    return(diss.names);
}
