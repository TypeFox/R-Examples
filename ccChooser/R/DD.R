DD <-
function(CC,EC){

        if (ncol(CC)!=ncol(EC))
                warning("Erorr: Diferent numbers of traits")

        dist_EC <- mean(daisy(EC))
        dist_CC <- mean(daisy(CC))

        DD_index <- (abs(dist_EC - dist_CC)/dist_EC)*100

        return(DD_index)
}
