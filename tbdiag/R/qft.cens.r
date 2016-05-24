# This function censors QFT result values greater than 10 IU/mL to 10 IU/mL,
# per Cellestis' instructions

qft.cens <- function(x){
    # If any values are greater than 10, censor those values to 10
    if(any(x > 10, na.rm = TRUE)){

        x.cens <- x
        x.cens[x.cens > 10] <- 10

        warning("One or more values were greater than 10 IU/mL and have been censored to 10 IU/mL.")

        return(x.cens)

    } else return(x)
}
