# This function censors TSPOT result values greater than 20 spots to 20 spots,
# per Oxford Immunotec's instructions

tspot.cens <- function(x){
    if(any(x > 20, na.rm = TRUE)){

        x.cens <- x
        x.cens[x.cens > 20] <- 20

        warning("One or more values were greater than 20 spots and have been censored to 20 spots.")

        return(x.cens)

    } else return(x)
}
