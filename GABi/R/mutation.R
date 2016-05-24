mutation <- function(pop,mfreq){

        changes <- array(as.vector(runif(pop)<mfreq,mode="numeric"),dim=dim(pop));
        array(as.vector(xor(pop,changes),mode="numeric"),dim=dim(pop));

}

