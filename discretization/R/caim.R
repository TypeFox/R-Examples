caim <-
function(tb){
        nr <- dim(tb)[1]
        nc <- dim(tb)[2]
        maxr <- apply(tb,1,max) ## max_j of n_ij, i=1,...,nr
        Mr <- apply(tb,1,sum) ## sum_j n_ij
        sum(maxr^2/Mr)/nr
}
