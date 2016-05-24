membercheck <- function(member){
    member <- as.matrix(member)
    if(any(member < 0) || any(member > 1))
        stop("membership values (posterior probabilities)", " must be in [0,1]")
    if (!isTRUE(all.equal(as.vector(rowSums(member)), rep.int(1, nrow(member)))))
        stop("Membership values (posterior probabilities)", " must sum up to 1")
}
