reorgRIgenoprob <- 
function (cross) 
{
    crosses <- cross$cross
    flag <- 0
    for (i in 1:ncol(crosses)) {
        if (any(crosses[, i] != i)) {
            flag <- 1
            break
        }
    }
    if (!flag) 
        return(cross)
    crosstype <- class(cross)[1]
    if (crosstype != "ri4sib" && crosstype != "ri4self" && crosstype != 
        "ri8sib" && crosstype != "ri8self") 
        stop("reorgRIgenoprob not appropriate for cross type ", 
            crosstype)
    n.str <- as.numeric(substr(crosstype, 3, 3))
    n.ind <- nind(cross)
    for (i in names(cross$geno)) {
        chrtype <- class(cross$geno[[i]])
        if (chrtype == "X") 
            warning("reorgRIgenoprob not working properly for the X chromosome.")
        if (!("prob" %in% names(cross$geno[[i]]))) {
            warning("No QTL genotype probabilities within cross.")
            return(cross)
        }
        prob <- cross$geno[[i]]$prob
        att <- attributes(prob)
        n.mar <- dim(prob)[2]
        if (dim(prob)[1] != n.ind) 
            stop("Mismatch between no. individuals in cross and in genoprobs.")
        if (dim(prob)[3] != n.str) {
            warning("Odd no. columns in genoprobs for chromosome ", 
                i)
            next
        }
        prob <- .C("R_reorgRIgenoprob", as.integer(n.ind), as.integer(n.mar), 
            as.integer(n.str), prob = as.double(prob), as.integer(crosses), 
            PACKAGE = "qtl")$prob
        prob <- array(prob, dim = c(n.ind, n.mar, n.str))
        for (j in names(att)) attr(prob, j) <- att[[j]]
        cross$geno[[i]]$prob <- prob
    }
    cross
}

