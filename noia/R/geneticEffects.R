geneticEffects <-
function (obj, reference = "P1", ref.genotype = NULL) 
{
    if (!is.null(ref.genotype)) {
        warning("The use of ref.genotype is obsolete. Use reference= instead")
        reference = ref.genotype
    }
    if (class(obj) == "noia.linear") {
        new.smat <- genZ2S(obj$genZ, reference = reference)
        new.smat <- solve(new.smat)
        new.smat <- new.smat[colnames(obj$smat), ]
        T <- new.smat %*% obj$smat
        effects <- T %*% obj$E
        std.err <- sqrt((T * T) %*% (obj$std.err * obj$std.err))
        ans <- cbind(effects, std.err)
        colnames(ans) <- c("Effects", "Std.err")
        return(ans)
    }
    else if (class(obj) == "noia.multilinear") {
        stop("Change of reference for the multilinear model: not implemented.")
    }
    else {
        stop("Object of class ", class(obj), " unknown.")
    }
}
