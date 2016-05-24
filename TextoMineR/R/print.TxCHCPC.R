print.TxCHCPC <-
function (x, file = NULL, sep = ";", ...) 
{
    res.TxCHCPC <- x
    if (!inherits( res.TxCHCPC, "TxCHCPC")) 
        stop("non convenient data")
    cat("**Results for the Constrained hierarchical clustering**\n")
    res <- array("", c(24, 2), list(1:24, c("name", "description")))
    res[1, ] <- c("$data.clust", "dataset with the cluster of the individuals")
    res[2, ] <- c("$HierWord", "characteristic words for every node of the hierarchy")
    res[3, ] <- c("$desc.var", "description of the clusters by the variables")
   indice <- 4
    if (!is.null(res.TxCHCPC$desc.var$quanti.var)) {
        res[indice, ] <- c("$desc.var$quanti.var", "description of the cluster var. by the continuous var.")
        res[indice + 1, ] <- c("$desc.var$quanti", "description of the clusters by the continuous var.")
        indice <- indice + 2
    }
    if (!is.null(res.TxCHCPC$desc.var$test.chi2)) {
        res[indice, ] <- c("$desc.var$test.chi2", "description of the cluster var. by the categorical var.")
        res[indice + 1, ] <- c("$desc.axes$category", "description of the clusters by the categories.")
        indice <- indice + 2
    }
    res[indice, ] <- c("$desc.axes", "description of the clusters by the dimensions")
    indice <- indice + 1
    if (!is.null(res.TxCHCPC$desc.axes$quanti.var)) {
        res[indice, ] <- c("$desc.axes$quanti.var", "description of the cluster var. by the axes")
        res[indice + 1, ] <- c("$desc.axes$quanti", "description of the clusters by the axes")
        indice <- indice + 2
    }
    res[indice, ] <- c("$desc.ind", "description of the clusters by the individuals")
    res[indice + 1, ] <- c("$desc.ind$para", "parangons of each clusters")
    res[indice + 2, ] <- c("$desc.ind$dist", "specific individuals")
    res[indice + 3, ] <- c("$call", "summary statistics")
    res[indice + 4, ] <- c("$call$t", "description of the tree")
    indice <- indice + 4
    print(res[1:indice, ])
    if (!is.null(file)) {
        write.infile(res.TxCHCPC, file = file, sep = sep)
        print(paste("All the results are in the file", file))
    }
}
