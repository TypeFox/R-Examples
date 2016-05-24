###########################################################################################
##  Funcion de impresión elementos tipo 'FactoClass'                                     ##
##                                                                                       ##
## Elaborado por: Pedro Cesar del Campo Neira                                            ##
## Universidad Nacional de Colombia                                                      ##
##                                                                                       ##
## requiere:ade4      library(ade4)                                                      ##
##                                                                                       ##
## print.FactoClass  ( x = elemento tipo 'FactoClass')                                   ##
##                                                                                       ##
###########################################################################################


print.FactoClass <- function(x, ...) 
{   
    cat("\n")
    cat(" FactoClass: combination of factorial methods and cluster analysis\n\n")

    cat("--------------------------------------------------------------------------\n")

    cat("Object $dudi (Factorial analyis) \n\n")

    cat("Duality diagramm\n")
    cat("class: ")
    cat(class(x$dudi))
    cat("\n$dudi$call: ")
    print(x$dudi$call)
    cat("\n$dudi$nf:", x$dudi$nf, "axis-components saved")
    cat("\n$dudi$rank: ")
    cat(x$dudi$rank)
    cat("\neigen values: ")
    l0 <- length(x$dudi$eig)
    cat(signif(x$dudi$eig, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n")
    else cat("\n")
    sumry <- array("", c(3, 4), list(1:3, c("vector", "length", 
        "mode", "content")))
    sumry[1, ] <- c("$dudi$cw",  length(x$dudi$cw),  mode(x$dudi$cw),  "column weights")
    sumry[2, ] <- c("$dudi$lw",  length(x$dudi$lw),  mode(x$dudi$lw),  "row weights")
    sumry[3, ] <- c("$dudi$eig", length(x$dudi$eig), mode(x$dudi$eig), "eigen values")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(5, 4), list(1:5, c("data.frame", "nrow", 
        "ncol", "content")))
    sumry[1, ] <- c("$dudi$tab", nrow(x$dudi$tab), ncol(x$dudi$tab), "modified array")
    sumry[2, ] <- c("$dudi$li",  nrow(x$dudi$li),  ncol(x$dudi$li),  "row coordinates")
    sumry[3, ] <- c("$dudi$l1",  nrow(x$dudi$l1),  ncol(x$dudi$l1),  "row normed scores")
    sumry[4, ] <- c("$dudi$co",  nrow(x$dudi$co),  ncol(x$dudi$co),  "column coordinates")
    sumry[5, ] <- c("$dudi$c1",  nrow(x$dudi$c1),  ncol(x$dudi$c1),  "column normed scores")
    class(sumry) <- "table"
    print(sumry)
    cat("other elements: ")
    if (length(names(x$dudi)) > 11) 
        cat(paste("$dudi$",names(x$dudi)[12:(length(x$dudi))],sep=""), "\n")
    else cat("NULL\n")
    cat("--------------------------------------------------------------------------\n")
    
    cat(" Number of axes for cluster:  ")
    cat(x$nfcl)
    
    cat("\n Number of clusters:  ")
    cat(x$k,"\n\n")
    
    sumry <- array("", c(7,2), list(1:7, c("Object", "Description")))
    sumry[1, ] <- c("$indices"    ,"level indices for Hierarchical Clustering (WARD)")
    sumry[2, ] <- c("$cor.clus"   ,"centroid cluster coordinates")
    sumry[3, ] <- c("$clus.summ"  ,"partition changes due to consolidation process")
    sumry[4, ] <- c("$cluster"    ,"a vector indicating the cluster in which each point is allocated")
    sumry[5, ] <- c("$carac.cate" ,"cluster characterization by qualitative variables")
    sumry[6, ] <- c("$carac.cont" ,"cluster characterization by  quantitative variables")
    sumry[7, ] <- c("$carac.frec" ,"cluster characterization by  frequency active variables")

    class(sumry) <- "table"
    print(sumry)
   
}

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
