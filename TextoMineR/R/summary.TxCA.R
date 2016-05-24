summary.TxCA <-
function (object, nb.dec = 3, nEig=5, ordDim=1, order=FALSE ,nbelements = 10, nbind=nbelements, ncp = 3, align.names = TRUE, nword=10,
    file = "", ...) 
{
 cat("\nCorrespondence analysis summary\n") 

    print2 <- function(mat, file = "") {
        if (file == "") 
            print(mat, quote = FALSE, right = TRUE)
        else {
            mat <- cbind(format(rownames(mat)), mat)
            mat <- rbind(colnames(mat), mat)
            mat2 <- cbind(format(mat[, 1], justify = "right"), 
                format(mat[, 2], justify = "right"))
            for (k in 3:ncol(mat)) mat2 <- cbind(mat2, format(mat[, 
                k], justify = "right"))
            mat2 <- cbind(mat2, "\n")
            for (i in 1:nrow(mat2)) cat(mat2[i, ], file = file, 
                append = TRUE)
        }
    }
    print3 <- function(obj, file = "", ncp, width.row = 0, nbelements = nbelements) {
        list.obj <- match.arg(names(obj), c("dist", "inertia", 
            "coord", "cos2", "contrib", "v.test", "vtest"), several.ok = TRUE)
        nb.col <- sum(c("coord", "cos2", "contrib", "v.test") %in% 
            list.obj)
        nbelements <- min(nbelements, nrow(obj$coord))
        mat <- matrix(NA, nbelements, sum(c("dist", "inertia") %in% 
            list.obj) + nb.col * ncp)
        if(order){
		if("contrib" %in% list.obj){
			obj$coord<-obj$coord[(rev(order(obj$contrib[,ordDim] ))),]
			obj$cos2<-obj$cos2[(rev(order(obj$contrib[, ordDim] ))),]
			obj$contrib<-obj$contrib[(rev(order(obj$contrib[,ordDim]))),]
        	}}
		colnames(mat) <- paste("v", 1:ncol(mat))
        	rownames(mat) <- format(rownames(obj$coord)[1:nbelements, drop = FALSE], width = width.row)
        	indice <- 1
		if ("coord" %in% list.obj) {
            	mat[, indice + nb.col * (0:(ncp - 1))] <- obj$coord[1:nbelements, 1:ncp, drop = FALSE]
            	colnames(mat)[indice + nb.col * (0:(ncp - 1))] <- paste("Dim.", 1:ncp, sep = "")
            	indice <- indice + 1
        	}
		if ("cor" %in% list.obj) {
            	mat[, indice + nb.col * (0:(ncp - 1))] <- obj$cor[1:nbelements, 1:ncp, drop = FALSE]
            	colnames(mat)[indice + nb.col * (0:(ncp - 1))] <- paste("cor", 1:ncp, sep = "")
            	indice <- indice + 1
        	}
        	if ("contrib" %in% list.obj) {
            	mat[, indice + nb.col * (0:(ncp - 1))] <- obj$contrib[1:nbelements, 1:ncp, drop = FALSE]
            	colnames(mat)[indice + nb.col * (0:(ncp - 1))] <- "ctr"
            	indice <- indice + 1
        	}
        	if ("cos2" %in% list.obj) {
            	mat[, indice + nb.col * (0:(ncp - 1))] <- obj$cos2[1:nbelements, 1:ncp, drop = FALSE]
            	colnames(mat)[indice + nb.col * (0:(ncp - 1))] <- "cos2"
            	indice <- indice + 1
        	}
            mat <- format(round(mat, nb.dec))
     		mat2 <- "|"
        	for (k in 1:ncp) mat2 <- cbind(mat2, mat[, (1 + nb.col * (k - 1)):(nb.col * k), drop = FALSE], "|")
        	colnames(mat2)[1] <- ""
        	print2(as.matrix(mat2), file = file)
    	}
     res1 <- object
     if (!inherits(res1, "TxCA")) 
        stop("non convenient object")
       cat("\n", file = file, append = TRUE)
	 print2(res1$TableSummary, file = file)
    res<-res1$res.ca
    cat(paste("\nCall:\n"), file = file)
    cat(paste(deparse(res$call$call), "\n"), file = file, append = TRUE)
    cat("\n", file = file, append = TRUE)
    cat("\nEigenvalues\n", file = file, append = TRUE)
    eige <- format(t(round(res$eig[, 1:3], nb.dec)), justify = "right")
    rownames(eige) <- c("Variance", "% of var.", "Cumulative % of var.")
    nEig<- min(nEig,nrow(res$eig))
    colnames(eige[,1:nEig]) <- paste("Dim", 1:nEig, sep = ".")
    	print2(eige[,1:nEig], file = file)
      cat("\n", file = file, append = TRUE)
	print(res1$Inertia.VCr, file = file)
	 width.row <- 0
    if (align.names == TRUE) {
        aux <- match.arg(names(res), c("ind", "ind.sup", "freq", 
            "freq.sup", "var", "quanti.var", "quanti.var.sup", 
            "quali.var", "quali.var.sup", "quanti.sup", "quali.sup", 
            "group", "row", "row.sup", "col", "col.sup"), several.ok = TRUE)
        width.row = max(nchar(rownames(res[aux[1]][[1]]$coord)))
        for (k in 1:length(aux)) width.row = max(width.row, nchar(rownames(res[aux[k]][[1]]$coord)
           [1:min(nrow(res[aux[k]][[1]]$coord), nbelements)]))
    }
    ncp <- min(res$call$ncp, ncp)
   if (nbind > 0) {
    cat("\nRows", file = file, append = TRUE)
    if (nrow(res$row$coord) > nbind)
   if(order){ 
        cat(paste(" (the ", nbind, " first most contributed to the)", sep = ""), 
            file = file, append = TRUE)
      cat("Dimension",ordDim,"\n")
     }else{
 	cat(paste(" (the ", nbind, " first)", sep = ""), 
            file = file, append = TRUE)
	}
    cat("\n", file = file, append = TRUE)
    print3(res$row, file = file, ncp = ncp, width.row = width.row, 
        nbelements = nbind)
   }
    cat("\nColumns", file = file, append = TRUE)
    if (nrow(res$col$coord) > nbelements)
	if(order){ 
        cat(paste(" (the ", nbelements, " first most contributed to the)", sep = ""), 
            file = file, append = TRUE)
	cat("Dimension",ordDim,"\n")
      }else{
          cat(paste(" (the ", nbelements, " first)", sep = ""), 
            file = file, append = TRUE)
      }
    cat("\n", file = file, append = TRUE)
    print3(res$col, file = file, ncp = ncp, width.row = width.row, 
        nbelements = nbelements)
    if (!is.null(res$row.sup)) {
        cat("\nSupplementary row", file = file, append = TRUE)
        if (nrow(res$row.sup$coord) > 1) 
            cat("s", file = file, append = TRUE)
        if (nrow(res$row.sup$coord) > nbelements) 
            cat(paste(" (the ", nbelements, " first)", sep = ""),  file = file, append = TRUE)
        cat("\n", file = file, append = TRUE)
        print3(res$row.sup, file = file, ncp = ncp, width.row = width.row, 
            nbelements = nbelements)
    }
    if (!is.null(res$col.sup)) {
        cat("\nSupplementary column", file = file, append = TRUE)
        if (nrow(res$col.sup$coord) > 1) 
            cat("s", file = file, append = TRUE)
        if (nrow(res$col.sup$coord) > nbelements) 
            cat(paste(" (the ", nbelements, " first)", sep = ""), 
                file = file, append = TRUE)
        cat("\n", file = file, append = TRUE)
        print3(res$col.sup, file = file, ncp = ncp, width.row = width.row, 
            nbelements = nbelements)
    }
    if (!is.null(res$quanti.var)) {
        cat("\nSupplementary continuous variable", file = file, 
            append = TRUE)
        if (nrow(res$col.sup$coord) > 1) 
            cat("s", file = file, append = TRUE)
        cat("\n", file = file, append = TRUE)
        print3(res$quanti.var, file = file, ncp = ncp, width.row = width.row, 
            nbelements = nbelements)
     }
   if (!is.null(res1$res.agg)) {

    cat("\nAggregation of documents according to the categorical variable\n", file = file, append = TRUE)
	print2(res1$res.agg$TableLexCat, file = file)
    cat("\nDistribution of the words into the categories\n", file = file, append = TRUE)
	print2(res1$res.agg$TableDistForm, file = file)
    }
	cat("\nGlossary of the ",nword," most frequent words after selection\n")
    nword<-min(nword,nrow(res1$Glossary)) 
     print(res1$Glossary[c(1:nword),])
   
}
