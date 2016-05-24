correctVarNames <- function(tab, rowcol = TRUE, cols = 1:ncol(tab)){

# row- and col-names
if (rowcol == TRUE){
for (c in 1:dim(tab)[1]){dimnames(tab)[[1]][c] <- gsub("_", "\\\\_", as.character(dimnames(tab)[[1]][c]))}
for (c in 1:dim(tab)[2]){dimnames(tab)[[2]][c] <- gsub("_", "\\\\_", as.character(dimnames(tab)[[2]][c]))}
}

## to replace a "." with a blank: colnames(vars) <- gsub("\\.", " ", colnames(vars))   

tab2 <- matrix(NA, ncol = ncol(tab), nrow = nrow(tab))
tab2 <- data.frame(tab2)
dimnames(tab2) <- dimnames(tab)

repl.cols <- (1:ncol(tab))[((1:ncol(tab)) %in% cols) == FALSE]
tab2[, repl.cols] <- tab[, repl.cols]

# entries
if (identical(cols, FALSE)){
    for (j in cols){
        col.j <- as.character(tab[, j])
    
        for (i in 1:nrow(tab)){
            if (is.na(col.j[i]) == FALSE){col.j[i] <- gsub("_", "\\\\_", as.character(col.j[i]))}
        } # end i

    tab2[, j] <- col.j
    } 
}

return(tab2)
}
