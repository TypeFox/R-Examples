 # function to create a list which will contain the solution(s) and the essential prime implicants

`writeSolution` <- 
function(sol, pic) {
    solution <- output <- list()
    row.matrix <- matrix(FALSE, nrow=nrow(pic), ncol=ncol(sol))
    rownames(row.matrix) <- rownames(pic)
    
    for (i in seq(ncol(sol))) {
        aa <- sol[, i]
        aa <- aa[!is.na(aa)]
        row.matrix[aa, i] <- TRUE
    }
    
    ess.PIs <- rownames(pic)[rowSums(row.matrix) == ncol(row.matrix)]
    
    for (i in seq(ncol(sol))) {
        aa <- sol[, i]
        aa <- aa[!is.na(aa)]
        if (length(ess.PIs) > 0) {
            solution[[i]] <- c(ess.PIs, aa[!(aa %in% ess.PIs)])
        }
        else {
            solution[[i]] <- aa
        }
    }
    output[[1]] <- lapply(solution, as.vector)
    output[[2]] <- as.vector(ess.PIs)
    return(output)
}

