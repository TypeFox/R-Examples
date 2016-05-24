displayKbyC <- function(v1, v2, percentage = c("none", "row", "col", "total")[1], names = c("v1", "v2"), cap = "", lab = "", row.nam = NA, col.nam = NA){
    if (is.factor(v1)){v1 <- factor(v1, exclude = NULL)}
    if (is.factor(v2)){v2 <- factor(v2, exclude = NULL)}
    mat <- table(v1, v2)
    if (is.na(row.nam[1]) == TRUE) {row.nam <- dimnames(mat)[[1]]}
    if (is.na(col.nam[1]) == TRUE) {col.nam <- dimnames(mat)[[2]]}
    s1 <- as.vector(apply(mat, 2, sum))
    s2 <- as.vector(apply(mat, 1, sum))
    s3 <- sum(mat)
    n.v1 <- length(unique(v1[is.na(v1) == FALSE]))
    n.v2 <- length(unique(v2[is.na(v2) == FALSE]))
    mat <- rbind(col.nam, mat, s1)
    mat <- cbind(c(NA, row.nam, "Total"), mat)
    mat <- cbind(mat, c("Total", s2, s3))
    mat <- rbind(c(NA, names[2], rep(NA, n.v2)), mat)
    mat <- cbind(c(NA, NA, names[1], rep(NA, n.v1)), mat)

    ## add row-wise percentages
    if (identical(percentage, "row")){
        for (r in 3:nrow(mat)){
            ro <- as.numeric(mat[r, 3:ncol(mat)])
            ro <- paste(ro, " (", disp(100 * ro / max(ro), 1), "\\%)", sep = "")
            ro <- sub("\\( ", "\\(", ro)
            mat[r, 3:ncol(mat)] <- ro
        }
    }
    
    ## add col-wise percentages
    if (identical(percentage, "col")){
        for (r in 3:ncol(mat)){
            ro <- as.numeric(mat[3:nrow(mat), r])
            ro <- paste(ro, " (", disp(100 * ro / max(ro), 1), "\\%)", sep = "")
            ro <- sub("\\( ", "\\(", ro)
            mat[3:nrow(mat), r] <- ro
        }
    }
    
    ## add col-wise percentages with respect to total number of observations
    if (identical(percentage, "total")){
        n <- sum
        for (r in 3:ncol(mat)){
            ro <- as.numeric(mat[3:nrow(mat), r])
            ro <- paste(ro, " (", disp(100 * ro / s3, 1), "\\%)", sep = "")
            ro <- sub("\\( ", "\\(", ro)
            mat[3:nrow(mat), r] <- ro
        }
    }
    
    mat <- data.frame(mat)

    ali <- "lll|"
    for (i in 1:n.v2) {ali <- paste(ali, "c", sep = "")}
    xtab2 <- xtable(mat, align = paste(ali, "|c", sep = ""), caption = cap, label = lab)
    print(xtab2, include.rownames = FALSE, include.colnames = FALSE, floating = FALSE, hline.after = c(2, 2 + n.v1), type = "latex", 
        size = "footnotesize", sanitize.text.function = function(x){x}, table.placement = "h!", tabular.environment = "longtable")        
}
