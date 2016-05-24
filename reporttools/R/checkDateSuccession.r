checkDateSuccession <- function(d1, d2, pat, names = NA, lab = "", typ = c("R", "tex")[2]){

    if (is.na(names[1]) == TRUE){names <- c("Observation", "v1", "v2")}
    
    dat <- data.frame(cbind(d1, d2))
    ind <- (apply(abs(dat), 1, sum) > 0)
    ind <- ind * 1:length(d1)
    ind <- ind[is.na(ind) == FALSE]
    dat <- cbind(pat, dat)[ind, ]
    dat2 <- data.frame(cbind(pat, cbind(as.character(d1), as.character(d2)))[ind, ])
    
    if (length(ind) == 1){dat2 <- t(dat2)}
    
    dimnames(dat2)[[2]] <- names
    res <- (dat2[dat[, 2] > dat[, 3], ])

    if (is.vector(res) == TRUE){dim(res) <- c(1, 3)}

    if (dim(res)[1] == 0){cat("No violations!\n")}

    if (dim(res)[1] > 0){
        
        colnames(res) <- names
        if (identical(typ, "R")){return(res)}
    
        if (identical(typ, "tex")){
        xtab1 <- xtable(res, align = "lrrr", caption = paste("Observations with ", names[2], " > ", names[3], ".", sep = ""), label = lab)
        xtab2 <- print(xtab1, include.rownames = FALSE, floating = FALSE, hline.after = c(-1, 0), type = "latex", size = "footnotesize", 
            sanitize.text.function = function(x){x}, tabular.environment = "longtable")
        }
    }
}



