write.dl <-
function (x, file = NULL, type = c("nodelist", "fullmat")) 
{
    if (isTRUE(is.na(dim(x)[3]) == TRUE) == FALSE) {
        if (isTRUE(is.null(dimnames(x)[[3]])) == TRUE) 
            dimnames(x)[[3]] <- 1:dim(x)[3]
    }
    if (isTRUE(is.null(dimnames(x)[[1]])) == TRUE) 
        dimnames(x)[[1]] <- dimnames(x)[[2]] <- 1:dim(x)[1]
    switch(match.arg(type), fullmat = {
        cat(paste("DL", collapse = "\n"), file = file, sep = "\n", 
            append = TRUE)
        cat(paste("N=", dim(x)[1], sep = ""), file = file, append = TRUE)
        cat(paste(" NM=", dim(x)[3], sep = ""), file = file, 
            append = TRUE)
        cat("", file = file, sep = "\n", append = TRUE)
        cat(paste("FORMAT = FULLMATRIX DIAGONAL PRESENT", collapse = "\n"), 
            file = file, sep = "\n", append = TRUE)
        cat(paste("ROW LABELS:", collapse = "\n"), file = file, 
            sep = "\n", append = TRUE)
        for (i in 1:dim(x)[1]) {
            cat(paste("\"", dimnames(x)[[1]][i], "\"", sep = "", 
                collapse = "\n"), file = file, sep = "\n", append = TRUE)
        }
        cat(paste("COLUMN LABELS:", collapse = "\n"), file = file, 
            sep = "\n", append = TRUE)
        for (i in 1:dim(x)[2]) {
            cat(paste("\"", dimnames(x)[[2]][i], "\"", sep = "", 
                collapse = "\n"), file = file, sep = "\n", append = TRUE)
        }
        cat(paste("LEVEL LABELS:", collapse = "\n"), file = file, 
            sep = "\n", append = TRUE)
        for (i in 1:dim(x)[3]) {
            cat(paste("\"", dimnames(x)[[3]][i], "\"", sep = "", 
                collapse = "\n"), file = file, sep = "\n", append = TRUE)
        }
        cat(paste("DATA:", collapse = "\n"), file = file, sep = "\n", 
            append = TRUE)
        for (i in 1:dim(x)[3]) {
            for (j in 1:dim(x)[1]) {
                cat(paste(x[j, , i]), file = file, collapse = "\n", 
                  append = TRUE)
            }
        }
        rm(i, j)
    }, nodelist = {
        tmp <- x
        dimnames(tmp)[[1]] <- dimnames(tmp)[[2]] <- 1:dim(x)[1]
        rs <- rel.sys(tmp)
        rm(tmp)
        vec <- vector()
        for (i in 1:length(rs$ties)) vec <- append(vec, rs$ties[[i]])
        rm(rs)
        vec <- sort(vec)
        lst <- list()
        for (i in 1:dim(x)[1]) lst[[i]] <- i
        for (j in 1:length(vec)) lst[[as.numeric(strsplit(vec[[j]][1], 
            ", ")[[1]][1])]] <- append(lst[[as.numeric(strsplit(vec[[j]][1], 
            ", ")[[1]][1])]], as.numeric(strsplit(vec[[j]][1], 
            ", ")[[1]][2]))
        rm(vec)
        cat(paste("dl", " n=", dim(x)[1], " format=nodelist1", 
            sep = "", collapse = "\n"), file = file, sep = "\n", 
            append = TRUE)
        cat(paste("labels: ", collapse = "\n"), file = file, 
            sep = "\n", append = TRUE)
        for (i in 1:dim(x)[1]) {
            cat(paste("\"", dimnames(x)[[1]][i], "\" ", sep = ""), 
                file = file, sep = " ", append = TRUE)
        }
        cat("", file = file, sep = "\n", append = TRUE)
        cat(paste("data:", collapse = "\n"), file = file, sep = "\n", 
            append = TRUE)
        for (i in 1:dim(x)[1]) {
            cat(paste(lst[[i]]), file = file, collapse = "\n", 
                append = TRUE)
        }
        rm(i, j)
    })
}
