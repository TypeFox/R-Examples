write.hclust <- function(x, file, prefix, h = NULL, k = NULL, append = FALSE,
                         dec = ","){
    grp <- cutree(x, h=h, k=k)
    grp.lab <- formatC(grp, digits = 0, width = max(nchar(as.character(grp))), 
                       format = "f", flag = "0")
    res <- data.frame("Sample" = names(grp),
                      "Cluster" = grp, 
                      "Cluster ID" = paste(prefix, "_H", h, "_", grp.lab, sep = ""),
                      check.names = FALSE)
    if(append){
        write.table(res, file = file, append = append, dec = dec, 
                    sep ="\t", row.names = FALSE, col.names = FALSE)
    }else{
        write.table(res, file = file, append = append, dec = dec, 
                    sep ="\t", row.names = FALSE)
    }
}
