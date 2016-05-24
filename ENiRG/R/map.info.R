map.info <-
function(map.name, format = "grass") {
    if (format == "r") {
        result <- gmeta()[c("GISDBASE", "LOCATION_NAME", "MAPSET", "proj4")]
        obj1 <- execGRASS("r.info", map = map.name, flags = c("g"), 
            intern = TRUE)[1:9]
        obj1 <- data.frame(matrix(unlist(strsplit(obj1, split = "=")), ncol = 2, 
            byrow = T))
        output <- data.frame(value = obj1[, 2])
        rownames(output) <- obj1[, 1]
        result[["extension"]] <- output
        obj1 <- execGRASS("r.univar", map = map.name, flags = c("g"), 
            intern = TRUE)
        obj1 <- data.frame(matrix(unlist(strsplit(obj1, split = "=")), ncol = 2, 
            byrow = T))    
        output <- data.frame(value = obj1[, 2])
        rownames(output) <- obj1[, 1]
        result[["statistics"]] <- output
        result$class <- "metadata"
        return(result)
    }
    if (format == "grass") {
        execGRASS("r.info", map = map.name)
    }
}
