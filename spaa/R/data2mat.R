data2mat <-
function(data = data)
{ 
    if (!any(colnames(data) == "abundance"))
        stop("A column named \"abundance\" must be speciefied.")
    if (!any(is.integer(data$abundance)))
        stop("Number of individuals must be integer!")
    col <- which(colnames(data) == "abundance")
    data1 <- data[,-col]
    abundance <- as.numeric(data[,col])
    result1 <- data.frame(rep(NA, sum(abundance)))
    colnames(result1) <- "plots"
    for (i in 1:(ncol(data)-1)){
        result1[, i] <-  rep(as.character(data[, i]), abundance)
    }
    result <- table(result1)
    return(result)
}

