`check.datasets` <- 
function(x, y) {
    factor.count <- 0
    factor.species <- NULL
    na.count <- 0
    na.species <- NULL
    neg.count <- 0
    neg.species <- NULL
    all.good <- TRUE
    for (i in 1:ncol(x)) {
        if(is.factor(x[,i])) {
            factor.count <- factor.count+1
            factor.species <- c(factor.species, colnames(x)[i])
        }
        if(any(is.na(x[,i]))) {
            na.count <- na.count+1
            na.species <- c(na.species, colnames(x)[i])
        }
        if(any(x[,i] < 0, na.rm=T)) {
            neg.count <- neg.count+1
            neg.species <- c(neg.species, colnames(x)[i])
        }
    }
    if (factor.count > 0) {
        all.good <- FALSE
        cat("Warning:", factor.count, "variable(s) of the community dataset ( out of a total of", ncol(x), ") are factors\n")
        print(factor.species)
    }
    if (na.count > 0) {
        all.good <- FALSE
        cat("Warning:", na.count, "variable(s) of the community dataset ( out of a total of", ncol(x), ") with missing values\n")
        print(na.species)
    }
    if (neg.count > 0) {
        all.good <- FALSE
        cat("Warning:", neg.count, "variable(s) of the community dataset ( out of a total of", ncol(x), ") with negative values\n")
        print(neg.species)
    }
    if(nrow(x)!=nrow(y)){
        all.good <- FALSE
        cat("Warning: community and environmental datasets have different numbers of rows\n")
    }else{
        if(any(rownames(x)!=rownames(y))){
            all.good <- FALSE
            cat("Warning: rownames for community and environmental datasets are different\n")
        }
    }
    if (all.good == TRUE) {cat("OK\n")}
}
