# Verified 1.3.18
# Version 3.0
est.union = function(collection, fun = mean, return.matrix = FALSE){
        estado.pr = collection$data[[1]]
        for ( i in 2:length(collection$data) ) {
                estado.pr = ts.union(estado.pr, collection$data[[i]])
        }
        colnames(estado.pr) <- paste0("C", 1:length(collection$data))
        if(return.matrix){ return(estado.pr) }
        st = start(estado.pr)
        estado.pr = apply(estado.pr, 1, fun, na.rm = TRUE)
        estado.pr = ts(estado.pr, start=st, frequency=12)
        col = collection
        col$union = estado.pr
        class(col) <- "Catalog"
        return(col)
}
