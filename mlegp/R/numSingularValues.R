`numSingularValues` <-
function(Y, cutoff = 99) {
        p = singularValueImportance(Y)
        np = 1
        while (p[np] < cutoff && np < dim(p)[1]) {
           np = np+1
        }
        return (np)
}

