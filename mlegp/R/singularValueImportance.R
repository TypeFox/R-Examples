`singularValueImportance` <-
function(Y) {
        s = svd(Y)
        p = matrix(s$d)
        total = sum(p**2)
        p[1] = p[1]**2/total * 100
        for (i in 2:dim(p)[1]) {
                p[i] = p[i-1] + p[i]**2 / total * 100
        }
        return (p)
}

