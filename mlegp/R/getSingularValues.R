`getSingularValues` <-
function(Y) {
        Y = t(Y)
        s = svd(Y)
        p = matrix(s$d)
        return (p)
}

