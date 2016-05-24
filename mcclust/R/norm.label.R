`norm.label` <-
function(cl){
    old.label <- names(table(cl))
    apply(matrix(cl), 1, function(x) which(old.label==as.character(x)))
    }

