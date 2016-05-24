`Table.N.Per` <-
function(var, dep, subset = !is.na(var))
    {
        var <- as.factor(var)
        dep <- as.factor(dep)
        ta <- table(var[subset], dep[subset])
        dimnames(ta) <- list(levels(var), levels(dep))
        per <- matrix(nrow = dim(ta)[1], ncol = dim(ta)[2])
        dimnames(per) <- list(levels(var[subset]), c("%", "%"))
        for(i in 1.:dim(ta)[1.]) {
                for(j in 1.:dim(ta)[2]) {
                        per[i, j] <- cbind(as.matrix(round((ta[i, j] * 100)/sum(ta[, j]), 1)))
                }
        }
        tp <- cbind(ta, per)
        tp <- tp[, order(rep(1:2, 2))]
        list(tp = tp)
    }

