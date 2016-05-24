print.mpmcorrelogram <-
function(x,...){
nclas <- length(x$clases)
print(data.frame(class = 1:nclas, distance.range=x$clases, rM = x$rM,
                       p =x$pvalues, p.Bonferroni=x$pval.Bonferroni,...))
}

