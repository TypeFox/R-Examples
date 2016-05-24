# Verified 1.3.18
est.sort <-
function(collection, ascending = T, by.year.only = F) {
        pr = collection$data
        catalogo = collection$catalog
        t = plyr::ldply(pr, function(x) { c(start(x), end(x)) }  )
        if ( ascending == T ) {
                if ( by.year.only == T ) {
                        x = sort(t$V1, index.return=T)
                } else {
                        x = sort(t$V1+t$V2/12, index.return=T)
                }
        } else {
                if ( by.year.only == T ) {
                        x = sort(t$V3, index.return=T)
                } else {
                        x = sort(t$V3+t$V4/12, index.return=T)
                }
        }
        t.o = t
        t.o[,] = t[x$ix,]
        p.o = pr
        c.o = catalogo
        for (i in 1 : length(catalogo)) { 
                c.o[[i]] = catalogo[[x$ix[i]]]
                p.o[[i]] = pr[[x$ix[i]]]
        }
        col = list(catalogo = c.o, data = p.o)
        class(col) <- "Catalog"
        return( col )
}
