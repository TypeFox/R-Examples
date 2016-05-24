# Verified 1.3.18
est.rm <-
function(collection, list){
        pr = collection$data
        catalogo = collection$catalog
        disp = 1 : length(pr)
        disp = disp[-list]
        pr.c.new = list()
        cat.new = list()
        k = 0
        for(j in disp){
                k = k + 1
                pr.c.new [[k]] = pr [[j]]
                cat.new [[k]] = catalogo [[j]]
        }
        col = list(data = pr.c.new, catalog = cat.new)
        class(col) <- "Catalog"
        return(col)
}
