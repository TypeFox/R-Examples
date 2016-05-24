spurr.ker <-
function(imarks, jmarks, dists, dranks, par = list(type=1, smark=1)) {
# Spurr's competition kernel, from Burkhart & Tome (2012) Sec. 9.2.2.1
    with(as.list(par),
         if(type == 1) x <- 0.5
         else if(type ==2) x <- -0.5
         else stop('spurr.ker: type must be 1 or 2')
         (2500  / length(dranks)) * (dranks - x) * (jmarks[smark] / dists)
    )
}
