cuth <- function(heig, ind, models, m.full, cont){
    y <- m.full$y
    sp <- list()
    form <- c()
    d <- m.full$model
    if (heig[ind] == max(heig)){
        form <- c(form, 'no cont')
    } else {
        if (length(models) > 0){
            sp <- lapply(models, function(x) cutree(x, h = heig[ind]))
            names(sp) <- colnames(d)[2:(length(models)+1)]
        } else {
            sp <- list()
        }
        if (cont > 0){
            form <- c(form, names(which((heig > heig[ind]) & (names(heig) != 'fac') )))
        }
    }
    return(list(n_cont = form,  Crit = 0, SPart = sp))
}
