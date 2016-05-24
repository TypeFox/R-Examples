genPattern<-function (th, it, model = NULL, D = 1, seed = NULL) 
{
    it <- rbind(it)
    if (!is.null(seed)) 
        set.seed(seed)
res<-matrix(NA,length(th),nrow(it))
for (iter in 1:length(th)){
    if (is.null(model))
        res[iter,] <- rbinom(nrow(it), 1, Pi(th[iter], it, model = model, 
            D = D)$Pi)
    else {
        pr <- Pi(th[iter], it, model = model, D = D)$Pi
        RES <- NULL
        for (i in 1:nrow(pr)) {
            pp <- pr[i, ][!is.na(pr[i, ])]
            vec <- rmultinom(n = 1, size = 1, prob = pp)
            RES[i] <- (1:nrow(vec))[vec[, 1] == 1] - 1
        }
      res[iter,]<-RES
    }
}
if (length(th)==1) res<-as.numeric(res)
    return(res)
}
