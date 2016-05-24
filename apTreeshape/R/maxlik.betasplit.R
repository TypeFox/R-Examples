maxlik.betasplit<-function (phylo, up = 10, remove.outgroup = FALSE, confidence.interval = "none",
    conf.level = 0.95, size.bootstrap = 100)
{
    vrais.aldous.fin <- function(i, n, b) {
       aux<-beta(b + i + 1, b + n - i + 1)/beta(i + 1, n -
            i + 1)
	if(is.na(aux) | (aux==Inf))
		aux<-(i/n)^b*(1-i/n)^(b)
	return(aux)
    }
    bbalance <- function(phylo) {
        return(t(apply(balance(phylo), FUN = function(x) {
            c(x[1], x[1] + x[2])
        }, MARGIN = 1)))
    }
    renorm.aldous <- function(n, beta) {
        return(sum(sapply(1:(n - 1), FUN = vrais.aldous.fin,
            n = n, b = beta)))
    }
    vrais.aldous.renorm <- function(i, n, beta) {
        return(vrais.aldous.fin(i, n, beta)/renorm.aldous(n,
            beta))
    }
    logvrais.aldous.phylo <- function(b, phylo, remove.outgroup = TRUE) {
            if (class(phylo) == "treeshape")
                bal <- bbalance(as.phylo(phylo))
            if (class(phylo) == "phylo")
                bal <- bbalance(phylo)
            if (remove.outgroup) {
                if ((bal[1, 1] <= 2) || ((bal[1, 2] - bal[1,
                  1]) <= 2))
                  bal <- bal[-1, ]
            }
        return(sum(log(apply(bal, FUN = function(x, b) {
            return(vrais.aldous.renorm(x[1], x[2], b))
        }, b = b, MARGIN = 1))))
    }
     logvrais.aldous.bal <- function(b,bal) {
        return(sum(log(apply(bal, FUN = function(x, b) {
            return(vrais.aldous.renorm(x[1], x[2], b))
        }, b = b, MARGIN = 1))))
    }
    if (class(phylo) == "treeshape")
        bal <- bbalance(as.phylo(phylo))
    if (class(phylo) == "phylo")
        bal <- bbalance(phylo)
    if (class(phylo) != "phylo" && class(phylo) != "treeshape") {
        print("The phylogeny shall be of class phylo or treeshape")
        return
    }
    if (remove.outgroup) {
        if ((bal[1, 1] <= 2) || ((bal[1, 2] - bal[1, 1]) <= 2))
            bal <- bal[-1, ]
    }
    nb.tip <- max(bal[1, ])
    optim_lik_aldous <- function(phylo, remove.outgroup) {
        optimize(f = function(x) {
            logvrais.aldous.phylo(x, phylo, remove.outgroup)
        }, lower = -2, upper = up, maximum = TRUE)
    }
    optim_lik_aldous_bal <- function(bal) {
        optimize(f = function(x) {
            logvrais.aldous.bal(x,bal)
        }, lower = -2, upper = up, maximum = TRUE)
    }
    res <- optim_lik_aldous(phylo, remove.outgroup)
    if (confidence.interval == "bootstrap") {
    	nb_b<-dim(bal)[1]
    	fun_aux<-function()
    	{
    		optim_lik_aldous_bal(bal[sample(1:nb_b,size=nb_b,replace=TRUE),])$maximum
    		
    	}
        thebeta <- replicate(size.bootstrap,fun_aux())
        up.conf <- 1 - ((1 - conf.level)/2)
        low.conf <- (1 - conf.level)/2
        conf_interval <- quantile(thebeta, c(low.conf, up.conf))
    }
    if (confidence.interval == "profile") {
    if( (res$objective-1.92)-logvrais.aldous.phylo(-2, phylo, remove.outgroup) <0)
    	low.conf <-(-2)
    else
      low.conf <-uniroot(f = function(x) {
            (res$objective-1.92)-logvrais.aldous.phylo(x, phylo, remove.outgroup)
        }, lower = -2, upper = res$maximum)$root
    if( (res$objective-1.92)-logvrais.aldous.phylo(up, phylo, remove.outgroup) <0)
    	up.conf <-up
    else
      up.conf <-uniroot(f = function(x) {
            (res$objective-1.92)-logvrais.aldous.phylo(x, phylo, remove.outgroup)
        }, lower = res$maximum, upper = up)$root
       conf_interval <-c(low.conf,up.conf)
    }
    if (confidence.interval == "none") {
        conf_interval <- NULL
    }
    return(list(max_lik = res$maximum, conf_interval = conf_interval))
}