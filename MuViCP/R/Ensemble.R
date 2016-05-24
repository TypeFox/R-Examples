##See file ensemble.org for documentation
ensemble <- function(dat, train, test, labs, bel.type = c('kde', 'knn', 'jit'), bel_options)
    {
        cols <- ncol(dat)
        bpamats <- NULL
        nmats <- 0
        bel_fn <- switch(bel.type,
                         kde = kde_bel.builder(labs, test, train, options = bel_options),
                         knn = knn_bel.builder(labs, test, train, options = bel_options),
                         jit = jit_bel.builder(labs, test, train, options = bel_options))

        try.matrices <- function(n, mc.cores = detectCores() - 2)
            {
                tmp <- vector('list' , n)
                
                for(i in 1:n)
                    {
                        ## seed <- .Random.seed
                        F <- basis_random(cols)
                        tmp[[i]] <- bpamat(info = F)
                    }

                l <- mclapply(tmp, function(x) bel_fn(P = dat %*% x$get.info()) , mc.cores = mc.cores)

                for(i in 1:n)
                    tmp[[i]]$assign.mat(l[[i]])
                
                if(nmats == 0)
                    bpamats <<- tmp
                else
                    bpamats <<- c(bpamats, tmp)

                nmats <<- nmats + n
            }

        get.class.mat <- function()
            {
                aa <- do.call('rbind', lapply(bpamats, function(x) t(x$get.mat())))
                aa[order(row.names(aa)),] ##This works
            }

        return(list(get.bpamats = function() bpamats,
                    try.matrices = try.matrices,
                    get.class.mat = get.class.mat))
    }
