#'Generating cluster data from the Chinese Restaurant Process
#'
#'@param n number of observations
#'
#'@param alpha concentration parameter
#'
#'@param hyperG0 base distribution hyperparameter
#'
#'@param verbose logical flag indicating whether info is written in the console.
#'
#'@importFrom stats runif rgamma rmultinom rnorm
#'
#'@export rCRP
#'
#'@examples
#'
#'rm(list=ls())
#'
#' d=2
#' hyperG0 <- list()
#' hyperG0[["NNiW"]] <- list()
#' hyperG0[["NNiW"]][["b_xi"]] <- rep(0,d)
#' hyperG0[["NNiW"]][["b_psi"]] <- rep(0,d)
#' hyperG0[["NNiW"]][["D_xi"]] <- 100
#' hyperG0[["NNiW"]][["D_psi"]] <- 8
#' hyperG0[["NNiW"]][["nu"]] <- d+1
#' hyperG0[["NNiW"]][["lambda"]] <- diag(c(1,1))
#'
#' hyperG0[["scale"]] <- list()
#'
#' set.seed(4321)
#' N <- 200
#' alph <- runif(n=1,0.2,2)
#' GvHD_sims <- rCRP(n=2*N, alpha=alph, hyperG0=hyperG0)
#' library(ggplot2)
#' q <- (ggplot(data=cbind.data.frame("D1"=GvHD_sims$data[1,],
#'                                   "D2"=GvHD_sims$data[2,],
#'                                   "Cluster"=GvHD_sims$cluster),
#'              aes(x=D1, y=D2))
#'       + geom_point(aes(colour=Cluster), alpha=0.6)
#'       + theme_bw()
#'       )
#' q
#'#q + stat_density2d(alpha=0.15, geom="polygon")
#'
#'\dontrun{
#' MCMCy1 <- DPMGibbsSkewT(z=GvHD_sims$data[,1:N],
#'                         hyperG0$NNiW, a=0.0001, b=0.0001, N=5000,
#'                         doPlot=TRUE, nbclust_init=64, plotevery=500,
#'                         gg.add=list(theme_bw()), diagVar=FALSE)
#'  s1 <- summary(MCMCy1, burnin=4000, thin=5,
#'                posterior_approx=TRUE)
#'  F1 <- FmeasureC(ref=GvHD_sims$cluster[1:N], pred=s1$point_estim$c_est)
#'
#'  # s <- summary(MCMCy1, burnin=4000, thin=5,
#'  #               posterior_approx=TRUE, K=1)
#'  # s2 <- summary(MCMCy1, burnin=4000, thin=5,
#'  #               posterior_approx=TRUE, K=2)
#'  # MCMCy2_seqPost<- DPMGibbsSkewT(z=GvHD_sims$data[,(N+1):(2*N)],
#'  #                                  hyperG0=s1$param_post$parameters,
#'  #                                  a=s1$param_post$alpha_param$shape,
#'  #                                  b=s1$param_post$alpha_param$rate,
#'  #                                  N=5000, doPlot=TRUE, nbclust_init=64, plotevery=500,
#'  #                                  gg.add=list(theme_bw()), diagVar=FALSE)
#'
#'  MCMCy2_seqPost <- DPMGibbsSkewT_SeqPrior(z=GvHD_sims$data[,(N+1):(2*N)],
#'                                            prior=s1$param_post, hyperG0=hyperG0$NNiW, , N=1000,
#'                                            doPlot=TRUE, nbclust_init=10, plotevery=100,
#'                                            gg.add=list(theme_bw()), diagVar=FALSE)
#'  s2_seqPost <- summary(MCMCy2_seqPost, burnin=600, thin=2)
#'  F2_seqPost <- FmeasureC(ref=GvHD_sims$cluster[(N+1):(2*N)], pred=s2_seqPost$point_estim$c_est)
#'
#'  MCMCy2 <- DPMGibbsSkewT(z=GvHD_sims$data[,(N+1):(2*N)],
#'                          hyperG0$NNiW, a=0.0001, b=0.0001, N=5000,
#'                          doPlot=TRUE, nbclust_init=64, plotevery=500,
#'                          gg.add=list(theme_bw()), diagVar=FALSE)
#'  s2 <- summary(MCMCy2, burnin=4000, thin=5)
#'  F2 <- FmeasureC(ref=GvHD_sims$cluster[(N+1):(2*N)], pred=s2$point_estim$c_est)
#'
#'  MCMCtot <- DPMGibbsSkewT(z=GvHD_sims$data,
#'                           hyperG0$NNiW, a=0.0001, b=0.0001, N=5000,
#'                           doPlot=TRUE, nbclust_init=10, plotevery=500,
#'                           gg.add=list(theme_bw()), diagVar=FALSE)
#'  stot <- summary(MCMCtot, burnin=4000, thin=5)
#'  F2tot <- FmeasureC(ref=GvHD_sims$cluster[(N+1):(2*N)], pred=stot$point_estim$c_est[(N+1):(2*N)])
#'
#'  c(F1, F2, F2_seqPost, F2tot)
#'}
#'
rCRP <- function(n=1000, alpha=2, hyperG0, verbose=TRUE){

    d <- length(hyperG0$NNiW[[1]])
    theta <- list()
    cluster <- numeric(n)
    z <- matrix(NA, ncol=n, nrow=d)

    for (c in 1:n){
        p0 <- alpha/(c-1+alpha)
        u <- runif(n=1, min = 0, max = 1)

        if (u<p0){
            # Accept: sample new value
            # cat("acceptation:", u, "<", p0, "\n")
            cluster[c] <- max(cluster)+1
            theta[[cluster[c]]] <- rNNiW(hyperG0$NNiW, diagVar=FALSE)
            theta[[cluster[c]]][["nu"]] <- 1 + rgamma(n=1, shape = 2, rate=1/10)

        }else{
            # Reject: sample old value
            # cat("rejection:", u, ">=", p0, "\n")
            u1  <-  u - p0
            weights <- summary(factor(cluster))[-1]
            cluster[c] <- which(rmultinom(n=1, size=1, prob=weights)==1)
        }
        w <- rgamma(n=1, shape=theta[[cluster[c]]]$nu/2, rate=theta[[cluster[c]]]$nu/2)
        ltnz <- rtruncnorm(n=1, a=0, sd=1/sqrt(w))
        eps <- matrix(rnorm(d), ncol=d)%*%chol(theta[[cluster[c]]]$S/w)
        z[1,c] <- theta[[cluster[c]]]$xi[1]+theta[[cluster[c]]]$psi[1]*ltnz+eps[,1]
        z[2,c] <- theta[[cluster[c]]]$xi[2]+theta[[cluster[c]]]$psi[2]*ltnz+eps[,2]

        if(verbose){
            cat(c,"/", n," sim\n", sep="")
        }
    }

    return(list("theta"=theta, "cluster"=as.factor(cluster), "data"=z))
}