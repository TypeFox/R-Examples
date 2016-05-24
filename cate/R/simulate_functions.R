#' Generate simulation data set
#'
#' @description \code{gen.sim.data} generates data from the following model 
#' Y = X_0 Beta_0^T + X_1 Beta_1^T + Z Gamma^T + E Sigma^{1/2},
#' Z|X_0, X_1 = X_0 Alpha_0^T + X_1 Alpha_1^T + D,
#' cov(X_0, X_1) ~ Sigma_X
#'
#' @param n number of observations
#' @param p number of observed variables
#' @param r number of confounders
#' @param d0 number of nuisance regression covariates
#' @param d1 number of primary regression covariates
#' @param X.dist the distribution of X, either "binary" or "normal"
#' @param alpha association of X and Z, a r*d vector (d = d0 + d1)
#' @param beta treatment effects, a p*d vector
#' @param beta.strength strength of beta
#' @param beta.nonzero.frac if beta is not specified, fraction of nonzeros in beta
#' @param Gamma confounding effects, a p*r matrix
#' @param Gamma.strength strength of Gamma, more precisely the mean of square entries of Gamma * alpha
#' @param Gamma.beta.cor the "correlation" (proportion of variance explained) of beta and Gamma
#' @param Sigma noise variance, a p*p matrix or p*1 vector or a single real number
#' @param seed random seed
#'
#' @return a list of objects
#' \describe{
#' \item{X0}{matrix of nuisance covariates}
#' \item{X1}{matrix of primary covariates}
#' \item{Y}{matrix Y}
#' \item{Z}{matrix of confounders}
#' \item{alpha}{regression coefficients between X and Z}
#' \item{beta}{regression coefficients between X and Y}
#' \item{Gamma}{coefficients between Z and Y}
#' \item{Sigma}{noise variance}
#' \item{beta.nonzero.pos}{the nonzero positions in beta}
#' \item{r}{number of confounders}
#' }
#' 
#' @import stats
#' 
#' @export
#'
gen.sim.data <- function(n,
                         p,
                         r,
                         d0 = 0,
                         d1 = 1,
                         X.dist = c("binary", "normal"),
                         alpha = matrix(0.5, r, d0 + d1),
                         beta = NULL,
                         beta.strength = 1,
                         beta.nonzero.frac = 0.05,
                         Gamma = NULL,
                         Gamma.strength = sqrt(p),
                         Gamma.beta.cor = 0,
                         Sigma = 1,
                         seed = NULL) {

    ## check arguments
    X.dist <- match.arg(X.dist, c("binary", "normal"))

    ## set a random seed, if necessary
    if (!is.null(seed)) {
        set.seed(seed)
    }

    d <- d0 + d1

### comment: we don't need first d0 columns of beta to be sparse

    ## generate a beta if it is NULL
    if (is.null(beta)) {
        beta <- matrix(0, p, d)
        beta.nonzero.pos <- sample(1:p, p*beta.nonzero.frac)
### comment:different columns can have different non-zero positions
        beta[beta.nonzero.pos, ] <- beta.strength
        ## an alternative strategy:
        ## beta[beta.nonzero.pos, ] <- beta.strength * rnorm(p*beta.nonzero.frac)
    }

    ## generate a Gamma if it is NULL
    if (is.null(Gamma)) {
        Gamma <- matrix(rnorm(p*r), p, r)
### add the next line if you want to make Gamma very correlated with beta
        if(sum(beta^2) != 0) {
            Gamma[,1] <- beta * Gamma.beta.cor / sqrt(sum(beta^2)) * sqrt(p) + sqrt(1 - Gamma.beta.cor^2) * rnorm(p)
        }
        Gamma <- t(t(qr.Q(qr(Gamma))) * sample(c(-1, 1), r, replace = T))
        Gamma <- t(t(Gamma) * Gamma.strength)
    }

    ## make sure Sigma is a p*p matrix
    SS <- Sigma
    if (length(Sigma) == p || length(Sigma) == 1) {
        #Sigma <- diag(Sigma)
		if (length(Sigma) == 1)
			Sigma <- rep(Sigma, p)
    } else if (length(Sigma) != p*p) {
        stop("The noise variance Sigma can either be a p*p matrix, a p*1 vector, or a single real number.")
    }


    ## simulate X
    X <- switch(X.dist,
                binary = matrix(2 * rbinom(n*d, 1, 0.5) - 1, n, d),
                normal = matrix(rnorm(n*d), n, d))
                                        #X <- scale(X)

    ## simulate Z
    Z <- X %*% t(alpha) + matrix(rnorm(n*r), n, r)

    ## simulate Y
	if (length(Sigma) == p) {
		Y <- X %*% t(beta) + Z %*% t(Gamma) + t(t(matrix(rnorm(n*p), n, p)) * sqrt(Sigma))
	} else {
		 Sigma <- matrix(Sigma, p, p)
    	Sigma.svd <- svd(Sigma)
    	Sigma.half <- Sigma.svd$u %*% diag(sqrt(Sigma.svd$d)) %*% t(Sigma.svd$v)
		Y <- X %*% t(beta) + Z %*% t(Gamma) + matrix(rnorm(n*p), n, p) %*% Sigma.half
	}


    if (d0 >= 1) {
        X0 = X[, 1:d0, drop = FALSE]
    } else {
        X0 = NULL
    }

    return(list(X0 = X0,
                X1 = X[, (d0+1):d, drop = FALSE],
                Y = Y,
                Z = Z,
                alpha = alpha,
                beta = beta,
                Gamma = Gamma,
                Sigma = SS,
                beta.nonzero.pos = beta.nonzero.pos,
                beta.zero.pos = setdiff(1:(p), beta.nonzero.pos),
                r = r))

}

# 
# #' Simulation
# #'
# cate.simulation.example <- function() {
# 
#     power <- function(p.value, sig.level = 0.05) {
#         sum(p.value < sig.level) / length(p.value)
#     }
# 
#     library(corpcor)
# 
#     settings <- expand.grid(n = c(100),
#                             p = c(1000),
#                             r = c(2, 5),
#                             alpha.strength = c(5),
#                             beta.zero.frac = c(0.02),
#                             Gamma.beta.cor = c(0.5))
#     methods <- expand.grid(adj.method = c("rr", "nc", "nc.no.correction"),
#                            fa.method = c("ml", "pc", "esa"),
#                            calibration = c(FALSE, TRUE))
#     methods <- rbind(methods,
#                      expand.grid(adj.method = "naive",
#                                  fa.method = "ml",
#                                  calibration = c(FALSE, TRUE)))
#     methods <- rbind(methods,
#                      expand.grid(adj.method = c("RUV4", "RUVinv", "leapp", "sva"),
#                                  fa.method = "ml",
#                                  calibration = TRUE))
# 
#     niter <- 10
# 
#     df <- data.frame(matrix(0, nrow(settings) * (nrow(methods) + 1) * niter,
#                             ncol(settings) + ncol(methods) + 3))
#     names(df) <- c(names(settings), names(methods), "type.I.error", "power", "fdp")
# 
#     i <- 0
#     t0 <- proc.time()[3]
#     for (setting.ind in 1:nrow(settings)) {
# 
#         setting <- settings[setting.ind,]
# 		print(setting)
# 
#         for (iter in 1:niter) {
# 
#             print(iter)
# 
#             i <- i + 1
# 
#             n <- setting$n
#             p <- setting$p
#             r <- setting$r
#             data <- gen.sim.data(n = n, p = p, r = r, d0 = 0, d1 = 1,
#                                  alpha = matrix(sqrt(setting$alpha.strength/r), r, 1),
#                                  Gamma.strength = c(seq(3, 1, length = r)),
#                                  Gamma.beta.cor = setting$Gamma.beta.cor,
#                                  beta.nonzero.frac = setting$beta.zero.frac,
#                                  beta.strength = 3 * sqrt(1 + setting$alpha.strength) / sqrt(n),
#                                  #Sigma = rexp(p),
# 								 #Sigma = 1,
# 								 Sigma = 1 / rgamma(p, 3, rate = 2),
#                                  seed = (iter *  274876858367) %% 1000000)
#             data$X0 <- cbind(rep(1, n), data$X0)
#             nc <- sample(data$beta.zero.pos, 30)
# 
#             #save(data, file = "puzzle.rda")
# 
#                                         # oracle
#             output <- lm(data$Y ~ data$X1 + data$Z - 1)
# 
#             s <- summary(output)
#             t <- sapply(1:ncol(data$Y), function(j) {s[[j]]$coefficients[1,3]})
#             p.value <- 2 * (1 - pnorm(abs(t)))
#             type.I.error <- power(p.value[data$beta.zero.pos])
#             pow <- power(p.value[data$beta.nonzero.pos])
#             discoveries <- which(p.adjust(p.value, "BH") < 0.2)
#             fdp <- length(setdiff(discoveries, data$beta.nonzero.pos)) / max(length(discoveries), 1)
# 
#             df[i, ] <- c(setting, 0, 0, 0, type.I.error, pow, fdp)
# 
#             for (method.ind in 1:nrow(methods)) {
#                 method <- methods[method.ind, ]
# 
#                 i <- i + 1
# 
#                 if ((i %% 100) == 0) {
#                     print(paste("Time left:",
#                                 (proc.time()[3] - t0) / i * (nrow(df) - i)))
#                 }
# 
#                 adj.method <- method$adj.method
#                 if (as.character(method$adj.method) == "nc.no.correction") {
#                     adj.method <- "nc"
#                     correction <- FALSE
#                 } else {
#                     correction <- TRUE
#                 }
#                 if (length(grep("ruv", as.character(adj.method), ignore.case = TRUE)) > 0) {
#                     output <- ruv.wrapper(data$X1, data$Y, data$X0,
#                                           r = r,
#                                           nc = nc,
#                                           ruv.method = as.character(adj.method))
#                 } else if (adj.method == "leapp") {
#                     output <- leapp.wrapper(data$X1, data$Y, data$X0,
#                                           r = r)
#                 } else if (adj.method == "sva") {
#                     try(output <- sva.wrapper(data$X1, data$Y, data$X0,
#                                           r = r))
#                 } else {
#                     output <- cate(data$X1, data$Y, data$X0,
#                                    fa.method = as.character(method$fa.method),
#                                    adj.method = as.character(adj.method),
#                                    r = r,
#                                    psi = psi.bisquare,
#                                    nc = nc,
#                                    nc.var.correction = correction,
#                                    calibrate = method$calibration)
#                 }
#                 if (adj.method == "nc") {
#                     type.I.error <- power(output$beta.p.value[setdiff(data$beta.zero.pos, nc)])
#                 } else {
#                     type.I.error <- power(output$beta.p.value[data$beta.zero.pos])
#                 }
#                 pow <- power(output$beta.p.value[data$beta.nonzero.pos])
#                 discoveries <- which(p.adjust(output$beta.p.value, "BH") < 0.2)
#                 fdp <- length(setdiff(discoveries, data$beta.nonzero.pos)) / max(length(discoveries), 1)
# 
#                 df[i, ] <- c(setting, method, type.I.error, pow, fdp)
#             }
#         }
#     }
# 
#     library(reshape2)
#     dff <- melt(df, id = names(df)[1:9])
#     library(ggplot2)
#     dff$adj.method <- factor(dff$adj.method)
#     dff$fa.method <- factor(dff$fa.method)
#     levels(dff$adj.method) <- c("oracle", "rr", "nc", "nc2", "naive", "RUV4", "RUVinv", "leapp", "sva")
#     levels(dff$fa.method) <- c("", "ml", "pc", "esa")
#     dff$method <- factor(paste(dff$adj.method, dff$fa.method, dff$calibration, sep = "."))
#     cate.plot <- ggplot(subset(dff, alpha.strength == 1 & Gamma.beta.cor == 0 & beta.zero.frac == 0.05)) + aes(x = method) + geom_boxplot(aes(y = value)) + facet_grid(r~variable) + geom_hline(yintercept = 0.05, linetype = "dashed", col = "red") + geom_hline(yintercept = 0.2, linetype = "dashed", col = "red") + theme(axis.text.x = element_text(angle = 15))
#     cate.plot
#     ggsave("cate.pdf", cate.plot)
# 
# }
# 
# 
# cate.demo <- function(data.set = c("gender", "alzheimers", "COPD")) {
# 
#     require(MASS)
#     library(robustbase)
#     library(ggplot2)
#     library(GEOquery)
# 
#     if (data.set == "gender") {
#         library(ruv.data.gender.sm)
#         data(gender.sm)
#         data <- gender.sm
#     } else if (data.set == "alzheimers") {
#         library(ruv.data.alzheimers.sm)
#         data(alzheimers.sm)
#         data <- alzheimers.sm
#         data$Z <- matrix(1, nrow(data$Y), 1)
#     } else if(data.set == "COPD") {
#     	gse22148 <- getGEO("GSE22148")
#         covariates <- pData(phenoData(gse22148[[1]]))[, c(10:11, 13, 15:16)]
#         colnames(covariates) <- c("Age", "Gender", "Disease", "Batch", "Date")
#         covariates$Age <- sapply(covariates$Age, function(s) as.numeric(strsplit(as.character(s), " ")[[1]][2]))
#         for (i in 2:4) {
#             covariates[, i] <- as.factor(sapply(covariates[, i],
#                                                 function(s) (strsplit(as.character(s), " ")[[1]][2])))
#         }
#         covariates[, 5] <- as.factor(sapply(covariates[, 5], function(s) (strsplit(as.character(s), " ")[[1]][3])))
#         data <- list()
#         X <- cbind(Age = covariates$Age, lm(Age~ ., data = covariates, x = T)$x[, -c(1, 11)])
#         data$X <- X[, 3, drop = F]
#         Y <- t(exprs(gse22148[[1]]))
#         data$Y <- Y[, apply(Y, 2, function(x) sum(is.na(x)) == 0)]
#         data$Z <- cbind(rep(1, nrow(Y)), X[, -3])
#         data$pctl <- rep(FALSE, nrow(Y))
#     }
# 
#     if (data.set == "gender") {
#         output.naive <- cate(data$X, data$Y,
#                              matrix(1, nrow(data$Y), 1),
#                              fa.method = "pc",
#                              adj.method = "naive",
#                              r = 0,
#                              calibrate = FALSE)
# 
#         xlow <- switch(data.set,
#                        gender = -1,
#                        alzheimers = -4,
#                        COPD = -8)
#         xup <- switch(data.set,
#                       gender = 1,
#                       alzheimers = 4,
#                       COPD = 8)
#         xtext <- switch(data.set,
#                         gender = 0.5,
#                         alzheimers = 0,
#                         COPD = 5)
#         ytext <- switch(data.set,
#                         gender = 5,
#                         alzheimers = 0.6,
#                         COPD = 0.13)
# 
#         plot.naive <- ggplot() + aes(x = output.naive$beta.t, y = ..density..) + geom_histogram(binwidth = 0.05 * xup / 2, colour = "black", fill = "white") + xlim(c(xlow, xup)) + geom_line(aes(x = seq(xlow, xup, 0.01), y = dnorm(seq(xlow, xup, 0.01), median(output.naive$beta.t), mad(output.naive$beta.t))), col = "violetred3", size = 1) + geom_text(aes(x = xtext, y = ytext, label = paste0("N(", signif(median(output.naive$beta.t), 2),",", signif(mad(output.naive$beta.t), 2), "^2)")), col = "violetred3", show_guide= FALSE, size = 7) + xlab("t-statistics") + ylab("density") + theme_bw(base_size = 20)
#         plot.naive
#         ggsave(paste0("../../plots/", data.set, "-naive.pdf"), plot.naive)
#     }
# 
#     output.batch <- cate(data$X, data$Y,
#                          data$Z,
#                          fa.method = "ml",
#                          adj.method = "naive",
#                          r = 0,
#                          calibrate = FALSE)
# 
#     xlow <- switch(data.set,
#                    gender = -1,
#                    alzheimers = -4,
#                    COPD = -8)
#     xup <- switch(data.set,
#                   gender = 1,
#                   alzheimers = 4,
#                   COPD = 8)
#     xtext <- switch(data.set,
#                     gender = 0.6,
#                     alzheimers = 0,
#                     COPD = 5)
#     ytext <- switch(data.set,
#                     gender = 1.3,
#                     alzheimers = 0.6,
#                     COPD = 0.13)
# 
#     plot.batch <- ggplot() + aes(x = output.batch$beta.t, y = ..density..) + geom_histogram(binwidth = 0.05 * xup / 2, colour = "black", fill = "white") + xlim(c(xlow, xup)) + geom_line(aes(x = seq(xlow, xup, 0.01), y = dnorm(seq(xlow, xup, 0.01), median(output.batch$beta.t), mad(output.batch$beta.t))), col = "violetred3", size = 1) + geom_text(aes(x = xtext, y = ytext, label = paste0("N(", signif(median(output.batch$beta.t), 2),",", signif(mad(output.batch$beta.t), 2), "^2)")), col = "violetred3", show_guide= FALSE, size = 7) + xlab("t-statistics") + ylab("density") + theme_bw(base_size = 20)
#     plot.batch
#     ggsave(paste0("../../plots/", data.set, "-batch.pdf"), plot.batch)
# 
#                                         # confounder adjustment
#     r.max <- switch(data.set,
#                     gender = 40,
#                     alzheimers = 20,
#                     COPD = 60)
#     summary.statistics <- matrix(0, r.max + 1, 9)
# 
#     library(moments)
# 
#     for (r in 0:r.max) {
#         print(r)
#         if (r == 0) {
#             output <- output.batch
#         } else {
#             output <- cate(data$X, data$Y,
#                            data$Z,
#                            fa.method = "ml",
#                            adj.method = "rr",
#                            r = r,
#                            psi = psi.bisquare,
#                            nc = data$hkctl,
#                            nc.var.correction = TRUE,
#                            calibrate = FALSE)
#         }
# 
#         summary.statistics[r+1, ] <- c(
#              length(which( 1 - pnorm(abs(output$beta.t - median(output$beta.t)) / mad(output$beta.t)) < 0.005)),
#             length(intersect(which( 1 - pnorm(abs(output$beta.t - median(output$beta.t)) / mad(output$beta.t)) < 0.005), which(data$pctl))),
#             length(intersect(which(rank(1 - pnorm(abs(output$beta.t - median(output$beta.t)) / mad(output$beta.t))) <= 100), which(data$pctl))),
#             mean(output$beta.t), median(output$beta.t),
#             sd(output$beta.t), mad(output$beta.t),
#             skewness(output$beta.t), mc(output$beta.t))
# 
#     }
# 
#     summary.statistics <- cbind(0:r.max, summary.statistics)
#     colnames(summary.statistics) <- c("r", "#sig.", "X/Y", "top 100", "mean", "median", "sd", "mad", "skewness", "medcouple")
# 
#     library(xtable)
#     rows <- switch(data.set,
#                    gender = c(1:6,8,11,16,21,26,31,41),
#                    alzheimers = c(1:6,8,10,12,14,16,21),
#                    COPD = c(1:6,8,11,16,21,31,41,50,61))
#     cols <- c(1,5:10,2:4)
#     print(xtable(summary.statistics[rows, cols],
#                  display = c(rep("d", 2), rep("fg", 6), rep("d", 3)),
#                  digits = c(0,0,2,2,3,3,3,2,0,0,0),
#                  caption = data.set), include.rownames = FALSE)
# 
#     r <- switch(data.set,
#                 gender = 25,
#                 alzheimers = 11,
#                 COPD = 49)
#     output <- cate(data$X, data$Y,
#                    data$Z,
#                    fa.method = "ml",
#                    adj.method = "rr",
#                    r = r,
#                    psi = psi.bisquare,
#                    nc = data$hkctl,
#                    nc.var.correction = TRUE,
#                    calibrate = FALSE)
#     plot.output <- ggplot() + aes(x = as.vector(output$beta.t), y = ..density..) + geom_histogram(binwidth = 0.2, colour = "black", fill = "white") + xlim(c(-6, 6)) + ylim(c(0, 0.35)) + geom_line(aes(x = seq(-6, 6, 0.01), y = dnorm(seq(-6, 6, 0.01), median(output$beta.t), mad(output$beta.t))), col = "violetred3", size = 1) + geom_text(aes(x = 3, y = 0.3, label = paste0("N(", signif(median(output$beta.t), 2),",", signif(mad(output$beta.t), 3), "^2)")), col = "violetred3", show_guide= FALSE, size = 7) + xlab("statistics") + ylab("density") + theme_bw(base_size = 20)
#     plot.output
#     ggsave(paste0("../../plots/", data.set, "-output.pdf"), plot.output)
# 
# }
# 
# 
# bcv.plot <- function() {
#     load("bcvResult.rda")
#     df <- data.frame(r = as.numeric(names(results$bcvErr.gender)),
#                      error = results$bcvErr.gender/results$bcvErr.gender[1])
#     plot.gender <- ggplot(df) + aes(x = r, y = error) + geom_point(size = 3, shape = 1) + geom_line() + theme_bw(base_size = 18) + scale_y_log10(breaks = c(0.1, 0.2, 0.5, 1)) + ylab("relative BCV MSE") + xlab("r") + geom_point(aes(x = r[which.min(error)], y = min(error)), col = "violetred3", size = 3)
#     ggsave("../../plots/gender-bcv.pdf", plot.gender)
# 
#     df <- data.frame(r = as.numeric(names(results$bcvErr.alzheimer)),
#                      error = results$bcvErr.alzheimer/results$bcvErr.alzheimer[1])
#     plot.alzheimers <- ggplot(df) + aes(x = r, y = error) + geom_point(size = 3, shape = 1) + geom_line() + theme_bw(base_size = 18) + scale_y_log10(breaks = c(0.1, 0.2, 0.5, 1)) + ylab("relative BCV MSE") + xlab("r") + geom_point(aes(x = r[which.min(error)], y = min(error)), col = "violetred3", size = 3)
#     plot.alzheimers
#     ggsave("../../plots/alzheimers-bcv.pdf", plot.alzheimers)
# 
#     df <- data.frame(r = as.numeric(names(results$bcvErr.COPD)),
#                      error = results$bcvErr.COPD/results$bcvErr.COPD[1])
#     plot.COPD <- ggplot(df) + aes(x = r, y = error) + geom_point(size = 3, shape = 1) + geom_line() + theme_bw(base_size = 18) + scale_y_log10(breaks = c(0.4, 0.5, 0.8, 1)) + ylab("relative BCV MSE") + xlab("r") + geom_point(aes(x = r[which.min(error)], y = min(error)), col = "violetred3", size = 3)
#     plot.COPD
#     ggsave("../../plots/COPD-bcv.pdf", plot.COPD)
# 
# }
