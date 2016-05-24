kde_bel.builder <- function(labs, test, train,
                            options = list(coef = 0.90))
    {
        stopifnot(length(train) == length(labs)) ##Why are you giving me extra labels?
        cl <- levels(labs)
        nclass <- length(cl)
        ntest <- length(test)

        kde.bel <- function(P)
            {
                B <- matrix(0, nrow = nclass, ncol = ntest)

                eval.bel <- function(cl)
                    sm.density(P[train[labs == cl],],
                               eval.points = P[test,],
                               eval.grid = FALSE, display = 'none')$estimate

                for(i in 1:nclass)
                    B[i,] <- eval.bel(cl[i])

                B <- B / rep(colSums(B, na.rm = TRUE), each = nclass)
                B <- options$coef * B
                bel.v <- rbind(B, 1 - colSums(B, na.rm = TRUE))
                row.names(bel.v) <- c(cl, 'Inf')
                colnames(bel.v) <- test
                return(bel.v)
            }
        
        return(kde.bel)
    }

knn_bel.builder <- function(labs, test, train,
                            options = list(k = 3, p = FALSE,
                                dist.type = c('euclidean', 'absolute', 'mahal'),
                                out = c('var', 'cv'), coef = 0.90))
    {
        stopifnot(length(train) == length(labs)) ##Why are you giving me extra labels?
        cl <- levels(labs)
        nclass <- length(cl)
        ntest <- length(test)
        knn.bel <- function(P)
            {
                B <- matrix(0, nrow = nclass, ncol = ntest)
                nn1 <- get.NN(P, k = options$k, p = options$p, dist.type = options$dist.type, nn.type = 'which', test = test, train = train)
                nn1 <- rbind(test, nn1)

                eval.bel <- function(cl)
                    {
                        M <- sm.density(P[train[labs == cl],], eval.points = P[nn1,], eval.grid = FALSE, display = 'none')$estimate
                        M2 <- matrix(M, ncol = ntest)
                        r <- apply(M2, MARGIN = 1, FUN = function(x) (x - M2[1,])^2)
                        result <- switch(options$out,
                                         var = rowSums(r),
                                         cv  = sqrt(rowSums(r)) / M2[1,])
                        return(result)
                    }
                
                for(i in 1:nclass)
                    B[i,] <- eval.bel(cl[i])

                B <- B / rep(colSums(B, na.rm = TRUE), each = nclass)
                B <- options$coef * B
                bel <- rbind(B, 1 - colSums(B, na.rm = TRUE))
                row.names(bel) <- c(cl, 'Inf')
                colnames(bel) <- test
                return(bel)
            }
        return(knn.bel)
    }

jit_bel.builder <- function(labs, test, train,
                            options = list(k = 3, p = FALSE, s = 5,
                                dist.type = c('euclidean', 'absolute', 'mahal'),
                                out = c('var', 'cv'), coef = 0.90))
    {
        stopifnot(length(train) == length(labs)) ##Why are you giving me extra labels?
        cl <- levels(labs)
        nclass <- length(cl)
        ntest <- length(test)
        
        jit.bel <- function(P)
            {
                B <- matrix(0, nrow = nclass, ncol = ntest)

                fake_points <- function()
                    {
                        dst <- get.NN(P, k = options$k, p = options$p, dist.type = options$dist.type, nn.type = 'max', test = test, train = train)
                        theta <- runif(n = ntest*options$s, min = 0, max = 2*pi)
                        R <- runif(ntest*options$s, 0, dst)
                        X <- R * cos(theta)
                        Y <- R * sin(theta)
                        Z <- cbind(P[test,1] + X, P[test,2] + Y)
                        return(Z)
                    }

                eval.bel <- function(cl)
                    {
                        M <- sm.density(P[train[labs == cl],], eval.points = nn1, eval.grid = FALSE, display = 'none')$estimate
                        y <- (M[(ntest+1):length(M)] - M[1:ntest]) ^ 2
                        r <- matrix(y, nrow = ntest)
                        result <- switch(options$out,
                                         var = rowSums(r),
                                         cv  = sqrt(rowSums(r)) / M[1:ntest])
                        return(result)
                    }        
                
                nn1 <- rbind(P[test,], fake_points())                
                for(i in 1:nclass)
                    B[i,] <- eval.bel(cl[i])
                B <- B / rep(colSums(B, na.rm = TRUE), each = 2)
                B <- options$coef * B
                bel.v <- rbind(B, 1 - colSums(B, na.rm = TRUE))
                row.names(bel.v) <- c(cl, 'Inf')
                colnames(bel.v) <- test
                return(bel.v)
            }
        return(jit.bel)
    }
