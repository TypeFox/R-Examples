#######################################################################
# rEMM - Extensible Markov Model (EMM) for Data Stream Clustering in R
# Copyrigth (C) 2011 Michael Hahsler
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

synthetic_stream <- function(k=10, d=2, n_subseq=100, p_transition=.5, p_swap=0,
        n_train=5000, n_test=1000, p_outlier = .01, rangeVar=c(0,0.005)) {

    test <- NA
    train <- NA
    outlier_position <- NA
    sequence_train <- NA
    sequence_test <- NA
    swap_train <- NA
    swap_test <- NA

    get_mu <- function(k) runif(k, min=.1, max=.9)
    get_outlier <- function(k) runif(k, min=0, max=1)

    mu <- replicate(d, get_mu(k))

    Sigma <- replicate(k, genPositiveDefMat("unifcorrmat", 
                    rangeVar=rangeVar, dim=d)$Sigma, 
            simplify=FALSE)

    #subseq <- sample(1:k, replace=TRUE)
    subseq <- integer(n_subseq)
    subseq[1] <- 1L
    for(i in 2:n_subseq) {
        if(runif(1) < p_transition) subseq[i] <- sample(1:k, 1)
        else subseq[i] <- subseq[i-1]
        }

    model <- list(k=k, d=d, mu=mu, Sigma=Sigma, subseq=subseq)

    if(n_train>0) {
        sequence_train <- rep(subseq, n_train/n_subseq)


        ## randomly mess up sequence missing
        if(p_swap>0) {
            swap_train <- which(runif(n_train-1)<p_swap)+1L
            for(i in swap_train) {
                sequence_train[(i-1):i] <-
                rev(sequence_train[(i-1):i])
            }
            #tmp <- sequence_train[swap_train]
            #sequence_train[swap_train] <- sequence_train[swap_train-1]
            #sequence_train[swap_train-1] <- tmp
        }

        train <- t(sapply(sequence_train, FUN = function(i)
                        mvrnorm(1, mu=mu[i,], Sigma=Sigma[[i]])))
    }

    if(n_test>0) {
        sequence_test <- rep(subseq, n_test/n_subseq)

        ## randomly mess up sequence missing
        if(p_swap>0) {
            swap_test <- which(runif(n_test-1)<p_swap)+1L
            for(i in swap_test) {
                sequence_test[(i-1):i] <-
                rev(sequence_test[(i-1):i])
            }
            #swap_test <- which(runif(n_test-1)<p_swap)+1L
            #tmp <- sequence_test[swap_test]
            #sequence_test[swap_test] <- sequence_test[swap_test-1]
            #sequence_test[swap_test-1] <- tmp
        }

        test <- t(sapply(sequence_test, FUN = function(i)
                        mvrnorm(1, mu=mu[i,], Sigma=Sigma[[i]])))

        ## outliers are random points
        outlier_position <- runif(n_test)<p_outlier
        n_outliers <- sum(outlier_position)
        if(d>1) test[outlier_position, ] <- replicate(d, 
                get_outlier(n_outliers))
        else test[outlier_position] <- get_outlier(n_outliers)

    }   

    list(test=test, train=train, 
            sequence_test=sequence_test, sequence_train=sequence_train, 
            swap_test=swap_test, swap_train = swap_train, 
            outlier_position=outlier_position, 
            model = model)
}

