#' @useDynLib rankdist
#' @importFrom Rcpp sourceCpp
#' @import stats
NULL


#' Calculate Kendall distance matrix between rankings
#' 
#' @param ranking a matrix of rankings
#' @return  Kendall distance matrix between rankings. The value in ith row and jth column is the
#' Kendall distance between the ith and jth rankings.
#' @export
DistanceMatrix <- function(ranking){
  n = ncol(ranking)
  w = c(rep(0,n-2),1)
  distcalc = function(x){
    fai.coeff = matrix(CWeightGivenPi(ranking,x),nrow=nrow(ranking))
    as.numeric(fai.coeff %*% w)
  }
  apply(ranking,1,distcalc)
}

#' Calculate Kendall distance between one ranking and a matrix of rankings
#' @param mat a matrix of rankings. Each row stores a ranking.
#' @param r a single ranking
#' @return a vector of Kendall distances
#' @export
DistanceBlock <- function(mat,r){
  n = ncol(mat)
  w = c(rep(0,n-2),1)
  fai.coeff = matrix(CWeightGivenPi(mat,r),nrow=nrow(mat))
  as.numeric(fai.coeff %*% w)
}

#' Calculate Kendall distance between a pair of rankings
#' @param r1 a single ranking
#' @param r2 a single ranking
#' @return Kendall distance value
#' @export
DistancePair <- function(r1,r2){
  n = length(r1)
  w = c(rep(0,n-2),1)
  fai.coeff = matrix(CWeightGivenPi(r1,r2),nrow=1)
  as.numeric(fai.coeff %*% w)
}

#' Print a brief summary of the fitted model
#' 
#' Print a brief summary of the fitted model. This includes information about goodness
#' of fit as well as parameter estimation.
#' @param model a ranking model returned by a call to RankDistanceModel function.
#' @export
ModelSummary <- function(model){
  clus = rev(order(model$p))
  nobj = length(model$modal_ranking.est[[1]])
  mod_name = deparse(substitute(model))
  cat("Summary of",mod_name,"\n")
  cat(rep("=",nchar(mod_name)+10),"\n",sep="")
  cat("Goodness of Fit\n")
  cat("SSR:\t",model$SSR,"\n")
  cat("BIC:\t",model$BIC,"\n")
  cat("dof:\t",model$free_params,"\n")
  cat("Parameter Estimation\n")
  cat("Cluster",LETTERS[1:nobj],"p","Parameters\n",sep="\t")
  for (i in 1:length(clus)){
    cat(i,' ',' ',' ',model$modal_ranking.est[[clus[i]]],round(model$p[clus[i]],digits=2),round(model$w.est[[clus[i]]],digits=2),"\n",sep="\t")
  }
}

#' Generate simple examples
#' 
#' This function generates simple examples for illustrative proposes.
#' The sample contains rankings of five objects and the underlying model is a Mallows' phi
#' model with default dispersion parameter set to 0.2
#' and modal ranking set to (1,2,3,4,5)
#' @param ranking TRUE if "ranking" representation is used in the output data; otherwise "ordering" representation is used.
#' @param central The modal ranking.
#' @param lambda The parameter in the model.
#' @export
GenerateExample <- function(ranking=TRUE, central=1:5, lambda=0.2){
  rankings <- rbind(1:5,permute::allPerms(5))
  # central <- 1:5
  Kdist <- DistanceBlock(rankings,central)
  # lambda <- 0.2
  prob <- exp(-lambda*Kdist)
  prob <- prob/sum(prob)
  indx <- sample(1:120, 2000,replace=TRUE,prob=prob)
  count <- as.numeric(table(indx))
  if (ranking){
    return(list(ranking=rankings,count=count))
  } else {
    return(list(ordering=OrderingToRanking(rankings),count=count))
  }
}

#' Generate simple examples of top-q rankings
#' 
#' This function generates simple examples for illustrative proposes.
#' The sample contains the top-3 rankings of five objects and the underlying model is a weighted Kendall distance model
#' with default weights set to (0.7,0.5,0.3,0)
#' and modal ranking set to (1,2,3,4,4)
#' @param central The modal ranking.
#' @param w The weights in the model.
#' @export
GenerateExampleTopQ <- function(central=c(1,2,3,4,4),w=c(0.7,0.5,0.3,0)){
  prankings <- rbind(1:5,permute::allPerms(5))
  prankings[prankings>3] <- 4
  prankings <- prankings[!duplicated(prankings),]
  fai <- wToparam(w)
  distmat <- matrix(CWeightGivenPi(prankings,central),ncol=nrow(prankings),byrow=TRUE)
  distvec <- as.numeric(fai%*%distmat)
  probs <- exp(-distvec)/sum(exp(-distvec))
  samples <- sample(1:nrow(prankings),10000,replace=TRUE,prob=probs)
  count <- table(samples)
  count <- as.numeric(count)
  return(list(ranking=prankings,count=count))
}

#' Create Hash Value for Ranking
#'
#' Sometimes it is handy to deal with rankings as a hash value. \code{RanktoHash} returns hash values for ranks. 
#' Maximum 52 objects are supported.
#' @param r  a vector or matrix of rankings. Each row of the matrix represents a ranking.
#'   The ranking should be a integer from one to number of objects. No NA is allowed
#' @return a vector of character strings representing the hash values.
#' @seealso \code{\link{HashtoRank}} for a reverse operation.
#' @export
RanktoHash <- function(r){
    char <- c(letters,LETTERS)
    if (!is.matrix(r)){
        a <- paste(char[r],collapse="")
    } else {
        a <- apply(r,1,function(x)paste(char[x],collapse=""))
    }
    as.vector(a,mode="character")
}


#' Obtain Ranking from Hash Value
#'
#' \code{HashToRank} returns rankings from given hash values. Maximum 52 objects are supported.
#' @param h  A vector of hash values. 
#' @return a matrix of rankings if input has more than one element or a single ranking (numeric vector) if input has only one element
#' @seealso \code{\link{RanktoHash}} for a reverse operation.
#' @export
HashtoRank <- function(h){
    char <- c(letters,LETTERS)
    transvec <- seq_along(char)
    names(transvec) <- char
    nobj <- nchar(h[1])
    numrank <- transvec[unlist(strsplit(h,split=""))]
    if (length(h)==1){
        as.vector(numrank,mode="numeric")
    } else {
        matrix(data=numrank,ncol=nobj,byrow=TRUE)
    }
    
}


#' Transformation between Rankings and Orderings
#' 
#' \code{OrderingToRanking} transforms between ranking representation and ordering representation.
#'
#'    Ranking representation encodes the position of objects. Ordering representation is an ordered sequence of objects.
#'    For example ranking (2 3 1 4) is equivalent to ordering (3 1 2 4), which means object 3 is first, object 1 is second, followed by object 2 and 4.
#'    Also note that we can use this function to transform rankings into orderings, and applying this function twice will not change the input value.
#' @param ordering  a matrix of orderings or rankings. Each row contains an observation.
#' @return  a matrix of transformed rankings or orderings. Each row contains an observation.
#' @export
OrderingToRanking <- function(ordering){
    if (is.matrix(ordering)){
        ranking <- t(apply(ordering,1,order))
    } else {
        ranking <- order(ordering)
    }
    ranking
}



#' Fit A Mixture of Distance-based Models
#' 
#' \code{RankDistanceModel} fits ranking models based on inputs
#' 
#' The procedure will estimate central rankings, the probability of each cluster and weights.
#' 
#' @param dat A \linkS4class{RankData} object.
#' @param init A \linkS4class{RankInit} object.
#' @param ctrl A \linkS4class{RankControl} object.
#' 
#' @return A list containing the following components:
#' \describe{
#' \item{\code{modal_ranking.est}}{the estimated modal ranking for each cluster.}
#' \item{\code{p}}{the marginal probability of each cluster.}
#' \item{\code{w.est}}{the estimated weights of each cluster.}
#' \item{\code{param.est}}{the phi parametrisation of weights of each cluster (for Weighted Kendall model only).}
#' \item{\code{SSR}}{the sum of squares of Pearson residuals}
#' \item{\code{log_likelihood}}{the fitted log_likelihood}
#' \item{\code{BIC}}{the fitted Bayesian Information Criterion value}
#' \item{\code{free_params}}{the number of free parameters in the model}
#' \item{\code{expectation}}{the expected value of each observation given by the model}
#' \item{\code{iteration}}{the number of EM iteration}
#' \item{\code{model.call}}{the function call}
#' }
#' @seealso \code{\link{RankData}}, \code{\link{RankInit}}, \code{\link{RankControl}}
#' @export
RankDistanceModel <- function(dat,init,ctrl){
    flag_transform = (min(dat@topq) != dat@nobj - 1)
    if (class(ctrl) == "RankControlWeightedKendall"){
        flag_transform = (ctrl@assumption == "equal-probability")
    }
    if (flag_transform){
        dat_org = dat
        dat = BreakTieEqualProb(dat)
    }
    # get parameters
    tt <- dat@nobj # number of objects
    n <- dat@nobs    # total observations
    distinctn <- dat@ndistinct	# total distinct observations
    clu <- init@clu
    count <- dat@count
    p_clu <- matrix(ncol = distinctn,nrow = clu)
    # setting up return values
    modal_ranking.est <- list()
    modal_ranking.est.last <- list()
    p <- rep(1/clu,clu)
    p.last <- numeric(clu)
    param <- list()
    param.last <- list()
    func.call <- match.call()
    loopind=0

    if (clu > 1L){  # use EM to fit multi-cluster model
        # further initialization
        modal_ranking.est <- init@modal_ranking.init
        param <- init@param.init
        init.clu = list()
        for (i in 1:clu){
            init.clu[[i]] <- new("RankInit",param.init=list(param[[i]]),modal_ranking.init=list(modal_ranking.est[[i]]),clu=1L)
        }
        
        while(TRUE){
            loopind <- loopind+1
            # E step
            z <- matrix(nrow = distinctn,ncol = clu)
            for (i in 1:clu){
                z[,i] <- p[i]*FindProb(dat,ctrl,modal_ranking.est[[i]],param[[i]])
            }
            sums <- rowSums(z)
            z <- z/sums
            
            # M step
            p <- t(z) %*% count
            p <- p/sum(p)
            if (any(p<0.03)){
                warning("One cluster has probability smaller than 3%; Try fewer clusters")
                return (list(p=p,modal_ranking.est=modal_ranking.est,param=param,iteration=loopind))
            }
            for ( i in 1:clu){
                dat.clu <- UpdateCount(dat, z[,i] * count)
                init.clu[[i]]@param.init <- list(param[[i]])
                init.clu[[i]]@modal_ranking.init <- list(modal_ranking.est[[i]])
                solve.clu <- SearchPi0(dat.clu,init.clu[[i]],ctrl)
                modal_ranking.est[[i]] <- solve.clu$pi0.ranking
                param[[i]] <- solve.clu$param.est
                # change here: conditional prob for each cluster
                # p_clu[i,] <- FindProb(dat,ctrl,modal_ranking.est[[i]],param[[i]])*p[i]
                p_clu[i,] <- FindProb(dat.clu,ctrl,modal_ranking.est[[i]],param[[i]])*p[i]
            }
            log_likelihood_clu <- sum(log(colSums(p_clu))%*%dat@count)
            # break?
            if(loopind != 1){
                if (loopind < ctrl@EM_limit) {
                    cond1 <- all( p.last - p < ctrl@EM_epsilon)
                    cond2 <- all( unlist(modal_ranking.est.last) - unlist(modal_ranking.est) < ctrl@EM_epsilon)
                    cond3 <- all( unlist(param.last) - unlist(param) < ctrl@EM_epsilon)
                    if (cond1 && cond2 && cond3){
                        break
                    }else if(abs(log_likelihood_clu.last - log_likelihood_clu)<0.01){
                        break
                    }
                } else {
                    print(paste("Algorithm did not converge in",ctrl@EM_limit,"iterations"))
                    return(NULL)
                }
            }
            p.last <- p
            modal_ranking.est.last <- modal_ranking.est
            param.last <- param
            log_likelihood_clu.last <- log_likelihood_clu
        } # inf loop
        
        est_prob <- colSums(p_clu)
    } else {  # fit single cluster model
        sigle_cluster_mod <- SearchPi0(dat,init,ctrl)
        modal_ranking.est <- list(sigle_cluster_mod$pi0.ranking)
        p <- 1
        param <- list(sigle_cluster_mod$param.est)
        log_likelihood_clu.last <- sigle_cluster_mod$log_likelihood
        est_prob <- FindProb(dat,ctrl,modal_ranking.est[[1]],param[[1]])*p[1]
    }
    # finishing up with parameter estimation and summary statistics
    p <- as.numeric(p)
    v <- length(p) - 1 + sum(unlist(param)!=0)
    if (flag_transform){
        full_prob = rep(0, length(dat_org@count))
        # aggregate probability
        # iterate through different q
        for (i in 1:length(dat_org@topq)){
            ind_start <- dat_org@q_ind[i]
            ind_end <- dat_org@q_ind[i+1] - 1
            this_q <- dat_org@topq[i]
            ranking_q <- dat@ranking
            ranking_q[ranking_q > this_q] <- this_q + 1
            hash_q = RanktoHash(ranking_q)
            # iterate through partial rankings
            for (ind_partial in ind_start:ind_end){
                this_partial <- dat_org@ranking[ind_partial, ]
                this_hash <- RanktoHash(this_partial)
                ind_agg<- which(hash_q == this_hash)
                prob_agg<- sum(est_prob[ind_agg])
                full_prob[ind_partial] <- prob_agg
            }
            # conditional
            p_q = dat_org@subobs[i] / dat_org@nobs
            full_prob[ind_start:ind_end] = full_prob[ind_start:ind_end]*p_q
        }
        log_like<- dat_org@count %*% log(full_prob)
        bic <- -2*log_like + v*log(n)
        expectation <- as.numeric(full_prob*n)
        SSR <- sum((expectation-dat_org@count)^2/expectation)
    } else {
        
        log_like<- dat@count %*% as.numeric(log(est_prob))
        log_likelihood_clu.last - log_like
        
        bic <- -2*log_like + v*log(n)
        expectation <- as.numeric(est_prob*n)
        SSR <- sum((expectation-dat@count)^2/expectation)
    }
    ret <- list(p=p,modal_ranking.est=modal_ranking.est,w.est=lapply(param,paramTow),
                param.est=param,SSR=SSR,log_likelihood=log_likelihood_clu.last,free_params=v,
                BIC=bic,expectation=expectation,iteration=loopind,model.call=func.call)
    ret
}


#' Find Initial Values of phi
#'
#' \code{MomentsEst} finds the initial values of phi which can be used in the subsequent optimization problems.
#' Linear model is fitted to the log odds of rankings.This function is only useful to the Weighted Kendall model.
#' 
#' @param dat a RankData object
#' @param size the number of samples to take in the linear model
#' @param pi0 an optional argument showing the location of central ranking. 
#' If not provided, Borda Count method is used to estimate the central ranking.
#' @return estimated phi
#' @examples MomentsEst(apa_obj,40)
#' @export
MomentsEst <- function(dat,size,pi0=NULL){
    # estimating pi0
    avg_rank <- dat@count %*% dat@ranking;
    if (is.null(pi0)){
        modal_ranking <- sort(OrderingToRanking(dat@ranking[1,]))[order(avg_rank)]
    } else {
        modal_ranking <- pi0
    }
    nparam <- dat@nobj - 1
    # construct equation set
    prob_obs <- dat@count/dat@nobs
    # row-wise
    distance_mat <- matrix(CWeightGivenPi(dat@ranking,modal_ranking),nrow = dat@ndistinct)
    pair <- sample(1:nrow(distance_mat),2*size,replace=TRUE)
    pair_mat <- matrix(pair,nrow=2)
    logodd <- numeric(size)
    design_mat <- matrix(ncol=nparam, nrow=size)
    for (i in 1:size){
        logodd[i] <- prob_obs[pair_mat[1,i]]/prob_obs[pair_mat[2,i]]
        design_mat[i,] <- distance_mat[pair_mat[1,i],]-distance_mat[pair_mat[2,i],]
    }
    param.est <- stats::lm(logodd~design_mat-1)
    param.est <- param.est$coefficients
    param.est[param.est<0] <- 0
    names(param.est) <- NULL
    param.est
}

#' American Psychological Association (APA) election data
#'
#' A dataset containing 5738 complete votes in APA election. There are 5 candidates in total.
#'
#' @format a RankData object
#' 
#' @references Marden, J. I. (1995). Analyzing and Modeling Rank Data (94-96). Chapman Hall, New York.
#' @seealso \code{\link{apa_partial_obj}}
"apa_obj"

#' American Psychological Association (APA) election data (partial rankings included)
#'
#' A dataset containing 5738 complete votes and 9711 partial votes in APA election. There are 5 candidates in total.
#'
#' @format a RankData object
#' 
#' @references Marden, J. I. (1995). Analyzing and Modeling Rank Data (94-96). Chapman Hall, New York.
#' @seealso \code{\link{apa_obj}}
"apa_partial_obj"

