setMethod("initialize", "RankData",
	function(.Object, ...){
        arg = list(...)
		fields = names(arg)
        # init ranking
        if("ranking" %in% fields){
            .Object@ranking = arg[["ranking"]]
        } else if ("ordering" %in% fields){
            .Object@ranking = OrderingToRanking(arg[["ordering"]])
        } else {
          stop("Either ordering or ranking matrix should be given")
        }
        # init nobj
        if("nobj" %in% fields){
            .Object@nobj = arg[["nobj"]]
        } else {
            .Object@nobj = max(.Object@ranking)
        }
        # init count
        if("count" %in% fields){
            .Object@count = arg[["count"]]
        } else {
            .Object@count = rep(1,nrow(.Object@ranking))
        }
        # init nobs
        if("nobs" %in% fields){
            .Object@nobs = arg[["nobs"]]
            stopifnot(.Object@nobs==sum(.Object@count))
        } else {
            .Object@nobs = sum(.Object@count)
        }
        # init ndistinct
        .Object@ndistinct = nrow(.Object@ranking)
        # handle topq
        if ("topq" %in% fields){  # the field is given
            if (min(arg[["topq"]]) < .Object@nobj-1){  # true topq case
                if (any(arg[["topq"]]>=.Object@nobj)){
                    warning("topq value should range between 1 and nobs-1")
                    arg[["topq"]][arg[["topq"]]>=.Object@nobj] = .Object@nobj-1
                }
                max_rank = apply(.Object@ranking,1,max)
                max_rank = as.numeric(max_rank) - 1
                if (!setequal(max_rank,arg[["topq"]])){
                    warning("The supplied top-q vector is not valid")
                    .Object@topq = unique(max_rank)
                } else {
                    .Object@topq = arg[["topq"]]
                }
                .Object@q_ind = c(1,cumsum(as.numeric(table(max_rank)[as.character(.Object@topq)]))+1)
                .Object@subobs = numeric(length(arg[["topq"]]))
                for (i in 1:length(arg[["topq"]])){
                    .Object@subobs[i] = sum(.Object@count[ .Object@q_ind[i]: (.Object@q_ind[i+1]-1) ])
                }
            } else {
                .Object@q_ind = -1
                .Object@subobs = -1
            }
        } else {
            .Object@topq = .Object@nobj - 1
            .Object@q_ind = -1
            .Object@subobs = -1
        }
        .Object
	}
)




setGeneric("SingleClusterModel",
        def=function(dat,init,ctrl,modal_ranking){standardGeneric("SingleClusterModel")}
)
# 


# single cluster model method for Weighted Kendall Distance
setMethod(
  "SingleClusterModel",
  signature = c("RankData","RankInit","RankControlWeightedKendall"),
  definition = function(dat,init,ctrl,modal_ranking) {
    param_len = max(dat@topq)
    param.coeff = CWeightGivenPi(dat@ranking,modal_ranking)
    param.coeff = matrix(param.coeff,ncol = dat@ndistinct,byrow = TRUE) %*%
      dat@count
    param.coeff = as.numeric(param.coeff)[1:param_len]
    
    if (length(dat@topq) == 1 && dat@topq == dat@nobj-1) {
        
      obj = function(param) {
        a = -1 * param %*% param.coeff - dat@nobs * LogC(c(param,rep(0,dat@nobj -
                                                                       1 - param_len)))
        as.numeric(-1 * a)
      }
      tt = t.gen(param_len)
      gradiant = function(param) {
        grad = GHC(param,tt)
        dat@nobs * grad + param.coeff
      }
      opt_res = optimx::optimx(
        par = init@param.init[[init@clu]][1:param_len],fn = obj,gr = gradiant,lower =
          rep(0,param_len),upper = rep(Inf,param_len),method = "L-BFGS-B",control =
          ctrl@optimx_control
      )
      param.est = unlist(opt_res[1:param_len])
      log_likelihood = -1 * opt_res[[param_len + 1]]
    } else {
      obj = function(param) {
        cond_prob = dat@subobs/dat@nobs
        norm_vec = numeric(length(dat@topq))
        for (i in 1:length(dat@topq)) {
          j = dat@topq[i]
          norm_vec[i] = LogC(c(param[1:j],rep(0,dat@nobj-1 - j))) - lgamma(dat@nobj-j+1) - log(cond_prob[i])
        }
        a = -1 * param %*% param.coeff - dat@subobs %*% norm_vec
        as.numeric(-1 * a)
      }
      opt_res = optimx::optimx(
        par = init@param.init[[init@clu]][1:param_len],fn = obj,lower = rep(0,param_len),upper =
          rep(Inf,param_len),method = "L-BFGS-B",control = ctrl@optimx_control
      )
      param.est = unlist(opt_res[1:param_len])
      log_likelihood = -1 * obj(param.est)
    }
    param.est = c(param.est,rep(0,dat@nobj - 1 - param_len))
    list(
      param.est = param.est,w.est = paramTow(param.est),log_likelihood = log_likelihood
    )
  }
)

# single cluster model method for Kendall distance
setMethod("SingleClusterModel",
    signature = c("RankData","RankInit","RankControlKendall"),
    definition = function(dat,init,ctrl,modal_ranking){
        param.coeff <- FindV(dat@ranking,modal_ranking)
        param.coeff <- rowSums(param.coeff)%*%dat@count
        param.coeff <- as.numeric(param.coeff)
        
        obj <- function(param){
            param*param.coeff + dat@nobs*LogC_Component(rep(param,dat@nobj-1))
        }
        
        opt_res <- stats::optimize(f=obj,interval =c(0,100))
        list(param.est=opt_res$minimum,log_likelihood=-1*opt_res$objective)
    }
)

# single cluster model method for Phi Component Model
setMethod("SingleClusterModel",
    signature = c("RankData","RankInit","RankControlPhiComponent"),
    definition = function(dat,init,ctrl,modal_ranking){
        param.coeff <- FindV(dat@ranking,modal_ranking)
        param.coeff <- t(param.coeff)%*%dat@count
        param.coeff <- as.numeric(param.coeff)
        param_len <- dat@nobj-1
        obj <- list()
        opt_res <- list()
        for (i in 1:param_len){
            obj[[i]] <- function(param){
                rhs <- exp(-param)/(1-exp(-param))-(dat@nobj+1-i)*exp(-(dat@nobj+1-i)*param)/(1-exp(-(dat@nobj+1-i)*param))
                ret <- param.coeff[i] - dat@nobs* rhs
                ret
            }
            opt_res[[i]] <- stats::uniroot(f=obj[[i]],interval=c(0,100),f.lower=-Inf)
        }
        param.est <- vapply(opt_res,function(x)x$root,numeric(1))
        log_likelihood <- -param.coeff%*%param.est - dat@nobs*LogC_Component(param.est)
        list(param.est=param.est,log_likelihood=log_likelihood)
    }
)

setMethod("SingleClusterModel",
          signature = c("RankData","RankInit","RankControlWtau"),
          definition = function(dat, init, ctrl, modal_ranking){
              # the same order as param.expand
              param_len <- dat@nobj
              param.coeff <- Wtau(dat@ranking, modal_ranking)
              param.coeff <- t(param.coeff)%*%dat@count
              param.coeff <- as.numeric(param.coeff)
              allperm <- AllPerms(dat@nobj)
              all_coeff <- Wtau(allperm, modal_ranking)
              obj <- function(param){
                  param.expand <- outer(param,param)[upper.tri(diag(length(param)))]
                  LogC <- log(sum(exp(-1*all_coeff %*% param.expand)))
                  loglike <- param.coeff %*% param.expand + dat@nobs*LogC
                  as.numeric(loglike)
              }
              opt_res = optimx::optimx(
                  par = init@param.init[[init@clu]][1:param_len],fn = obj,
                    method = "Nelder-Mead",control = ctrl@optimx_control
              )
              param.est = unlist(opt_res[1:param_len])
              log_likelihood = -1 * opt_res[[param_len + 1]]
              list(
                  param.est = param.est,log_likelihood = log_likelihood
              )
          }
)

# single cluster model method for Kendall distance
setMethod("SingleClusterModel",
          signature = c("RankData","RankInit","RankControlSpearman"),
          definition = function(dat,init,ctrl,modal_ranking){
              param.coeff <- colSums((t(dat@ranking) - modal_ranking)^2)
              param.coeff <- as.numeric(param.coeff)%*%dat@count
              param.coeff <- 0.5*as.numeric(param.coeff)
              allperm <- AllPerms(dat@nobj)
              all_coeff <- 0.5*as.numeric(colSums((t(allperm) - modal_ranking)^2))
              
              obj <- function(param){
                  LogC <- log(sum(exp(-1 * all_coeff * param)))
                  param*param.coeff + dat@nobs*LogC
              }
              
              opt_res <- stats::optimize(f=obj,interval =c(0,100))
              list(param.est=opt_res$minimum,log_likelihood=-1*opt_res$objective)
          }
)

setMethod("SingleClusterModel",
          signature = c("RankData","RankInit","RankControlFootrule"),
          definition = function(dat,init,ctrl,modal_ranking){
              param.coeff <- apply(dat@ranking, 1, function(x){ sum(abs(x-modal_ranking))} )
              param.coeff <- as.numeric(param.coeff)%*%dat@count
              param.coeff <- as.numeric(param.coeff)
              allperm <- AllPerms(dat@nobj)
              all_coeff <- as.numeric(apply(allperm, 1, function(x){ sum(abs(x-modal_ranking))} ))
              
              obj <- function(param){
                  LogC <- log(sum(exp(-1 * all_coeff * param)))
                  param*param.coeff + dat@nobs*LogC
              }
              
              opt_res <- stats::optimize(f=obj,interval =c(0,100))
              list(param.est=opt_res$minimum,log_likelihood=-1*opt_res$objective)
          }
)

setMethod("SingleClusterModel",
          signature = c("RankData","RankInit","RankControlHamming"),
          definition = function(dat,init,ctrl,modal_ranking){
              param.coeff <- apply(dat@ranking, 1, function(x){sum(x != modal_ranking)} )
              param.coeff <- as.numeric(param.coeff)%*%dat@count
              param.coeff <- as.numeric(param.coeff)
              allperm <- AllPerms(dat@nobj)
              all_coeff <- as.numeric(apply(allperm, 1, function(x){sum(x != modal_ranking)} ))
              
              obj <- function(param){
                  LogC <- log(sum(exp(-1 * all_coeff * param)))
                  param*param.coeff + dat@nobs*LogC
              }
              
              opt_res <- stats::optimize(f=obj,interval =c(0,100))
              list(param.est=opt_res$minimum,log_likelihood=-1*opt_res$objective)
          }
)

setMethod("SingleClusterModel",
          signature = c("RankData","RankInit","RankControlCayley"),
          definition = function(dat,init,ctrl,modal_ranking){
              param.coeff <- FindCayley(dat@ranking, modal_ranking)%*%dat@count
              param.coeff <- as.numeric(param.coeff)
              allperm <- AllPerms(dat@nobj)
              all_coeff <- as.numeric(FindCayley(allperm, modal_ranking))
              
              obj <- function(param){
                  LogC <- log(sum(exp(-1 * all_coeff * param)))
                  param*param.coeff + dat@nobs*LogC
              }
              
              opt_res <- stats::optimize(f=obj,interval =c(0,100))
              list(param.est=opt_res$minimum,log_likelihood=-1*opt_res$objective)
          }
)

setGeneric("FindProb",
        def=function(dat,ctrl,modal_ranking,param){standardGeneric("FindProb")}
)


setMethod("FindProb",
        signature=c("RankData","RankControlWeightedKendall"),
        definition = function(dat,ctrl,modal_ranking,param){
            distance = param %*% matrix(CWeightGivenPi(dat@ranking,modal_ranking),ncol = dat@ndistinct,byrow = TRUE)
            if (length(dat@topq) == 1 && dat@topq == dat@nobj-1) {
                C = exp(LogC(param))
                prob = exp(-1*distance)/C
                if(dat@topq>0){
                    prob = prob*factorial(dat@nobj-dat@topq)
                }
            } else {
                cond_prob = dat@subobs/dat@nobs
                prob = exp(-1*distance)
                for (i in 1:length(dat@topq)) {
                    j = dat@topq[i]
                    norm_c = exp(LogC(c(param[1:j],rep(0,length(param) - j))) - lgamma(dat@nobj-j+1))
                    prob[dat@q_ind[i]:(dat@q_ind[i+1]-1)] = prob[dat@q_ind[i]:(dat@q_ind[i+1]-1)]/norm_c*cond_prob[i]
                }
                
            }
            prob
        }
)

setMethod("FindProb",
        signature=c("RankData","RankControlKendall"),
        definition = function(dat,ctrl,modal_ranking,param){
            param = param[1]
            distance = FindV(dat@ranking,modal_ranking) %*% rep(param,dat@nobj-1)
            C = exp(LogC_Component(rep(param,dat@nobj-1)))
            prob = exp(-1*distance)/C
            prob
        }
)


setMethod("FindProb",
        signature=c("RankData","RankControlPhiComponent"),
        definition = function(dat,ctrl,modal_ranking,param){
            distance = FindV(dat@ranking,modal_ranking) %*% param
            C = exp(LogC_Component(param))
            prob = exp(-1*distance)/C
            prob
        }
)

setMethod("FindProb",
          signature=c("RankData","RankControlWtau"),
          definition<- function(dat,ctrl,modal_ranking,param){
              allperm <- AllPerms(dat@nobj)
              all_coeff <- Wtau(allperm, modal_ranking)
              param.expand <- outer(param,param)[upper.tri(diag(length(param)))]
              C <- sum(exp(-1*all_coeff %*% param.expand))
              distance<- Wtau(dat@ranking, modal_ranking) %*% param.expand
              prob<- exp(-1*distance)/C
              prob
          }
)

setMethod("FindProb",
          signature=c("RankData","RankControlSpearman"),
          definition<- function(dat,ctrl,modal_ranking,param){
              allperm <- AllPerms(dat@nobj)
              all_coeff <- 0.5*as.numeric(colSums((t(allperm) - modal_ranking)^2))
              C <- sum(exp(-1*all_coeff * param))
              param.coeff <- colSums((t(dat@ranking) - modal_ranking)^2)
              param.coeff <- 0.5*as.numeric(param.coeff)
              distance <- param.coeff * param
              prob <- exp(-1*distance)/C
              prob
          }
)


setMethod("FindProb",
          signature=c("RankData","RankControlFootrule"),
          definition<- function(dat,ctrl,modal_ranking,param){
              allperm <- AllPerms(dat@nobj)
              all_coeff <- as.numeric(apply(allperm, 1, function(x){ sum(abs(x-modal_ranking))} ))
              C <- sum(exp(-1*all_coeff * param))
              param.coeff <- apply(dat@ranking, 1, function(x){ sum(abs(x-modal_ranking))} )
              param.coeff <- as.numeric(param.coeff)
              distance <- param.coeff * param
              prob <- exp(-1*distance)/C
              prob
          }
)


setMethod("FindProb",
          signature=c("RankData","RankControlHamming"),
          definition<- function(dat,ctrl,modal_ranking,param){
              allperm <- AllPerms(dat@nobj)
              all_coeff <- as.numeric(apply(allperm, 1, function(x){sum(x != modal_ranking)} ))
              C <- sum(exp(-1*all_coeff * param))
              param.coeff <- apply(dat@ranking, 1, function(x){sum(x != modal_ranking)} )
              param.coeff <- as.numeric(param.coeff)
              distance <- param.coeff * param
              prob <- exp(-1*distance)/C
              prob
          }
)

setMethod("FindProb",
          signature=c("RankData","RankControlCayley"),
          definition<- function(dat,ctrl,modal_ranking,param){
              allperm <- AllPerms(dat@nobj)
              all_coeff <- as.numeric(FindCayley(allperm, modal_ranking))
              C <- sum(exp(-1*all_coeff * param))
              param.coeff <- FindCayley(dat@ranking, modal_ranking)
              distance <- param.coeff * param
              prob <- exp(-1*distance)/C
              prob
          }
)