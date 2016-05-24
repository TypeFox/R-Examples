#' Bivariate correlations including randomized response variables
#' 
#' \code{RRcor} calculates bivariate Pearson correlations of variables measured with or without RR.
#' @param x a numeric vector, matrix or data frame.
#' @param y \code{NULL} (default) or a vector, matrix or data frame with compatible dimensions to \code{x}. 
#' @param models a vector defining which RR design is used for each variable. Must be in the same order as variables appear in \code{x} and \code{y} (by columns). Available discrete models: \code{Warner}, \code{Kuk}, \code{FR}, \code{Mangat}, \code{UQTknown}, \code{UQTunknown}, \code{Crosswise}, \code{SLD} and \code{direct} (i.e., no randomized response design). Available continuous models: \code{mix.norm}, \code{mix.exp}.
#' @param p.list a \code{list} containing the randomization probabilities of the RR models defined in \code{models}. Either, all \code{direct}-variables (i.e., no randomized response) in \code{models} can be excluded in \code{p.list}; or, if specified, randomization probabilities \code{p} are ignored for \code{direct}-variables. See \code{\link{RRuni}} for a detailed specification of p.
#' @param group a matrix defining the group membership of each participant (values 1 and 2) for all multiple group models(\code{SLD}, \code{UQTunknown}). If only one of these models is included in \code{models}, a vector can be used. For more than one model, each column should contain one grouping variable 
#' @param bs.n number of samples used to get bootstrapped standard errors
#' @param bs.type to get boostrapped standard errors, use \code{"se.p"} for the parametric and/or \code{"se.n"} for the nonparametric bootstrap. Use \code{"pval"} to get p-values from the parametric bootstrap (assuming a true correlation of zero). Note that \code{bs.n} has to be larger than 0. The parametric bootstrap is based on the assumption, that the continuous variable is normally distributed within groups defined by the true state of the RR variable. For polytomous forced response (FR) designs, the RR variable is assumed to have equally spaced distances between categories (i.e., that it is interval scaled)
#  @param pval if \code{TRUE}, p-values are obtained from a separate parametric bootstrap using the null hypothesis that the true correlation is zero (i.e., cor=0; note that \code{bs.n} must be larger than zero)
#' @param nCPU number of CPUs used for the bootstrap
#' @return \code{RRcor} returns a list with the following components:: 
#' 
#' \code{r} estimated correlation matrix
#' 
#' \code{rSE.p}, \code{rSE.n} standard errors from parametric/nonparametric bootstrap
#' 
#' \code{prob} two-sided p-values from parametric bootstrap
#' 
#' \code{samples.p}, \code{samples.n} sampled correlations from parametric/nonparametric bootstrap (for the standard errors)
#' 
#' @details Correlations of RR variables are calculated by the method of Fox & Tracy (1984) by interpreting the variance induced by the RR procedure as uncorrelated measurement error. Since the error is independent, the correlation can be corrected to obtain an unbiased estimator.
#' 
#' Note that the continuous RR model \code{mix.norm} with the randomization parameter \code{p=c(p.truth, mean, SD)} assumes that participants respond either to the sensitive question with probability \code{p.truth} or otherwise to a known masking distribution with known mean and SD. The estimated correlation only depends on the mean and SD and does not require normality. However, the assumption of normality is used in the parametric bootstrap to obtain standard errors.
#' 
#' @seealso \code{vignette('RRreg')} or \url{https://dl.dropboxusercontent.com/u/21456540/RRreg/index.html} for a detailed description of the RR models and the appropriate definition of \code{p} 
#' @references Fox, J. A., & Tracy, P. E. (1984). Measuring associations with randomized response. \emph{Social Science Research, 13}, 188-197.
#' @examples 
#' # generate first RR variable
#' n <-1000
#' p1 <- c(.3,.7)
#' gData <- RRgen(n,pi=.3,model="Kuk",p1)
#' 
#' # generate second RR variable
#' p2 <- c(.8,.5)
#' t2 <- rbinom(n=n, size=1, prob=(gData$true+1)/2)
#' temp <- RRgen(model="UQTknown",p=p2, trueState=t2)
#' gData$UQTresp <- temp$response
#' gData$UQTtrue <- temp$true
#' 
#' # generate continuous covariate
#' gData$cov <- rnorm(n,0,4) + gData$UQTtrue + gData$true
#' 
#' # estimate correlations using directly measured / RR variables
#' cor(gData[,c("true","cov","UQTtrue")])
#' RRcor(x=gData[,c("response","cov","UQTresp")],
#'       models=c("Kuk","d","UQTknown"),p.list= list(p1,p2) )
#' @import doParallel
#' @import parallel
#' @import foreach
#' @export
RRcor <- function(x, y=NULL, models, p.list, group=NULL, bs.n=0, bs.type=c("se.n","se.p","pval"), nCPU=1){
  
  # input handling  (abkopiert von function 'cor' ; na.method funktioniert nicht
  #   na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", 
  #                              "everything", "na.or.complete"))   ## returns integer!
  #   if (is.na(na.method)) 
  #     stop("invalid 'use' argument")
  
  if (is.data.frame(y) || is.vector(y) || is.matrix(y) || is.numeric(y)){   
    y.name <- deparse(substitute(y))
    y <- as.matrix(y)
    my <- ncol(y)
    ynames <- colnames(y)
    if  ( is.null(ynames) & my==1){
      ynames <- y.name
    }else if( is.null(ynames)){
      ynames <- 1:my
      ynames <- paste( "y", ynames,sep="")
    }
  }
  if (is.data.frame(x)|| is.vector(x) || is.matrix(x) || is.numeric(y)) {
    x.name <- deparse(substitute(x))
    x <- as.matrix(x)
    mx <- ncol(x)
    xnames <- colnames(x)
    if  ( is.null(xnames) & mx==1){
      xnames <- x.name
    }else if( is.null(xnames)){
      xnames <- 1:mx
      xnames <- paste( "x", xnames,sep="")
    }
  }
  if (!is.matrix(x) && is.null(y)) 
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  if (!(is.numeric(x) || is.logical(x))) 
    stop("'x' must be numeric")
  stopifnot(is.atomic(x))
  if (!is.null(y)) {
    if (!(is.numeric(y) || is.logical(y))) 
      stop("'y' must be numeric")
    stopifnot(is.atomic(y))
  }
  
  # erzeuge gemeinsamen datensatz X
  if (!is.null(x) && !is.null(y)){
    if ( !is.null(mx) && !is.null(my) && mx+my <2){
      stop("Not enough variables specified")
    }
    X <- cbind(x,y)
    colnames(X) <- c(xnames,ynames)
    m <- mx+my
  }
  
  if (!is.null(x) && is.null(y)){
    X <- x
    colnames(X) <- xnames
    m <- mx
  }
  if (is.null(x) && !is.null(y)){
    X <- y
    colnames(X) <- ynames
    m <- my
  }
  modelNames <- c("Warner","UQTknown","UQTunknown","Mangat","Kuk","FR","Crosswise","direct","SLD", "mix.norm","mix.exp", "mix.unknown")
  models <- pmatch(models, modelNames, duplicates.ok=T )
  models <- modelNames[models]
  n <- nrow(X)
  
  # adjust length and entries of p.list
  if (length(p.list)<length(models) && length(p.list) == sum(models!="direct") ){
    cnt <- 1
    p.list2 <- list()
    for (i in 1:length(models)){
      if(models[i]=="direct"){
        p.list2[[i]]  <- 1
      } else {
        p.list2[[i]] <-  p.list[[cnt]]
      }
      cnt <- ifelse(models[i]=="direct",cnt,cnt+1)
    }
    p.list <- p.list2
  }

  if (!missing(group) && !is.null(group)) group <- as.matrix(group)
  group.2g <- group
  RRcheck.cor(X,m,models,p.list,colnames(X),group)
  
  # generate matrix with groups and calculate additional parameters like gamma, t, piUQ
  group.mat <- matrix( 1, nrow=n,ncol=m)
  par2 <- list() # rep(NA,m)   # additional parameters except pi
  pi <- rep(NA,m)
  piVar <- rep(NA,m)
  cnt <- 1
  for (i in 1: m){
    if (models[i] == "direct"){
      pi[i] <- mean(X[,i])
      piVar[i] <- var(X[,i])
    }else{
      if (is2group(models[i])){
        group.mat[,i] <-  as.numeric(  group[,cnt])  
        cnt <- cnt+1
      }
      est <- RRuni(X[,i],model=models[i],p=p.list[[i]],group=group.mat[,i],MLest=F)
      if(models[i] =="FR"){
        scale <- 1:length(est$pi)-1
        
         pi[i] <- sum(scale*est$pi)
         piVar[i] <- sum(est$pi*scale^2) - pi[i]^2
      }else{
        pi[i] <- est$pi
        piVar[i] <- est$pi * (1-est$pi)
      }
      switch (models[i],
              #               "CDM" = {par2[i] <- est$gamma},
              #               "CDMsym" = {par2[i] <- est$gamma},
              "UQTunknown" ={par2[[i]] <- est$piUQ},
              "SLD" = {par2[[i]] <- est$t},
              "mix.unknown" = {
                par2[[i]] <- c(est$piUQ, est$y.var)
                piVar[i] <- est$x.var
                })
    } 
  }
  
  #   if (ncol(group) != n) group <- group.mat
  group <- group.mat
  
  
  # x.est : estimated true values for respondents (dependent on group)
  x.est <- X
  # for SLD and UQTunknown: which group to select (always with larger p[i] !)
  group.sel <- rep(1,m)
  
  # used sample sizes for all pairwise comparisons
  n.mat <- diag(m) * n
  g.list <- list()
  quotient.list <- list()
  gN.list <- list()
  cnt <- 0
  
  # quotient contains the scaling factors for the observed correlation
  # rows: different models; second column for 2-group models
  quotient <- rep(0,m)
  for (i in 1:m){
    # pairwise comparison to get adequate subgroups!
    if (i != m){
      for (j in (i+1):m){
        # find all available subgroups
        subgroups <- group[,c(i,j)][!duplicated(group[,c(i,j)]),,drop=F]
        numsub <- nrow(subgroups)
        cnt <- cnt+1
        gN.list[[cnt]] <- rep(0, numsub)
        for (ss in 1:numsub){
          gN.list[[cnt]][ss] <- sum(apply(group[,c(i,j)], 1, function(x) all(x == subgroups[ss,])))
        }
#         # calculate information index per subgroup:
#         inf.idx <- rep(0, numsub)
#         for (ss in 1:numsub){
#           # get the subgroup sample size n[i,j]
#           subgroup.n <- sum(apply(group[,c(i,j)]  == subgroups[ss,], 1, any))
#           pi <- p.list[[i]][subgroups[ss,1]]
#           pj <- p.list[[j]][subgroups[ss,2]]
#           inf.idx[ss] <- subgroup.n * pi * pj
#         }
#         # g.list contains the most informative subgroups (take care of the right order!)
        g.list[[cnt]] <- subgroups #[inf.idx == max(inf.idx),,drop=F]  
      }
    }

    
    p <- p.list[[i]]
    x.est[,i] <- get.x.est(X[,i],models[i],p,
                           par2[[i]],group[,i])
        
    # for each RR variable: get quotient for 2-group models: 2 quotients
    quotient.list[[i]] <- 0
    for (gg in 1:length(unique(group[,i]))){
      sel <- group[,i] == gg
      quotient.list[[i]][gg] <- get.quotient(X[sel,i],models[i],p, par2[[i]], group=gg, pi[i], piVar[i])
    }
#   print(quotient.list)
    
    # for 2-group models: separately
    sel <- rep(T,n)
    if ( is2group(models[i])){
      if (p[2] > p[1]){
        group.sel[i] <- 2
        sel <- group[,i] == 2
      }else{
        sel <- group[,i] == 1
      }
    }
    quotient[i] <- get.quotient(X[sel,i],models[i],p, par2[[i]], group=group.sel[i], pi[i], piVar[i])
  } 
# print(g.list)
# print(gN.list)
# print(g.list)
## loop to compute weights for averaging
cnt <- 0
weights <- list()
for (i in 1: (m-1)){
  for (k in (i+1):m){
    cnt <- cnt +1
    numsubs <- nrow(g.list[[cnt]])
    weights[[cnt]] <- rep(0, numsubs)
    for (gg in 1:numsubs){
      g.i <- g.list[[cnt]][gg,1]
      g.j <- g.list[[cnt]][gg,2]
      sel.i <- group[,i] == g.i
      sel.j <- group[,k] == g.j
      sel <- sel.i & sel.j
      weights[[cnt]][gg] <- sum(sel)/((1+ quotient.list[[i]][g.i])*(1+ quotient.list[[j]][g.j]))
    }
  }
}
wsum <- unlist(lapply(weights, sum))
 
  r <-  diag(m) 
  cnt <- 0
  # calculate correlation separately for each combination of variables
  for (i in 1: (m-1)){
    for (k in (i+1):m){
      # for multiple group models: go through all subgroups listed in g.list!
      cnt <- cnt +1
      for (gg in 1:nrow(g.list[[cnt]])){
#         # select subgroup
        g.i <- g.list[[cnt]][gg,1]
        g.j <- g.list[[cnt]][gg,2]
        sel.i <- group[,i] == g.i
        sel.j <- group[,k] == g.j
        sel <- sel.i & sel.j
#         print(sum(sel))
#         
        r.obs <- cor(x.est[sel,i],x.est[sel,k])
#       varianz <-   c((1+ quotient.list[[i]]) %o% (1+ quotient.list[[j]]))
#       wsum <- sum (gN.list[[cnt]]/varianz)
#       weight <- sum(sel)/(wsum*(1+ quotient.list[[i]][g.i])*(1+ quotient.list[[j]][g.j]))
      weight <- weights[[cnt]][gg]/wsum[cnt]
      r[i,k] <- r[i,k] + weight* r.obs * sqrt( (1+ quotient.list[[i]][g.i])*(1+ quotient.list[[k]][g.j]) )
#         print(r.obs * sqrt( (1+ quotient.list[[i]][g.i])*(1+ quotient.list[[k]][g.j]) ))
#         r[i,j] <- r.sign( r[i,j], models[i], p.list[[i]], models[j], p.list[[j]]) 
        n.mat[i,k] <- n.mat[i,k] + sum(sel)
      }
      # all possible group combinations
      # for multiple group models
#       sel <- group[,i] == group.sel[i] & group[,k] == group.sel[k]      
#       n.mat[i,k] <- sum(sel)

      #         sel <- group[,i] == gg[g,i] &  group[,k] == gg[g,k]
#       r.obs <- cor(x.est[sel,i],x.est[sel,k])
#       r[i,k] <-   r.obs * sqrt( (1+ quotient[i])*(1+ quotient[k]) ) # +  r[i,k] ??
# print("cor:")
# print(r.obs * sqrt( (1+ quotient[i])*(1+ quotient[k]) ))
#       r[i,k] <- r.sign( r[i,k],models[i],p.list[[i]] , models[k],p.list[[k]]) 
      #       }
    }
  }

r[lower.tri(r)] <- t(r)[lower.tri(r)]
n.mat[lower.tri(n.mat)] <- t(n.mat)[lower.tri(n.mat)]
rownames(r) <- colnames(X)
colnames(r) <- colnames(X)

if (sum( r < -1, na.rm=T) >0){
#   warning("The corrected correlation was smaller than -1 and thus set to -1.")
  r[-r >1] <- -1
}
if (sum( r > 1, na.rm=T) >0){
#   warning("The corrected correlation was larger than 1 and thus set to 1.")
  r[r>1] <- 1
}

rownames(n.mat) <- colnames(X)
colnames(n.mat) <- colnames(X)
# character vector including randomization probabilities
p.char <- character(m)
for (i in 1:m){
  if (models[i]!="direct"){
    for (k in 1:length(p.list[[i]])){
      p.char[i] <- paste(p.char[i], round(p.list[[i]][k],4))
    }
  } 
}
  
  rSE.p <- matrix(NA,nrow=m, ncol=m, dimnames=list(colnames(X), colnames(X))) 
  r.Prob <- rSE.p ; rSE.n <- rSE.p;
  prob <- rSE.p ;
  # parametric bootstrap: SE
  if(bs.n>0){ # && "se.p" %in% bs.type){
    if (any(models %in% c("mix.norm","mix.exp")))
      warning("Parametric boostrap not available for continuous mixture RR models. Use the nonparametric bootstrap instead: bs.type='se.n'")
    samples.p <- array(NA, dim=c(m,m,nCPU*ceiling(bs.n/nCPU)))
    for (i in 1:(m-1)){
      for(j in (i+1):m){
#         samples.p[i,i,] <- 1
#         samples.p[j,j,] <- 1
#         cat("variables",1,"and",2,": ",paste( models[c(i,j)]))
        ci <- c(1,1); cj <- c(1,1);
        if (models[i] == "SLD") ci <- c(par2[[i]],1)
        if (models[j] == "SLD") cj <- c(par2[[j]],1)
        #         if (models[i] == "CDM") ci <- rep(1-par2[i], 2)
        #         if (models[j] == "CDM") cj <- rep(1-par2[j], 2)
        groupRatios <- 2 - colMeans(group.mat)
        
        # how many RR variables in pairwise bootstrap?
        if(sum(models[c(i,j)] == "direct") <2){
          if(models[i] != "direct")
            RRi <- RRuni(X[,i], model =models[i], p=p.list[[i]], group = group.mat[,i],MLest = T)
          if(models[j] != "direct")
            RRj <- RRuni(X[,j], model =models[j], p=p.list[[j]], group = group.mat[,j],MLest = T)
          # two RRs
          if(sum(models[c(i,j)] == "direct") ==0){
            if("se.p" %in% bs.type) 
              mcsim <- RRsimu(numRep=bs.n, n=n, pi=c(RRi$pi, RRj$pi), model = models[c(i,j)],
                            p=p.list[c(i,j)], cor=r[i,j], complyRates=list(ci, cj), sysBias=c(0,0),
                            groupRatio=groupRatios[c(i,j)], method="RRcor", MLest=T, nCPU=nCPU)
            if("pval" %in% bs.type){
              mcsim.h0 <- RRsimu(numRep=bs.n, n=n, pi=c(RRi$pi, RRj$pi), model = models[c(i,j)],
                              p=p.list[c(i,j)], cor=0, complyRates=list(ci, cj), sysBias=c(0,0),
                              groupRatio=groupRatios[c(i,j)], method="RRcor", MLest=T, nCPU=nCPU)
            }
          }else{
            # one RR, one nonRR variable
            sel <- grep("direct", models[c(i,j)])
            idx <- ifelse(sel==2,i,j)  # select RR variable
            if("se.p" %in% bs.type)
              mcsim <- RRsimu(numRep=bs.n, n=n, pi=ifelse(idx==i, RRi$pi, RRj$pi), model = models[idx],
                            p=unlist(p.list[idx]), cor=r[i,j], complyRates=ifelse(rep(idx,2)==i,ci, cj),
                            sysBias=c(0,0),
                            groupRatio=groupRatios[idx], method="RRcor", MLest=T, nCPU=nCPU)
            if("pval" %in% bs.type){
              mcsim.h0 <- RRsimu(numRep=bs.n, n=n, pi=ifelse(idx==i, RRi$pi, RRj$pi), model = models[idx],
                              p=unlist(p.list[idx]), cor=0, complyRates=ifelse(rep(idx,2)==i,ci, cj),
                              sysBias=c(0,0),
                              groupRatio=groupRatios[idx], method="RRcor", MLest=T, nCPU=nCPU)
            }
          }
          if ("se.p" %in% bs.type)
            rSE.p[i,j] <- mcsim$results["cor.RRcor","SE"]
          #           r.Prob[i,j] <- sum(mcsim$parEsts[,"cor.RRcor"]>r[i,j])/bootstrap 
        }
        if ("se.p" %in% bs.type){
          samples.p[i,j,] <- mcsim$parEsts[,"cor.RRcor"]
          samples.p[j,i,] <- mcsim$parEsts[,"cor.RRcor"]
        }
        if ("pval" %in% bs.type){
          pcum <- ecdf(mcsim.h0$parEsts[,"cor.RRcor"])(r[i,j])
          prob[i,j] <- 2*ifelse(pcum<.5, pcum, 1-pcum)
        }
#         rrrr <- cor(x)[i,j]
#         ttt <- pt(rrrr * sqrt(nrow(x)-2)/sqrt(1-rrrr^2), nrow(x)-2)
#         print(ifelse(ttt>.5, 1-ttt, ttt) )
      }
    }
    rSE.p[lower.tri(rSE.p)] <- t(rSE.p)[lower.tri(rSE.p)]
    prob[lower.tri(prob)] <- t(prob)[lower.tri(prob)]
  }


# nonparametric bootstrap
  if(bs.n>0 && "se.n" %in% bs.type){
    getBoot <- function(rep){
      bs.ests <- array(NA, dim=c(m,m,rep))
      for (i in 1: bs.n){
        sel <- sample(1:n, n,T)
        try(bs.ests[,,i] <- RRcor(x=X[sel,],models=models,p.list=p.list, group=group.2g[sel,] )$r)
      }
      bs.ests
    }
    if (nCPU == 1){
      bs.ests <- getBoot(rep=bs.n)   
    }else{
      #       require(doParallel, quietly=T)
      if (nCPU=="max"){
        try(nCPU <-  as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')))
        if (nCPU=="max") nCPU <- 2
      }
      cl.tmp =  makeCluster(nCPU) 
      registerDoParallel(cl.tmp, cores=nCPU)
      rep <- ceiling(bs.n/nCPU)
      bs.ests <- foreach(k=1:nCPU , .packages='RRreg') %dopar% {getBoot(rep) }
      bs.ests <- array(unlist(bs.ests), dim = c(m,m, rep*nCPU))
      stopCluster(cl.tmp)
    }
    #     print(apply(bs.ests,1:2,mean))
    bs.n.NA <- sum(is.na(bs.ests[2,1,]))
    rSE.n <- apply(bs.ests,1:2,sd, na.rm=T)
    diag(rSE.n) <- NA
    dimnames(rSE.n) <- dimnames(rSE.p)
  }

res <- list(r = r, call = "Randomized response correlation",
            models = models, p.list = p.list , p.char = p.char, varNames = colnames(X), n.mat=n.mat,
            n=n, bs.n=bs.n, bs.type=bs.type, rSE.p=rSE.p, rSE.n=rSE.n, prob=prob)
if (bs.n>0 && "se.n" %in% bs.type)
  res$samples.n <- bs.ests
if (bs.n>0 && "se.p" %in% bs.type)
  res$samples.p <- samples.p
class(res) <- "RRcor"
res
}


#' @export
print.RRcor<-function(x,...){
#   cat("Call: \n")
#   write(x$call,"")
  cat("Randomized response variables:\n")
  TAB <- cbind(Variable = x$varNames,
               RRmodel = x$models,
               p = x$p.char)
  rownames(TAB)=1:length(x$models)
  print(TAB, quote = FALSE)
  if (min(x$n.mat) != x$n){
    cat(paste("\nSample size N:\n"))
    print( x$n.mat)
  }else{
    cat(paste("\nSample size N =",x$n,"\n"))
  }
  cat("\nEstimated correlation matrix:\n")
  print( round(x$r,6))
  
  if(x$bs.n >0 &&  "se.p" %in% x$bs.type){
    cat("\nStandard errors from parametric bootstrap:\n")
    print( round(x$rSE.p,6), quote = FALSE) 
  }
  if(x$bs.n >0 &&  "se.n" %in% x$bs.type){
    cat("\nStandard errors from nonparametric bootstrap:\n")
    print( round(x$rSE.n,6), quote = FALSE)
  }
  if(x$bs.n >0 &&  "pval" %in% x$bs.type){
    cat("\nTwo-sided p-values from (separate) parametric bootstrap:\n")
    print( round(x$prob,6), quote = FALSE)
  }
  
}



# 
# 
# summary.RRcor<-function(object,...){
#   zval <- object$pi/object$rSE
#   TAB <- cbind(Estimate = object$r,
#                StdErr = object$rSE,
#                z.value=zval,
#                p.value=pnorm(-abs(zval)))
#   rownames(TAB) <- "r"
#   res <- list(call=object$call,n=object$n,
#               coefficients=TAB)
#   class(res) <- "summary.RRcor"
#   return(res)
# }
# 
# 
# #' @aliases RRcor
# #' @method print summary.RRcor
# #' @export
# print.summary.RRcor<-function(x,...){
#   cat("Call:\n")
#   write(x$call,"")
#   cat("Sample size: ")
#   write(x$n,"")
#   cat("\n")
#   printCoefmat(x$coefficients)
# }

