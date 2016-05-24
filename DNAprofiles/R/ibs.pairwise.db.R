#' Pairwise comparison of all database profiles on IBS alleles
#' 
#' Compares every database profile with every other database profile and keeps track of the number of pairs that match fully and partially on all numbers of loci.
#' @param db An integer matrix which is the database of profiles.
#' @param hit Integer; when > 0, the function keeps track of the pairs with at least this number of matching loci
#' @param showprogress Logical; show progress bar? (not available when \code{multicore=TRUE})
#' @param multicore Logical; use multicore implementation?
#' @param ncores Integer value, with \code{multicore=TRUE}, the number of cores to use or 0 for auto-detect.
#' @details Makes all pairwise comparisons of profiles in \code{db}. Counts the number of profiles that match fully/partially for each number of loci.
#' 
#' The number of pairwise comparisons equals \eqn{N*(N-1)/2}, where \eqn{N} equals the number of database profiles, so the computation time grows quadratically in \eqn{N}. The procedure using a single core takes a few minutes applied to a database of size 100.000 (Intel I5@@2.5GHz), but the time quadruples each time the database becomes twice as large.
#' 
#' A similar function with additional functionality is available in the \code{DNAtools} package. That function however does not handle large databases (about 70k is the maximum) and is a few times slower than the implementation used here. The \code{DNAtools} package comes with a specialized plotting function that can be used with the output of the \code{db.compare.pairwise} function after converting with \code{\link{as.dbcompare}}.
#' @return Matrix with the number of full/partial matches on 0,1,2,... loci.
#' @seealso \code{\link{as.dbcompare}}
#' @examples
#' data(freqsNLsgmplus)
#' 
#' # sample small db and make all pairwise comparisons
#' db <- sample.profiles(N=10^3,freqs=freqsNLsgmplus)
#' ibs.pairwise.db(db)
#' 
#' \dontrun{
#' # the multicore function has some overhead and is not faster when applied to small databases
#' db.small <- sample.profiles(N=10^4,freqs=freqsNLsgmplus)
#' 
#' system.time(Msingle <- ibs.pairwise.db(db.small))
#' system.time(Mmulti <- ibs.pairwise.db(db.small,multicore=T))
#' 
#' all.equal(Msingle,Mmulti)
#' 
#' # but significant speed gains are seen for large databases (46 vs 23 secs on my system)
#' 
#' db.large <- sample.profiles(N=5*10^4,freqs=freqsNLsgmplus)
#' 
#' system.time(Msingle <- ibs.pairwise.db(db.large))
#' system.time(Mmulti <- ibs.pairwise.db(db.large,multicore=T))
#' 
#' all.equal(Msingle,Mmulti)
#' }
#' @export
ibs.pairwise.db <- function(db,hit=0,showprogress=TRUE,multicore=FALSE,ncores=0){
  nloci <- ncol(db)/2
  db.t <- as.vector(t(db)) #the c++ function expects a vector (transpose is probably more efficient due to caching)
  
  if (!multicore){
    #single core function
    #actual work is done by (not exported) c++ function
    if (hit==0){
      ret <- list(M=Zdbcomparepairwise(db.t,nloci,showprogress))
    } else{
      ret <- Zdbcomparepairwisetrackhits(db.t,nloci,hit,showprogress)
    }
  }else{
    #multicore function
    #require(parallel)
    if (ncores==0) ncores <- detectCores()
    cl <- makeCluster(ncores)
    clusterExport(cl, c("ncores","nloci","db.t","hit","Zdbcomparepairwisemc","Zdbcomparepairwisemctrackhits"),envir=environment())
    if (hit==0){
      r <- parLapply(cl,1:ncores,function(j) Zdbcomparepairwisemc(db.t,nloci,njobs=ncores,job=j))
      ret <- list(M=Reduce("+",r))
    }else{
      r <- parLapply(cl,1:ncores,function(j) Zdbcomparepairwisemctrackhits(db.t,nloci,hit,njobs=ncores,job=j))
      ret <- list(M=Reduce("+",lapply(r,function(y) y$M)),
                  hits=Reduce("rbind",lapply(r,function(y) y$hits)))    
    }
    #ret <- Reduce("+",parLapply(cl,1:ncores,function(j) Zdbcomparepairwisemc(db.t,nloci,njobs=ncores,job=j)))
    stopCluster(cl)
  }
  
  dimnames(ret$M) <- list(match=as.character(0:nloci),partial=as.character(0:nloci))
  if (!is.null(ret$hits)){
    if (nrow(ret$hits)>0){
      #sort by full matches, then partial
      ret$hits <- ret$hits[order((nloci+1)*ret$hits[,3]+ret$hits[,4],decreasing=T),] 
    }else ret <- list(M=ret$M)  
  }else{
    ret <- list(M=ret$M)  
  }
  ret
}
NULL
#' Compute the probabilities that two profiles match a number of loci fully/partially
#' 
#' When two profiles are compared, this function computes the expected number of loci that match fully or partially. 
#' @param freqs1 List of allelic frequencies.
#' @param freqs2 List of allelic frequencies.
#' @param k IBD-probabilities, passed on to \code{\link{ibdprobs}}. Defaults to "UN", i.e. unrelated.
#' @details When all profiles in the database are compared pairwise, one can count the number of profiles that match fully/partially for each number of loci. Such a procedure is implemented as \code{\link{ibs.pairwise.db}}. The current function computes the probabilities that a single pair matches fully and partially at each possible number of loci. The two profiles can be assumed to originate from population with different allele frequencies (\code{freqs1} and \code{freqs2}), might be related through and might both be inbred.
#' 
#' @return Matrix with the expected number of full/partial matches on 0,1,2,... loci for a comparison between two profiles.
#' @seealso \code{\link{as.dbcompare}}
#' @examples
#' data(freqsNLsgmplus)
#' 
#' # sample small db and make all pairwise comparisons
#' 
#' N <- 1e3
#' db <- sample.profiles(N=N,freqs=freqsNLsgmplus)
#' 
#' O <- ibs.pairwise.db(db)
#' E <- N*(N-1)/2*ibs.pairwise.pr(freqs1 = freqsNLsgmplus,freqs2 = freqsNLsgmplus)
#' 
#' O # observed
#' E # expected
#' @export
ibs.pairwise.pr <- function(freqs1,freqs2=freqs1,k="UN"){  
  if (identical(names(freqs1),names(freqs2))&&
        identical(sapply(freqs1,length),sapply(freqs2,length))){
    # a single set of allelic freqs is supplied
    # compute the prob of 0,1,2 ibs alleles for all loci
    M.012 <- lapply(names(freqs1),function(L){
      Zibs.pairwise.pr.locus(freqs1[[L]],freqs2[[L]],ibdprobs(k))
    })
    # compute the matrix of full/partial match probabilities from the last of probs of 0,1,2 ibs
    M <-  ZMexp(M.012)
  }else{
    stop("Please supply proper allelic frequencies")
  }
  M
}
NULL
#' Compute expected number of profiles pairs in a database that match fully and partially at each number of loci
#' 
#' For a database comparison exercise, this function computes the expected number of pairs that match fully or partially at each number of loci in a heterogeneous database.
#' @param subpops List with allele frequencies in each subpopulation.
#' @param fractions Numeric 
#' @param N Total size of database.
#' @param ks IBD-probabilities for the relations that occur within subpopulations (with probabilities \code{alpha.w}) and between (with probabilities \code{alpha.b}). Passed on to \code{\link{ibdprobs}}.
#' @param alpha.w Numeric with same length as \code{ks}. The \code{i}'th element denotes the probability that a pair of profiles within a subpopulation is related as described by the \code{i}'th element of \code{ks}.
#' @details When all profiles in the database are compared pairwise, one can count the number of profiles that match fully/partially for each number of loci. Such a procedure is implemented as \code{\link{ibs.pairwise.db}}. The current function computes the expected value of this matrix. The database can be heterogeneous (consisting of subpopulations with different allele frequencies) and within-subpopulation inbreeding is supported.
#' 
#' @return Matrix with the expected number of full/partial matches on 0,1,2,... loci in the database.
#' @seealso \code{\link{as.dbcompare}}
#' @examples
#' data(freqsNLsgmplus)
#' 
#' # sample small db, make all pairwise comparisons and compute the expected number
#' 
#' N <- 1e3
#' db <- sample.profiles(N=N,freqs=freqsNLsgmplus)
#' 
#' O <- ibs.pairwise.db(db)
#' E <- ibs.pairwise.db.exp(subpops = list(freqsNLsgmplus),N = N)
#' 
#' O # observed
#' E # expected
#' @export
ibs.pairwise.db.exp <- function(subpops,fractions=rep(1/length(subpops),length(subpops)),N=2L,
                                                              ks=c("UN","FS","PO"),alpha.w = c(1,0,0)){

  alpha.b <- as.numeric(sapply(ks,function(k) identical(ibdprobs(k),c(1,0,0)))) # between -> unrelated
  
  
  if (length(fractions)!=length(subpops)) stop("The number of subpopulations is not equal to the number of subpopulation fractions")
  if (any((fractions<0)|(fractions>1))) stop("Fractions do not all lie between 0 and 1")    
  if ((round(sum(fractions),digits = 8)-1)!=0) stop("Fractions do not sum to one")
  
    if (!all.equal(length(alpha.w), length(alpha.b), length(ks))) stop("ks, alpha.w and alpha.b should be of equal length")
  if (((round(sum(alpha.w),digits = 8)-1)!=0)|((round(sum(alpha.b),digits = 8)-1)!=0)|any(alpha.w>1)|any(alpha.w<0)|any(alpha.b>1)|any(alpha.b<0)) warning("Implausible values for alpha.w or alpha.b: ",alpha.w,alpha.b)

  K <- length(subpops)
  p <- fractions
  
  comp <- Zcomb.pairs(nn = K) # enumerate pairwise subpop comparisons
  n0 <- N*p # no. profiles in each subpop

  comp.n <- ifelse(comp[,1]==comp[,2],
                   yes = choose(n0[comp[,1]],2), # no. within subpop comps
                   no = n0[comp[,1]]*n0[comp[,2]]) # no. between subpop comps
  
  
  # helper function to compute expectation for each comb. of subpops and kappa
  comp.exps <- function(comp,subpops,ks){
    ret <- list()
    for(j in seq_len(nrow(comp))){
      ret[[j]] <- list()
      k1 <- comp[j,1]; k2 <- comp[j,2] # compare subpop k1 with k2
      
      ret[[j]] <- lapply(ks, function(k)
        ibs.pairwise.pr(subpops[[k1]],subpops[[k2]],k=k))
    }
    ret
  }
  
  # computes grand total using expecteds and alphas (frac of comparisons according to kappas)
  total.exps <- function(a.w,a.b,expecteds,comp,comp.n){
    tmp <- list()
    for(j in seq_len(nrow(comp))){    
      k1 <- comp[j,1]; k2 <- comp[j,2]
      
      if (k1==k2){
        # within pop
        tmp[[j]] <- Reduce("+",mapply("*", expecteds[[j]], a.w[,k1]*comp.n[j], SIMPLIFY = FALSE))
      }else{
        #between pops
        tmp[[j]] <- Reduce("+",mapply("*", expecteds[[j]], a.b*comp.n[j], SIMPLIFY = FALSE))
      } 
    }
    Reduce("+",tmp)
  }
  
  # obtain full matrix for each comb of subpops and kappa
  expecteds <- comp.exps(comp = comp,subpops = subpops,ks = ks)
  # sum with alpha factors to obtain a grand total
  total.exps(a.w = matrix(rep(alpha.w,K),ncol=K),a.b = alpha.b,expecteds = expecteds,comp = comp,comp.n)
}
NULL
ZMexp <- function(M){
  # takes a list of 0,1,2 ibs probabilities and gives back the matrix with pr's of partial and full matches
  ret <- matrix(1,nrow=1,ncol=1)
  #  recursively compute the M matrix
  for (l in seq_along(M)){
    p <- M[[l]] # pr. of 0,1,2 ibs @ locus
    n <- nrow(ret) # size of the matrix so far
    ret <- cbind(rbind(ret*p[1],rep(0,n)),0) +
      cbind(rep(0,n+1),rbind(ret*p[2],rep(0,n))) +
      rbind(rep(0,n+1), cbind(ret*p[3],0))
  }
  dimnames(ret) <- list(match=as.character(0:length(M)),partial=as.character(0:length(M)))
  ret    
}
NULL
Zibs.pairwise.pr.locus <- function(p,q,ibdp=c(1,0,0)){
  # computes the pr. that two persons (x,y) have 0 or 2 alleles ibs
  # from which pr. for 1 IBS follows
  # p are allele freqs of x, q of y
  
  # f1 is inbreeding coefficient of x (Wrights Fis); f2 of y
  f1=0
  f2=0 # not used for now...
  
  k <- ibdp
  A <- length(p)
  # case that x is homozygous: aa
  aa <- seq_along(p) # genos
  x.aa.pr <- p[aa]^2*(1-f1)+p[aa]*f1
  # pr that allele from person x is not a
  y.nota.pr <- 1-q[aa]
  # pr. that person x and y share 0 IBS times indicator x is aa
  aa.p0 <- sum(x.aa.pr*k[1]*(y.nota.pr^2*(1-f2)+y.nota.pr*f2))
  y.aa.pr <- q[aa]^2*(1-f2)+q[aa]*f2
  y.a.pr <- q[aa]
  # pr that x,y share 2 IBS times indicator x is aa
  aa.p2 <- sum(x.aa.pr*(k[1]*y.aa.pr+k[2]*y.a.pr+k[3]))
  
  # case that x is heterozygous: ab
  ab <- t(combn(A,m = 2))
  x.ab.pr <- 2*p[ab[,1]]*p[ab[,2]]*(1-f1)
  y.ab.pr <- 2*q[ab[,1]]*q[ab[,2]]*(1-f2)
  y.a.pr <- q[ab[,1]]
  y.b.pr <- q[ab[,2]]
  # pr. that allele from person y is not a and not b
  y.notab.pr <- 1-q[ab[,1]]-q[ab[,2]]
  # pr. that person x and y share 0 IBS times indicator x is ab
  ab.p0 <- sum(x.ab.pr*k[1]*(y.notab.pr*f2+y.notab.pr^2*(1-f2)))
  # pr. that x,y share 2 IBS times indicator x is ab
  ab.p2 <- sum(x.ab.pr*((k[1]*y.ab.pr)+(k[2]*(0.5*(y.a.pr+y.b.pr)))+k[3]))
  p0 <- aa.p0+ab.p0
  p2 <- aa.p2+ab.p2
  c(p0,max(1-p0-p2,0),p2)  
}