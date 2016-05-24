#' Simulating a pathway with multiple SNPs.
#'
#' It gives a simulated SNPs consisting of multiple genes in a pathway. Each SNPs from a latent multivariate Gaussian variable with an AR1 correlation structure.
#'
#' @param nGenes The number of total genes.
#'
#' @param nGenes1 The number of causal genes.
#'
#' @param nSNPs A vector, length matched with total number of genes. Each elements of vector indicate the number of SNPs in the gene. Default is nSNPs = NULL, in this case the number of nSNPs randomly selected from nSNPlow to nSNPup.
#'
#' @param ncSNPs A vector, length matched with total number of genes. Each elements of vector indicate the number of causal SNPs in the gene. Default is ncSNPs = NULL, in this case the number of ncSNPs are randomly selected from nSNP0.
#'
#' @param nSNPlim If nSNPs = NULL, the number of SNPs in Gene randomly selected from Unif(nSNPlim[1], nSNPlim[2]).
#'
#' @param nSNP0 If ncSNPs = NULL, the number of causal SNPs in Gene randomly selected from nSNP0. Default is 1:3.
#'
#' @param LOR Association in log OR between a causal SNP and outcome.
#'
#' @param n # of cases (= # of controls).
#'
#' @param MAFlim MAF's of the SNPs are drawn from Unif(MAFlim[1], MAFlim[2]).
#'
#' @param rholim the SNPs in eahc gene are from a latent Normal
#'                    variable with a AR(rho) corr structure, rho's are drawn
#'                    from Unif(rholim[1], rholim[2]); the SNPs in diff genes are independant.
#'
#' @param p0 background disease prevalence;i.e. intercept=log(p0/(1-p0)).
#'
#' @param noncausal exclude causal SNPs if TRUE, it is the simulation set up d in the paper(Pan et al 2015).
#'
#' @export
#' @return a list of the binary outcome Y (=0 or 1) and SNPs (=0, 1 or 2);
#'               Y is a vector of length 2n; X is a matrix of 2n by nSNP.
#'
#' @examples
#'
#' # Simulation set up A a) in the paper (Pan et al 2015)
#' \dontrun{ simula <- simPathAR1Snp(nGenes=20, nGenes1=1, nSNPlim=c(1, 20), nSNP0=1:3,
#'                            LOR=.2, rholim=c(0,0),
#'                            n=100, MAFlim=c(0.05, 0.4), p0=0.05) }
#' \dontshow{ simula <- simPathAR1Snp(nGenes=20, nGenes1=1, nSNPlim=c(1, 20), nSNP0=1:3,
#'                            LOR=.2, rholim=c(0,0),
#'                            n=10, MAFlim=c(0.05, 0.4), p0=0.05) }
#'
#' # Simulation set up A b) in the paper
#' #simulb <- simPathAR1Snp(nGenes=20, nGenes1=1, nSNPlim=c(1, 100), nSNP0=1:3,
#' #                           LOR=.2, rholim=c(0,0),
#' #                           n=100, MAFlim=c(0.05, 0.4), p0=0.05)
#'
#'
#'
#' @seealso \code{\link{aSPUpath}}


simPathAR1Snp<-function(nGenes=10, nGenes1=5, nSNPs=NULL, ncSNPs = NULL,
                 nSNPlim=c(1,20), nSNP0=1:3, LOR=0.3,
                 n=100, MAFlim=c(0.05, 0.4), rholim=c(0, 0),
                 p0=0.05, noncausal = FALSE){

    if (is.null(nSNPs)) nSNPs=sample(nSNPlim[1]:nSNPlim[2], nGenes, replace=T)
    if (is.null(ncSNPs)) ncSNPs=pmin(nSNPs[1:nGenes1], sample(nSNP0, nGenes1, replace=T) )
    ## total # of SNPs:
    q<-sum(nSNPs)

    rhos=runif(nGenes, rholim[1], rholim[2])

###get the overall corr matrix for the laten var's:
    R<-matrix(0, nrow=q, ncol=q)
    for(iGene in 1:nGenes){
        if(iGene==1) Rstart=0 else Rstart=sum(nSNPs[1:(iGene-1)])
        for(i in 1:nSNPs[iGene])
            for(j in 1:nSNPs[iGene])
                R[Rstart+i, Rstart+j]<-rhos[iGene]^(abs(i-j))
    }
    svd.R<-svd(R)
    R1<-svd.R$u %*% diag(sqrt(svd.R$d))
      #R2<- diag(sqrt(svd.R$d)) %*% t(svd.R$v)
      #Note above: R1 %*% t(R1)= R

##background disease prev:
    b0<-log(p0/(1-p0))
##LOR:
    b1<-rep(0, q)
    for(i in 1:nGenes1)
        if (i==1) b1[1:ncSNPs[1]]=LOR  else b1[sum(nSNPs[1:(i-1)])+(1:ncSNPs[i])]=LOR

    MAFs<-runif(q, MAFlim[1], MAFlim[2])
    cutoff<-qnorm(MAFs)

    X<-matrix(0, nrow=n+n, ncol=q)
    Y<-rep(0, n+n); Y[(n+1):(2*n)]<-1
    i<-1
#sampling controls:
    while ( i <= n){
        X0<-rnorm(q, 0, 1) #: X0 ~ MVN(0, I)
        X1<-R1 %*% X0   #: X1 ~ MVN(0, R)
        X2<-ifelse(X1<cutoff, 1, 0)
        X0<-rnorm(q, 0, 1) #: X0 ~ MVN(0, I)
        X1<-R1 %*% X0   #: X1 ~ MVN(0, R)
        X3<-ifelse(X1<cutoff, 1, 0)
        X4<-X2+ X3
        pr<-1/(1 + exp(-(b0 + sum(b1 * X4))))
        Y1<-sample(c(0, 1), 1, prob=c(1-pr, pr))
        if (Y1==0){
            X[i, ]<-X4
            i<-i+1
        }
    }
#sampling cases:
    while ( i <= 2*n){
        X0<-rnorm(q, 0, 1) #: X0 ~ MVN(0, I)
        X1<-R1 %*% X0   #: X1 ~ MVN(0, R)
        X2<-ifelse(X1<cutoff, 1, 0)
        X0<-rnorm(q, 0, 1) #: X0 ~ MVN(0, I)
        X1<-R1 %*% X0   #: X1 ~ MVN(0, R)
        X3<-ifelse(X1<cutoff, 1, 0)
        X4<-X2+ X3
        pr<-1/(1 + exp(-(b0 + sum(b1 * X4))))
        Y1<-sample(c(0, 1), 1, prob=c(1-pr, pr))
        if (Y1==1){
            X[i, ]<-X4
            i<-i+1
        }
    }

    if ( noncausal == FALSE ) {
        snp.chrom=snp.loc=NULL
        for(i in 1:nGenes){
            snp.chrom=c(snp.chrom, rep(i, nSNPs[i]))
            snp.loc=c(snp.loc, 1:nSNPs[i])
        }
        snp.info = cbind(1:q, snp.chrom, snp.loc)
        gene.info=cbind(1:nGenes, 1:nGenes, rep(0, nGenes), nSNPs)
        pathway = list(pathway1=as.character(1:nGenes))

        return( list(Y=Y, X=X, snp.info=snp.info, gene.info=gene.info,
                     pathway=pathway, nSNPs=nSNPs) )
    } else { ## exclude causal SNPs
        X = X[, abs(b1)<1e-10]

        snp.chrom=snp.loc=NULL
        for(i in 1:nGenes1){
            snp.chrom=c(snp.chrom, rep(i, nSNPs[i]-ncSNPs[i]))
            if (nSNPs[i]-ncSNPs[i]>=1)
                snp.loc=c(snp.loc, 1:(nSNPs[i]-ncSNPs[i]))
        }
        if (nGenes1 < nGenes){
            for(i in (nGenes1+1):nGenes){
                snp.chrom=c(snp.chrom, rep(i, nSNPs[i]))
                snp.loc=c(snp.loc, 1:nSNPs[i])
            }
        }
        snp.info = cbind(1:(q-sum(ncSNPs)), snp.chrom, snp.loc)
        nSNPsNone0=nSNPs
        for(i in 1:nGenes1)
            nSNPsNone0[i] = nSNPs[i]-ncSNPs[i]
        gene.info=cbind(1:nGenes, 1:nGenes, rep(0, nGenes), nSNPsNone0)
        pathway = list(pathway1=as.character(1:nGenes))

        return ( list(Y=Y, X=X, snp.info=snp.info, gene.info=gene.info,
                      pathway=pathway, nSNPs=nSNPsNone0) )
    }


}
