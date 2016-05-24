#' Pathway based Sum of Powered Score tests (SPUpath) and adaptive SPUpath (aSPUpath) test for single trait - pathway association.
#'
#' It gives p-values of the SPUpath tests and aSPUpath test.
#'
#' @param Y Response or phenotype data. It can be a disease indicator; =0 for controls, =1 for cases.
#' Or it can be a quantitative trait. A vector with length n (number of observations).
#'
#' @param X Genotype or other data; each row for a subject, and each column
#'     for an SNP (or a predictor). The value of each SNP is the # of the copies
#'     for an allele. A matrix with dimension n by k (n : number of observation, k : number of SNPs (or predictors) ).
#'
#' @param cov Covariates. A matrix with dimension n by p (n :number of observation, p : number of covariates).
#'
#' @param model Use "gaussian" for a quantitative trait, and use "binomial" for a binary trait.
#'
#' @param snp.info SNP information matrix, the 1st column is SNP id, 2nd column is chromosome #, 3rd column indicates SNP location.
#'
#' @param gene.info GENE information matrix, The 1st column is GENE id, 2nd column is chromosome #, 3rd and 4th column indicate start and end positions of the gene.
#'
#' @param pow SNP specific power(gamma values) used in SPUpath test.
#'
#' @param pow2 GENE specific power(gamma values) used in SPUpath test.
#'
#' @param n.perm number of permutations.
#'
#' @param usePCs indicating whether to extract PCs and then use PCs of X.
#'
#' @param varprop the proportion of the variations explained (cutoff) that
#'                 determines how many top PCs to use.
#'
#' @export
#' @return P-values for SPUpath tests and aSPUpath test.
#'
#' @author Il-Youp Kwak and Wei Pan
#'
#' @references
#' Wei Pan, Il-Youp Kwak and Peng Wei (2015)
#' A Powerful and Pathway-Based Adaptive Test for Genetic Association With Common or Rare Variants, The American Journal of Human Genetics, 97, 86-98
#'
#' @examples
#'
#' \dontrun{dat1<-simPathAR1Snp(nGenes=20, nGenes1=5, nSNPlim=c(1, 20),
#' 	       nSNP0=1, LOR=.2, n=100, MAFlim=c(0.05, 0.4), p0=0.05 ) }
#' \dontshow{dat1<-simPathAR1Snp(nGenes=20, nGenes1=5, nSNPlim=c(1, 20),
#'             nSNP0=1, LOR=.2, n=40, MAFlim=c(0.05, 0.4), p0=0.05 ) }
#'
#' # p-values of SPUpath and aSPUpath tests.
#' \dontrun{p.pathaspu<- aSPUpath(dat1$Y, dat1$X, snp.info = dat1$snp.info,
#'          gene.info = dat1$gene.info,
#'          model = "binomial", pow=1:8, pow2=c(1, 2, 4, 8), n.perm=1000) }
#' \dontshow{p.pathaspu<- aSPUpath(dat1$Y, dat1$X, snp.info = dat1$snp.info,
#'          gene.info = dat1$gene.info,
#'          model = "binomial", pow=1:8, pow2=c(1, 2, 4, 8), n.perm=30) }
#' p.pathaspu
#' ## pow = 1:8 and pow2 = 1,2,4,8
#' ## So, there are 8*4 = 32 SPUpath p-values.
#' ## SPUpathi,j corresponds pow = i , pow2 = j
#' ## The last element, aSPUpath gives aSPUpath p-value.
#'
#' @seealso \code{\link{simPathAR1Snp}} 

aSPUpath <- function(Y, X, cov = NULL, model=c("binomial", "gaussian"),
                     snp.info, gene.info, pow=c(1:8, Inf), pow2=c(1,2,4,8),
                     n.perm=200, usePCs=F, varprop=0.95 ){

    model = match.arg(model)

    n.gene <- nrow(gene.info)
    GL <- list(0)
    i = 1
    for(g in 1:n.gene) { # g = 2
        snpTF <- ( snp.info[,2] == gene.info[g,2] &
                      gene.info[g,3] <= as.numeric(snp.info[,3]) &
                          gene.info[g,4] >= as.numeric(snp.info[,3]) )

        if( sum(snpTF) != 0){
            GL[[i]] <- which(snpTF)
            i = i + 1
        }
    }

    X = X[, unlist(GL)]
    nSNPs=unlist(lapply(GL,length))

    nSNPs0<-rep(0, length(nSNPs))
    if (usePCs){
        Xg<-NULL
        for(iGene in 1:length(nSNPs)){
            if (iGene==1) SNPstart=1 else SNPstart=sum(nSNPs[1:(iGene-1)])+1
            indx=(SNPstart:(SNPstart+nSNPs[iGene]-1))
            Xpcs<-extractPCs(X[, indx], cutoff=varprop)
            Xg<-cbind(Xg, Xpcs)
            if (is.null(ncol(Xpcs))) nSNPs0[iGene]=1
            else nSNPs0[iGene]=ncol(Xpcs)
        }
    } else { Xg=X; nSNPs0=nSNPs}

    n<-length(Y)
    k<-ncol(Xg)

#######construction of the score vector:
    if (is.null(cov)){
        ## NO nuisance parameters:
        XUs<-Xg
        r<-Y-mean(Y)
        U<-as.vector(t(Xg) %*% r)
    } else {
        tdat1<-data.frame(trait=Y, cov)
        fit1<-glm(trait~.,family=model,data=tdat1)
        pis<-fitted.values(fit1)
        XUs<-matrix(0, nrow=n, ncol=k)
        Xmus = Xg
        for(i in 1:k){
            tdat2<-data.frame(X1=X[,i], cov)
            fit2<-glm(X1~.,data=tdat2)
            Xmus[,i]<-fitted.values(fit2)
            XUs[, i]<-(X[,i] - Xmus[,i])
        }
        r<-Y - pis
        U<-t(XUs) %*% r
    }

    nGenes=length(nSNPs0)

                                        # test stat's:
    TsUnnorm<-Ts<-StdTs<-rep(0, length(pow)*nGenes)
    for(j in 1:length(pow))
        for(iGene in 1:nGenes){
            if (iGene==1) SNPstart=1 else SNPstart=sum(nSNPs0[1:(iGene-1)])+1
            indx=(SNPstart:(SNPstart+nSNPs0[iGene]-1))
            if (pow[j] < Inf){
                a= (sum(U[indx]^pow[j]))
                TsUnnorm[(j-1)*nGenes+iGene] = a
                Ts[(j-1)*nGenes+iGene] = sign(a)*((abs(a)) ^(1/pow[j]))
                StdTs[(j-1)*nGenes+iGene] = sign(a)*((abs(a)/nSNPs0[iGene]) ^(1/pow[j]))
        # (-1)^(1/3)=NaN!
        #Ts[(j-1)*nGenes+iGene] = (sum(U[indx]^pow[j]))^(1/pow[j])
            }
            else {
                TsUnnorm[(j-1)*nGenes+iGene] = Ts[(j-1)*nGenes+iGene] = StdTs[(j-1)*nGenes+iGene] =max(abs(U[indx]))
            }
        }

                                        # Permutations:
    T0sUnnorm=T0s = StdT0s = matrix(0, nrow=n.perm, ncol=length(pow)*nGenes)
    for(b in 1:n.perm){
        R0 <- sample(r, length(r))
#########Null score vector:
        U0<-t(XUs) %*% R0

                                        # test stat's:
        for(j in 1:length(pow))
            for(iGene in 1:nGenes){
                if (iGene==1) SNPstart=1 else SNPstart=sum(nSNPs0[1:(iGene-1)])+1
                indx=(SNPstart:(SNPstart+nSNPs0[iGene]-1))
                if (pow[j] < Inf){
                    a = (sum(U0[indx]^pow[j]))

#                    T0sUnnorm[b, (j-1)*nGenes+iGene] = a

#                    T0s[b, (j-1)*nGenes+iGene] = sign(a)*((abs(a)) ^(1/pow[j]))
                    StdT0s[b, (j-1)*nGenes+iGene] = sign(a)*((abs(a)/nSNPs0[iGene]) ^(1/pow[j]))
                }

                else StdT0s[b, (j-1)*nGenes+iGene] = max(abs(U0[indx]))
            }
    }

   #combine gene-level stats to obtain pathway-lelev stats:
    Ts2<-rep(0, length(pow)*length(pow2))
    T0s2<-matrix(0, nrow=n.perm, ncol=length(pow)*length(pow2))
    for(j2 in 1:length(pow2)){
        for(j in 1:length(pow)){
#            TsU[(j2-1)*length(pow) +j] = sum(TsUnnorm[((j-1)*nGenes+1):(j*nGenes)]^pow2[j2])
#            Ts1[(j2-1)*length(pow) +j] = sum(Ts[((j-1)*nGenes+1):(j*nGenes)]^pow2[j2])
            Ts2[(j2-1)*length(pow) +j] = sum(StdTs[((j-1)*nGenes+1):(j*nGenes)]^pow2[j2])
            for(b in 1:n.perm){
#                T0sU[b, (j2-1)*length(pow) +j] = sum(T0sUnnorm[b, ((j-1)*nGenes+1):(j*nGenes)]^pow2[j2])
#                T0s1[b, (j2-1)*length(pow) +j] = sum(T0s[b, ((j-1)*nGenes+1):(j*nGenes)]^pow2[j2])
                T0s2[b, (j2-1)*length(pow) +j] = sum(StdT0s[b, ((j-1)*nGenes+1):(j*nGenes)]^pow2[j2])
            }
        }
    }

   # permutation-based p-values:
    pPerm2 = rep(NA, length(pow)*length(pow2));
    pvs = NULL;

    for(j in 1:(length(pow)*length(pow2))) {
#        pPermU[j] = sum( abs(TsU[j]) < abs(T0sU[,j]))/n.perm
#        pPerm1[j] = sum( abs(Ts1[j]) < abs(T0s1[,j]))/n.perm
        pPerm2[j] = sum( abs(Ts2[j]) < abs(T0s2[,j]))/n.perm
    }

#    P0s1 = PermPvs(T0s1)
    P0s2 = PermPvs(T0s2)
#    P0sU = PermPvs(T0sU)
#    minP0s1 = apply(P0s1, 1, min)
    minP0s2 = apply(P0s2, 1, min)
#    minP0sU = apply(P0sU, 1, min)
#    minP1 =  sum( min(pPerm1) > minP0s1 )/n.perm
    minP2 =  sum( min(pPerm2) > minP0s2 )/n.perm
#    minPU =  sum( min(pPermU) > minP0sU )/n.perm
    minP1s<-minP2s<-minPUs<-rep(NA, length(pow2))
#    for(j2 in 1:length(pow2)){
#        minP0s1 = apply(P0s1[, ((j2-1)*length(pow)+1):(j2*length(pow))] , 1, min)
#        minP0s2 = apply(P0s2[, ((j2-1)*length(pow)+1):(j2*length(pow))], 1, min)
#        minP0sU = apply(P0sU[, ((j2-1)*length(pow)+1):(j2*length(pow))], 1, min)
#        minP1s[j2] =  sum( min(pPerm1[((j2-1)*length(pow)+1):(j2*length(pow))]) > minP0s1 )/n.perm
#        minP2s[j2] =  sum( min(pPerm2[((j2-1)*length(pow)+1):(j2*length(pow))]) > minP0s2 )/n.perm
#        minPUs[j2] =  sum( min(pPermU[((j2-1)*length(pow)+1):(j2*length(pow))]) > minP0sU )/n.perm
#    }
#    stdPs=c(pPerm2, minP2s, minP2)
    pvs=c(pPerm2, minP2)

    nmvec <- NULL;
    for(ii in pow2) {
    	   for(jj in pow) {
           	  nmvec <- c(nmvec, paste("SPUMpath",jj,",",ii,sep=""))
           }
    }

#    nmvec <- NULL;
#    for(nm in paste("SPUpath",pow,",", sep=""))
#        nmvec <- c(nmvec, paste(nm, pow2, sep="") )

    nmvec <- c(nmvec, "aSPUpath")
    names(pvs) <- nmvec
    pvs
}




