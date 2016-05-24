#' Multitraits - Pathway based Sum of Powered Score tests (MTSPUsSetPath) and adaptive MTSPUsSetpath (MTaSPUsSetPath) test for multi trait - pathway association with GWAS summary statistics.
#'
#' It gives p-values of the MTSPUsSetPath tests and MTaSPUsSetPath test with GWAS summary statistics.
#'
#' @param Zs Z-score matrix. row represent SNPs and column represent traits. It could be P-values if the Ps option is TRUE. 
#'
#' @param corSNP Correlation matirx of the SNPs to be tested; estimated from a
#' reference panel (based on the same set of the reference alleles as
#' used in calculating Z-scores).
#'
#' @param corPhe Correlation matirx of phenotypes to be tested; Estimated from Z-scores.
#' 
#' @param snp.info SNP information matrix, the 1st column is SNP id, 2nd column is chromosome #, 3rd column indicates SNP location.
#'
#' @param gene.info GENE information matrix, The 1st column is GENE id, 2nd column is chromosome #, 3rd and 4th column indicate start and end positions of the gene.
#'
#' @param pow1 SNP specific power(gamma values) used in MTSPUsSetPath test.
#'
#' @param pow2 GENE specific power(gamma values) used in MTSPUsSetPath test.
#'
#' @param pow3 Trait specific power(gamma values) used in MTSPUsSetPath test.
#'
#' @param n.perm number of permutations.
#'
#' @param Ps TRUE if input is p-value, FALSE if input is Z-scores. The default is FALSE.
#'
#' @export
#' @return P-values for MTSPUsSetpath tests and MTaSPUsSetpPath test.
#'
#' @author Il-Youp Kwak and Wei Pan
#'
#' @references
#' Il-Youp Kwak, Wei Pan (2016)
#' Gene- and pathway-based association tests for multiple
#'       traits with GWAS summary statistics
#'
#' @examples
#'
#' Zs <- cbind ( c( 0.3, 0.2, 0.1,0.5,1.2),
#'           c(-.1, .3,-.1,.1,1.2) )
#'
#' varSNP = cbind( c( 1, .1,0, 0, .11),
#'                c(.1,1, 0, 0, 0),
#'                c(0,0, 1, 0, 0),
#'                c(0,0, 0, 1, 0),
#'                c(.11,0,0,0,1) )
#'
#' varPhe = cbind( c( 1, -.1),
#'                c(-.1,1) )
#'
#' gene.info = data.frame( Gnm = c( "G1", "G2"), chr = c(1,3),
#'                     loc1 = c(0, 0), loc2 = c(10,10) )
#' snp.info = data.frame( rsid = c("rs1", "rs2", "rs3", "rs4", "rs5"),
#'                 chr = c(1,1,3,3,3), loc = c(1,2,1,2,3) )
#'  out <- MTaSPUsSetPath(Zs, corPhe = varPhe, corSNP=varSNP,
#'             n.perm = 100, snp.info = snp.info, gene.info = gene.info)
#' out
#'
#' @seealso \code{\link{MTaSPUsSetC}}, \code{\link{MTaSPUsSet}}


MTaSPUsSetPath <- function(Zs, corPhe, corSNP, pow1=c(1,2,4,8),
                           pow2 = c(1,2,4,8),
                           pow3 = c(1,2,4,8),
                           snp.info, gene.info, n.perm=1000,
                           Ps = FALSE) {

    nsnp <- dim(Zs)[1]
    nphe <- dim(Zs)[2]
    nGenes <- nrow(gene.info)
    GL <- list(0)
    GLch <- NULL
                                        #   GL.CovSsqrt <- list(0)
    i = 1
    for(g in 1:nGenes) { # g = 1
        snpTF <- ( snp.info[,2] == gene.info[g,2] &
                   gene.info[g,3] <= as.numeric(snp.info[,3]) &
                   gene.info[g,4] >= as.numeric(snp.info[,3]) )
        
        if( sum(snpTF) != 0){
            GL[[i]] <- which(snpTF)
            GLch <- c(GLch, gene.info[g,2])
            i = i + 1
        }
    }
    
    chrs <- unique(GLch)
    CH <- list(0)
    CH.CovSsqrt <- list(0)
    for( i in 1:length(chrs) ) { # i = 2
        c = chrs[i]
        CH[[i]] <- unlist(GL[which( GLch == c )])
        Covtemp <- corSNP[CH[[i]], CH[[i]]]
        eS <- eigen(Covtemp, symmetric = TRUE)
        ev <- eS$values
        k1 <- length(ev)
        CH.CovSsqrt[[i]] <- diag(sqrt(pmax(ev, 0)), k1) %*% t(eS$vectors)
    }

    Zs = Zs[unlist(GL),]
    nSNPs0=unlist(lapply(GL,length))


    U <- corPhe

    e.U <- eigen(U)
    e.U$values[ e.U$values < 0 ] = 0
    A <- e.U$vectors %*% diag(sqrt(e.U$values))

    if(Ps == TRUE)
        Zs <- qnorm(1 - Zs/2)

    ## Permutations:
    StdTs <- list(0)
    StdT0s <- list(0)
    for(np in 1:nphe) {
        StdTs[[np]] <- rep(0, length(pow1)*nGenes)
        StdT0s[[np]] = matrix(0, nrow=n.perm, ncol=length(pow1)*nGenes)
    }

    for(b in 1:n.perm){
        
        ## U00<-rnorm(nsnp, 0, 1)
        U00 <- matrix(rnorm(n=nsnp*nphe, 0,1), nphe, nsnp)
        U0 <- NULL;
        for( ss in 1:length(CH)) { # ss = 21
            U0 <- rbind(U0, t(A %*% U00[,CH[[ss]]] %*% CH.CovSsqrt[[ss]])  )
        }
        if(Ps == TRUE)
            U0 <- abs(U0)

        ## test stat's: SPUs
        for(np in 1:nphe) {
            for(j in 1:length(pow1))
                for(iGene in 1:nGenes){
                    if (iGene==1) SNPstart=1 else SNPstart=sum(nSNPs0[1:(iGene-1)])+1
                    indx=(SNPstart:(SNPstart+nSNPs0[iGene]-1))
                    if (pow1[j] < Inf){
                        a = (sum(U0[indx,np]^pow1[j]))
                        StdT0s[[np]][b, (j-1)*nGenes+iGene] = sign(a)*((abs(a)/nSNPs0[iGene]) ^(1/pow1[j]))
                        aa = (sum(Zs[indx,np]^pow1[j]))
                        StdTs[[np]][(j-1)*nGenes+iGene] = sign(aa)*((abs(aa)/nSNPs0[iGene]) ^(1/pow1[j]))

                    } else StdT0s[[np]][b, (j-1)*nGenes+iGene] = max(abs(U0[indx,np]))
                }
        }        
    }

    ##combine gene-level stats to obtain pathway-lelev stats:

    Ts2 <- list(0)
    T0s2 <- list(0)

    for(np in 1:nphe) {
        Ts2[[np]] <- rep(0, length(pow1)*length(pow2))
        T0s2[[np]] <- matrix(0, nrow=n.perm, ncol=length(pow1)*length(pow2))

        for(j2 in 1:length(pow2)){
            for(j in 1:length(pow1)){
                a = sum(StdTs[[np]][((j-1)*nGenes+1):(j*nGenes)]^pow2[j2])
                Ts2[[np]][(j2-1)*length(pow1) +j] = sign(a)*((abs(a)/nGenes) ^(1/pow2[j2]))
                for(b in 1:n.perm){
                    aa = sum(StdT0s[[np]][b, ((j-1)*nGenes+1):(j*nGenes)]^pow2[j2])
                    T0s2[[np]][b, (j2-1)*length(pow1) +j] = sign(aa)*((abs(aa)/nGenes) ^(1/pow2[j2]))
                }
            }
        }
    }
    
    Ts3 <- rep(0, length(pow1)*length(pow2)*length(pow3) )
    T0s3 <- matrix(0, nrow=n.perm, ncol=length(pow1)*length(pow2)*length(pow3) )

    for( p3 in 1:length(pow3) ) {
        for( p2 in 1:length(pow2) ) {
            for( p1 in 1:length(pow1) ) {
                aa <- 0
                for(np in 1:nphe) {
                    aa = aa + Ts2[[np]][(p2-1)*length(pow1) + p1] ^ pow3[p3]
                }
                Ts3[ (p3-1)*length(pow1)*length(pow2) + (p2-1)*length(pow1) + p1] = aa

                for(b in 1:n.perm) {
                    aaa <- 0
                    for(np in 1:nphe) {
                        aaa = aaa + T0s2[[np]][b,(p2-1)*length(pow1) + p1] ^ pow3[p3]
                    }

                    T0s3[b, (p3-1)*length(pow1)*length(pow2) + (p2-1)*length(pow1) + p1] = aaa
                }
            }
        }
    }

                                        # permutation-based p-values:
    pPerm2 = rep(NA, length(pow1)*length(pow2)*length(pow3));
    pvs = NULL;

    for(j in 1:(length(pow1)*length(pow2)*length(pow3) )) {
        pPerm2[j] = sum( abs(Ts3[j]) < abs(T0s3[,j]))/n.perm
    }
    P0s3 = PermPvs(T0s3)
    minP0s2 = apply(P0s3, 1, min)
    minP2 =  sum( min(pPerm2) > minP0s2 )/n.perm

    pvs=c(pPerm2, minP2)

    nmvec <- NULL;
    for(kk in pow3) {
        for(ii in pow2) {
            for(jj in pow1) {
                nmvec <- c(nmvec, paste("MTSPUsPath",jj,",",ii,",",kk,sep=""))
            }
        }
    }

    nmvec <- c(nmvec, "MTaSPUsPath")
    names(pvs) <- nmvec
    pvs

}
