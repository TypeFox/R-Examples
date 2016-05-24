#' HYST (Hybrid Set-based Test) for single trait - pathway association 
#'
#' Get a p-value using HYST.
#'
#' @param pvec p-values for each SNP.
#'
#' @param ldmatrix numeric. A correlation matrix of SNPs with dimensions matching the length of pvec (the number of SNPs).
#'
#' @param snp.info SNP information matrix, the 1st column is SNP id, 2nd column is chromosome #, 3rd column indicates SNP location.
#'
#' @param gene.info GENE information matrix, The 1st column is GENE id, 2nd column is chromosome #, 3rd and 4th column indicate start and end positions of the gene.
#'
#' @export
#' @return A p-value.
#'
#' @author Il-Youp Kwak and Wei Pan
#'
#' @references
#' Miao-Xin Li, Johnny S.H. Kwan and Pak C. Sham (2012)
#' HYST: A Hybrid Set-Based Test for Genome-wide Association Studies, with Application to Protein-Protein Interaction-Based Association Analysis
#' The American Journal of Human Genetics, 91, 478-488.
#'
#'
#' @examples
#'
#' simula <- simPathAR1Snp(nGenes=20, nGenes1=1, nSNPlim=c(1, 20), nSNP0=1:3,
#'                            LOR=.2, rholim=c(0,0),
#'                            n=100, MAFlim=c(0.05, 0.4), p0=0.05)
#' logitp <- getlogitp(simula$Y, simula$X)
#'
#' ## get correlation of SNPs using controls
#' ldmat <- cor(simula$X[ simula$Y == 0, ])
#' out <- Hyst(pvec = logitp, ldmatrix = ldmat, snp.info = simula$snp.info,
#'             gene.info = simula$gene.info)
#'
#' @seealso \code{\link{GatesSimes}} \code{\link{GATES2}}

Hyst <- function(pvec, ldmatrix, snp.info, gene.info) {
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

    n.gene <- length(GL)


    PGs <- NULL;
    keyGs <- NULL;
    gpos <- 0

    for ( i in 1:length(GL) ) { #i = 26
        if ( length(GL[[i]]) == 1 ) {
            Pg <- pvec[ GL[[i]] ]
            keyG <- GL[[i]]
        } else {
            pvecG <- pvec[ GL[[i]] ]
            ldmat <- ldmatrix[GL[[i]], GL[[i]]]
            o.pv <- order(pvecG)
            ldmat2 <- ldmat[o.pv,o.pv]
            out <- GATES2(ldmatrix = ldmat2, p = sort(pvecG) )
            Pg <- out[1]
#            keyG <- GL[[i]][o.pv][out[2]]
            keyG <- GL[[i]][out[2]]
        }
        PGs <- c(PGs, Pg)
        keyGs <- c(keyGs, keyG)
    }
#    cat("keyGs :", keyGs);
#    cat("\n");
#    cat("PGs :", PGs);

    Hyst <- -2 * sum( log(PGs) )

    sums <- 0; # i = 1 ; j = 2
    for(i in 1:(n.gene-1) )
        for(j in (i+1):n.gene) {
            r = abs(ldmatrix[keyGs[i], keyGs[j]])
#            r = abs(cor(X[,keyGs[i]], X[,keyGs[j]]) )
            sums = sums + r * (3.25 + 0.75 *r)
        }
    c = 1 + sums / (2 * n.gene)
    f = 2* n.gene / c

    Physt = 1 - pchisq( Hyst / c, df = f)
#    Physt2 = 1 - pchisq( Hyst , df = 2*n.gene)
}
