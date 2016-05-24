#' GATES-Simes test for single trait - pathway association.
#'
#' Get the p-value of GATES-Simes. It uses an extended Simes procedure to combine GATES p-values across multiple genes in a pathway.
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
#' Hongsheng Gui, Miaoxin Li, Pak C Sham and Stacey S Cherny (2011)
#' Comparisons of seven algorithms for pathway analysis using the WTCCC Crohn's Disease
#' BMC Research Notes, 4:386
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
#' out <- GatesSimes(pvec = logitp, ldmatrix = ldmat, snp.info = simula$snp.info,
#'                   gene.info = simula$gene.info)
#' out
#'
#'
#' @seealso \code{\link{Hyst}}  \code{\link{GATES2}}

GatesSimes <- function(pvec, ldmatrix, snp.info, gene.info) {

    n.gene <- nrow(gene.info)

    GL <- list(0)
    i = 1
    for(g in 1:n.gene) { # g = 2

        snpTF <- ( snp.info[,2] == gene.info[g,2] &
                      gene.info[g,3] <= as.numeric(snp.info[,3]) &
                          gene.info[g,4] >= as.numeric(snp.info[,3]) )

        if( sum(snpTF) != 0) {
            GL[[i]] <- which(snpTF)
            i = i + 1
        }
    }

    n.gene <- length(GL)

    PGs <- NULL;
    for ( i in 1:length(GL) ) { #i = 26
        if ( length(GL[[i]]) == 1 ) {
            Pg <- pvec[ GL[[i]] ]
        } else {
            pvecG <- pvec[ GL[[i]] ]
            ldmat <- ldmatrix[GL[[i]], GL[[i]]]
            o.pv <- order(pvecG)
            ldmat2 <- ldmat[o.pv,o.pv]
            Pg <- GATES2(ldmatrix = ldmat2, p = sort(pvecG) )[1]
        }
        PGs <- c(PGs, Pg)
    }

    PgatesSime <- min( PGs * n.gene / rank(PGs) )
    PgatesSime
}

