## ----knitr_options, echo=FALSE, results=FALSE----------------------------
library(knitr)
opts_chunk$set(fig.width = 12)

## ----loading-------------------------------------------------------------
library(aSPU)
data(kegg9)

## ----kegg1---------------------------------------------------------------
kegg9$gene.info

## ----kegg2---------------------------------------------------------------
## SNPs mapped on 3rd and 4th genes of 9th Kegg pathway
kegg9$PPs[3:4]

## ----kegg3---------------------------------------------------------------
length(kegg9$Ps)
kegg9$Ps[1:10]

## ----kegg4---------------------------------------------------------------
kegg9$nP[1:10]
kegg9$ldmatrix[1:10,1:10]

## ----kegg5---------------------------------------------------------------
kegg9$snp.info[1:10,]

## ----gene-wise p---------------------------------------------------------
Gps<-NULL;
gl <- NULL;
for( g in kegg9$gene.info[,1]) {
    snps <-
    which( ( kegg9$snp.info[,2] == kegg9$gene.info[kegg9$gene.info[,1]==g,2]) &
           (kegg9$snp.info[,3] > kegg9$gene.info[kegg9$gene.info[,1] == g, 3])&
           (kegg9$snp.info[,3] < kegg9$gene.info[kegg9$gene.info[,1] == g, 4]))

    newP <- kegg9$nP[snps] ;
    ldsub <- kegg9$ldmatrix[snps, snps];

    if( length(snps) > 1) {
        out <- aSPUs(newP, corrSNP=ldsub , pow=c(1,2,4,8, Inf),
                     n.perm=10000, Ps=TRUE)

        o.pvec = order(newP)
        ldmat <- ldsub[o.pvec, o.pvec]
        gatesp <- GATES2(ldmat, sort(newP))[1]
        Gps <- rbind(Gps, c(length(newP),out$pvs, gatesp))
        gl <- c(gl, g)
    } else if (length(snps) == 1) {
        out <- aSPUs(newP, corrSNP=ldsub , pow=c(1,2,4,8, Inf),
                     n.perm=10000, Ps=TRUE)
        gatesp <- newP
        Gps <- rbind(Gps, c(length(newP),out$pvs, gatesp) )
        gl <- c(gl, g)
    }
}
colnames(Gps)[1] <- "nSNP"
rownames(Gps) <- gl
Gps

## ----pathway2------------------------------------------------------------
g = "LIPA"
    snps <-
    which( ( kegg9$snp.info[,2] == kegg9$gene.info[kegg9$gene.info[,1]==g,2]) &
           (kegg9$snp.info[,3] > kegg9$gene.info[kegg9$gene.info[,1] == g, 3])&
           (kegg9$snp.info[,3] < kegg9$gene.info[kegg9$gene.info[,1] == g, 4]))

newP <- kegg9$nP[snps] ;
newP


## ----pathway p-----------------------------------------------------------
out.g <- GatesSimes(pvec = kegg9$nP, ldmatrix = kegg9$ldmatrix,
                    snp.info=kegg9$snp.info, gene.info = kegg9$gene.info)

out.h <- Hyst(pvec = kegg9$nP, ldmatrix = kegg9$ldmatrix,
              snp.info=kegg9$snp.info, gene.info = kegg9$gene.info)

out.a <- aSPUsPath(kegg9$nP, corrSNP = kegg9$ldmatrix, pow=c(1,2,4,8, Inf),
                   pow2 = c(1,2,4,8),
                   snp.info=kegg9$snp.info, gene.info = kegg9$gene.info,
                   n.perm=1000, Ps = TRUE)

out.g; out.h; out.a


