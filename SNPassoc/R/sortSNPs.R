`sortSNPs` <-
function (data, colSNPs, info) 
{
    o <- order(info[, 2], info[, 3])
    label.SNPs.o <- info[o, 1]
    label.SNPs <- names(data[, colSNPs, drop=FALSE])

#control
    ans <- match(label.SNPs, label.SNPs.o)
    if (sum(is.na(ans)) > 0) {
        cat("Warning: ")
        cat("the SNPs: ", as.character(label.SNPs[is.na(ans)]), 
            "\n  are not included in the file with the genomic positions and they are discarded \n")
    }

    ans <- match(label.SNPs.o, label.SNPs)


    out <- colSNPs[ans[!is.na(ans)]]
    out <- out[!is.na(out)]
    res <- list(pos=out, dataSorted=info[o,])
    res 
}

