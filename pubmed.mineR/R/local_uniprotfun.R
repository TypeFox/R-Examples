local_uniprotfun = function (y) 
{
    check11 = which(as.character(HGNCdata[, 2]) == y)
    check11a = as.character(HGNCdata[check11, 1])
    check12 = which(as.character(HGNC2UniprotID[, 1]) == check11a[1])
    check13 = as.character(HGNC2UniprotID[check12, 2])
    for (i in 1:length(check13)) {
        temp = getURL(paste("http://www.uniprot.org/uniprot/", 
            check13[i], ".txt", sep = ""))
        res1 = temp
        write(res1, file = paste("x", ".txt", sep = ""))
    }
}

