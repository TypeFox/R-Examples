seqauto <- function(seqdata, order=1, measure="cv"){

    available.measures <- c("cv")

    srs <- seqformat(seqdata, from="STS", to="SRS")
    nc <- ncol(srs)

    x <- list()
    aa <- matrix(NA,order,2)

    for (i in 1:order) {
        x[[i]] <- table(srs[,nc],srs[,nc-i])
        if (measure=="cv"){
            chi.test <- suppressWarnings(chisq.test(x[[i]], correct=FALSE))
            chisq <- chi.test$statistic
            names(chisq) <- ""
            aa[i,1] <- sqrt(chisq/(sum(x[[i]])*min(dim(x[[i]])-1)))
            aa[i,2] <- chi.test$p.value
        }
        else {
            stop(paste("Measure",measure,"not available, should be one of",available.measures))
        }
    }
    rownames(aa) <- paste(measure,"(",1:order,")",sep="")
    colnames(aa) <- c("auto-association","p-value")
    return(aa)

}
