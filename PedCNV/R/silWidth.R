
silWidth <- function(dat,thres_sil=0.05,thres_MAF=0.01){

    old.order <- row.names(dat)
    obsNo <- dim(dat)[1]
    n <- dim(dat)[2]
    sdat <- dat[with(dat,order(dat[,n])),]
    clusRes <- sdat[,n]
    clusNo <- length(levels(factor(clusRes)))
    x <- as.matrix(sdat[,-n])
    t <- as.matrix(table(clusRes))
    breakp <- append(0,cumsum(t))
    thres_n <- ceiling(obsNo*thres_MAF)
    res <- .Call("sil_inter", x_ = x, clusRes_ = clusRes, clusNo_ = clusNo, obsNo_ = obsNo, breakp_= breakp)
    sil <- matrix(res$sil)
    rownames(sil) <- rownames(sdat)
    d0 <- res$d0
    rownames(d0) <- rownames(sdat)

    clus_adjust <- rep(NA,obsNo)
    for(i in 1:obsNo){
        if(abs(sil[i])<thres_sil) clus_adjust[i] <- -1
        else  clus_adjust[i] <- which(d0[i,]==min(d0[i,]))
    }
    clus_adjust_selec <- as.matrix(table(clus_adjust))

    abandon_clus <- -1
    for(i in 2:length(clus_adjust_selec))
        if(clus_adjust_selec[i,1]<thres_n) abandon_clus <- c(abandon_clus,rownames(clus_adjust_selec)[i])
    abandon_clus <- as.numeric(abandon_clus)

    sil_adjust <- matrix(NA,obsNo,1)
    rownames(sil_adjust) <- rownames(sil)

    for(i in 1:obsNo){
        if(clus_adjust[i]%in%abandon_clus)  sil_adjust[i,] <- NA
        else sil_adjust[i,] <- abs(sil[i,])
    }


    abandon.id <- row.names(sil_adjust)[which(is.na(sil_adjust))]
    silRes.old <- data.frame(clus=clusRes-1,sil=sil)
    silRes.new <- data.frame(clus=clus_adjust-1,sil=sil_adjust)
    silRes.adjust<- silRes.new[!is.na(silRes.new$sil),]

    ## change the order to the oringinal order
    silRes.old <- silRes.old[old.order,]
    old.order.abandon <- old.order[!old.order%in%abandon.id]
    silRes.adjust <- silRes.adjust[old.order.abandon,]

        
    clusAvg.adjust<- tapply(silRes.adjust$sil,silRes.adjust$clus,mean)
    silMean.adjust<- mean(silRes.adjust$sil)
    clusNum.adjust <- length(clusAvg.adjust)

    clusAvg<- tapply(silRes.old$sil,silRes.old$clus,mean)
    silMean<- mean(silRes.old$sil)
    clusNum <- length(clusAvg)
    
    return(list(adjusted=list(clusNum.adjust=clusNum.adjust,silMean.adjust=silMean.adjust,clusAvg.adjust=clusAvg.adjust,silRes.adjust=silRes.adjust,abandon.id=abandon.id),unadjusted=list(clusNum=clusNum,silMean=silMean,clusAvg=clusAvg,silRes=silRes.old)))
}
