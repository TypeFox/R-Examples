correl <- function(measures, limit=0.85, plot.scatter=FALSE, keep=NA){
    
    to.keep <- which(names(measures) %in% keep)
    answer <- list()
    correliert <- which(abs(cor(measures, use="complete.obs"))>limit)-1
    cor.paare <- cbind(correliert%/%NCOL(measures),correliert%%NCOL(measures))+1
    cor.paare <- cor.paare[cor.paare[,1]!=cor.paare[,2],]
    cor.paare <- cor.paare[cor.paare[,1]<cor.paare[,2],]

    answer$pairs <- cor.paare
    answer$pairs.by.name <- cbind(names(measures)[cor.paare[,1]],names(measures)[cor.paare[,2]])


    if(plot.scatter){
        opar <- par(mfrow=c(5,4),mar=c(4,4,0,0)+0.1)
        for(i in 1:NROW(cor.paare)){
                a <- cor.paare[i,1]
                b <- cor.paare[i,2]

                plot(measures[,a],measures[,b], xlab=names(measures)[a], ylab=names(measures)[b])
                if(i%%20==0) par(ask=TRUE)
                if(i%%20==1) par(ask=FALSE)

#                a <- measures[,cor.paare[i,1]]
#                b <- measures[,cor.paare[i,2]]
#
#                a[is.infinite(a)] <- NA
#                b[is.infinite(b)] <- NA
#                browser()
#                print(plot(hexbin(a ~b),
#                xlab=names(measures)[cor.paare[i,1]],
#                ylab=names(measures)[cor.paare[i,2]]),
#                split=c(5,4,5,4), more=T)
        }
        par(opar)
    }


    i=1
    toDrop <- cor.paare
    toDrop[toDrop[,2] %in% to.keep,] <- toDrop[toDrop[,2] %in% to.keep,c(2,1)]
    
    while(i < NROW(toDrop)){
            a <- toDrop[i,2]
            toDrop <- toDrop[toDrop[,1]!=a,] #den gerade rausgeschm. aus
            i = i+1
    }
    current=-1
    excluded <- list()

    for(i in 1:NROW(toDrop)){
       if(toDrop[i,1]!=current){
           if(i!=1){
               excluded[[names(measures)[current]]] <- names(measures)[entfernt]
           }
           entfernt <- c()
           current <- toDrop[i,1]
       } 
       entfernt <- c(entfernt, toDrop[i,2])
    }
   excluded[[names(measures)[current]]] <- names(measures)[entfernt]

   answer$possible.exclusion <- excluded

    answer$to.drop <- toDrop <- unique(toDrop[,2])
    answer$to.keep <- (1:length(names(measures)))[-toDrop]
    answer$to.drop.by.name <-names(measures)[toDrop]
    answer$to.keep.by.name <-names(measures)[answer$to.keep]
    return(answer)
}
