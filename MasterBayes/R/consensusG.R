consensusG<-function (GdP, cat.levels = NULL, gmax=FALSE, het=FALSE) 
{
    Gid <- GdP$id[-duplicated(GdP$id) == FALSE]
    G <- lapply(GdP$G, function(x) {
        x[-which(duplicated(GdP$id) == TRUE)]
    })
    G<-lapply(G, function(x){
       x[which(is.na(x)==FALSE)]<-NA
       x}
       )
    if (is.null(cat.levels)) {
        cat.levels <- 1
        GdP$categories <- rep(1, length(GdP$id))
    }else{
        if (any(GdP$categories %in% cat.levels == FALSE)) {
            stop("some GdP$categories not appearing in cat.levels")
        }
    }
    for (i in 1:length(G)) {
       if(gmax){
         nG<-tapply(GdP$G[[i]],GdP$id, function(x){   
           ngs<-table(as.character(x))
           if(length(ngs)>0){
             if(length(which(ngs==max(ngs)))==1){
               names(ngs)[which.max(ngs)]
             }else{
               NA
             }
           }else{
             NA
           }
         })
         G[[i]][match(names(nG), Gid)]<-as.genotype(c(nG))
       }

       miss<-is.na(GdP$G[[i]])*(1e+10)
       ehet<-homozygote(GdP$G[[i]])*(het*100)
       ehet[which(is.na(ehet))]<-0
       ecat<-match(GdP$categories, cat.levels)

       nG <- GdP$G[[i]][order(miss+ehet+ecat)]
       nid <- GdP$id[order(miss+ehet+ecat)]

       nG <- nG[-duplicated(nid) == FALSE]
       nid <- nid[-duplicated(nid) == FALSE]
       nG <- nG[match(Gid, nid)]

      G[[i]][which(is.na(G[[i]]))] <- nG[which(is.na(G[[i]]))]
    }
    GdataPed(id = Gid, G = G, perlocus = GdP$perlocus, marker.type = GdP$marker.type)
}

