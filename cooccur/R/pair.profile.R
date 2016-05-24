pair.profile <-
function(mod){
  ptab <- mod$results
  spp <- unique(c(ptab$sp1,ptab$sp2))
  profiles <- data.frame(spp=c(), rand=c(),pos=c(), neg=c())
  for (i in as.character(spp)){
    spptab <- ptab[ptab$sp1 == i| ptab$sp2 == i,]
    
    
    pos <- nrow(spptab[spptab$p_gt < 0.05,])
    neg <- nrow(spptab[spptab$p_lt < 0.05,])
    rand <- nrow(spptab) - (pos + neg)
    total <- rand + pos + neg
    if (total == 0){
        profiles[i,"pos"] <- 0
        profiles[i,"neg"] <- 0
        profiles[i,"rand"] <- 100
        profiles[i,"spp"] <- i
    }else{
        profiles[i,"pos"] <- pos / total * 100
        profiles[i,"neg"] <- neg / total * 100
        profiles[i,"rand"] <- rand / total * 100
        profiles[i,"spp"] <- i
    }
  }
    profiles$cooccur <- profiles$pos + profiles$neg 
    profiles <- profiles[with(profiles,order(cooccur)),]

  ## Reorder fullname based on the the sum of the other columns
  profiles$spp <- factor(x=profiles$spp,ordered=T,levels=profiles$spp)
  profiles$cooccur <- NULL

  pos <- mod$positive
  neg <- mod$negative
  rand <- mod$random
  total <- pos + neg + rand
  pos <- mod$positive / total * 100
  neg <- mod$negative / total * 100
  rand <- mod$random / total * 100
  
  if ("sp1_name" %in% colnames(ptab)){
    spp_key <- unique(data.frame(sppnum=c(ptab$sp1,ptab$sp2),sppname=c(as.character(ptab$sp1_name),as.character(ptab$sp2_name))))
    
    names <- merge(x=data.frame(order=1:length(row.names(profiles)),profiles),y=spp_key,by.x="spp",by.y="sppnum",all.x=T)
    
    names <- names[with(names,order(order)),]  
    names$order <- NULL
    names$spp <- NULL
      names$sppname <- factor(x=names$sppname,ordered=T,levels=names$sppname)
    sidebar <- data.frame(pos=pos,neg=neg,rand=rand,sppname="All Species")

      names <- rbind(names,sidebar)
      md <- melt(names, id=(c("sppname")))
    
    sppname <- md$sppname
    value <- md$value
    variable <- md$variable  
    
  p <- ggplot(data=md, aes(x=sppname, y=value, fill=variable) ) + geom_bar(stat = "identity",position = "stack",width=1)
  p <- p + scale_fill_manual(values = c("light blue","#FFCC66","dark gray"), name = "", labels = c("positive","negative","random")) + theme(plot.title = element_text(vjust=2,size=20, face="bold"),legend.text=element_text(size=18),axis.title.y = element_text(size = 20),axis.text.y=element_text(size=18),axis.text.x=element_text(size=12,hjust=0,vjust=1,angle=-45), panel.background = element_rect(fill='white', colour='black'),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("") + ylab("Percent of pairings")
  p <- p + ggtitle("Species Association Profile") + scale_y_continuous(expand = c(0.005,0)) + scale_x_discrete(expand=c(0.005, 0))
      p <- p + geom_bar(stat = "identity",data=md[(md$sppname=="All Species"),], aes(sppname), alpha=0, size=0.1, color="white",width=1)
  p
  }else{
      sidebar <- data.frame(pos=pos,neg=neg,rand=rand,spp="All Species")

      profiles <- rbind(profiles,sidebar)
  md <- melt(profiles, id=(c("spp")))
  p <- ggplot(data=md, aes(x=spp, y=value, fill=variable) ) + geom_bar(stat = "identity",position = "stack",width=1)
  p <- p + scale_fill_manual(values = c("light blue","#FFCC66","dark gray"), name = "", labels = c("positive","negative","random")) + theme(plot.title = element_text(vjust=2,size=20, face="bold"),legend.text=element_text(size=18),axis.title.y = element_text(size = 20),axis.text.y=element_text(size=18),axis.text.x=element_text(size=12,hjust=0,vjust=1), panel.background = element_rect(fill='white', colour='black'),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("") + ylab("Percent of pairings")
  p <- p + ggtitle("Species Association Profile") + scale_y_continuous(expand = c(0.005,0)) + scale_x_discrete(expand=c(0.005, 0))
    p <- p + geom_bar(stat = "identity",data=md[(md$spp=="All Species"),], aes(spp), alpha=0, size=0.1, color="white",width=1) + theme(axis.text.x=element_text(hjust=0,vjust=1,angle=-45))
  p
  }

}
