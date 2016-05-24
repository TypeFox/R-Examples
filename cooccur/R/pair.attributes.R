pair.attributes <-
function(mod){
  ptab <- mod$results
  spp <- unique(c(ptab$sp1,ptab$sp2))
  profiles <- data.frame(spp=c(), rand=c(),pos=c(), neg=c(), num_rand=c(),num_pos=c(), num_neg=c())
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
        profiles[i,"num_pos"] <- pos
        profiles[i,"num_neg"] <- neg
        profiles[i,"num_rand"] <- rand
        profiles[i,"spp"] <- i
    }else{
        profiles[i,"pos"] <- round(pos / total * 100,2)
        profiles[i,"neg"] <- round(neg / total * 100,2)
        profiles[i,"rand"] <- round(rand / total * 100,2)
        profiles[i,"num_pos"] <- pos
        profiles[i,"num_neg"] <- neg
        profiles[i,"num_rand"] <- rand
        profiles[i,"spp"] <- i
    }
  }
    profiles$cooccur <- profiles$pos + profiles$neg 
    profiles <- profiles[with(profiles,order(cooccur)),]

  ## Reorder fullname based on the the sum of the other columns
  profiles$spp <- factor(x=profiles$spp,ordered=T,levels=profiles$spp)
  profiles$cooccur <- NULL
  
  if ("sp1_name" %in% colnames(ptab)){
    spp_key <- unique(data.frame(sppnum=c(ptab$sp1,ptab$sp2),sppname=c(as.character(ptab$sp1_name),as.character(ptab$sp2_name))))
    
    names <- merge(x=data.frame(order=1:length(row.names(profiles)),profiles),y=spp_key,by.x="spp",by.y="sppnum",all.x=T)
    
    names <- names[with(names,order(order)),]  
    names$order <- NULL
    names$spp <- NULL

names[order(match(names$sppname,mod$spp.names)),]
  }else{
profiles[order(match(profiles$spp,mod$spp.names)),]
  }

}
