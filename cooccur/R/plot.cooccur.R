plot.cooccur <-
function(x, ...){
  
  ##
  allargs <- match.call(expand.dots = TRUE)
  plotrand <- allargs$plotrand
    plotrand <- ifelse(test = is.null(plotrand),yes = FALSE,no = plotrand)
  randsummary<- allargs$randsummary
    randsummary <- ifelse(test = is.null(randsummary),yes = FALSE,no = randsummary)

  ##
  
  dim <- x$species
  comat_pos <- comat_neg <- matrix(nrow=dim,ncol=dim)
  
  co_tab <- x$result
  for (i in 1:nrow(co_tab)){
    comat_pos[co_tab[i,"sp1"],co_tab[i,"sp2"]] <- co_tab[i,"p_gt"]
    comat_pos[co_tab[i,"sp2"],co_tab[i,"sp1"]] <- co_tab[i,"p_gt"]
    
    row.names(comat_pos[co_tab[i,"sp2"],co_tab[i,"sp1"]])
    
  }
  for (i in 1:nrow(co_tab)){
    comat_neg[co_tab[i,"sp1"],co_tab[i,"sp2"]] <- co_tab[i,"p_lt"]
    comat_neg[co_tab[i,"sp2"],co_tab[i,"sp1"]] <- co_tab[i,"p_lt"]
  }
  comat <- ifelse(comat_pos>=0.05,0,1) + ifelse(comat_neg>=0.05,0,-1)
  colnames(comat) <- 1:dim
  row.names(comat) <- 1:dim
  
  if ("spp_key" %in% names(x)){
    
    sp1_name <- merge(x=data.frame(order=1:length(colnames(comat)),sp1=colnames(comat)),y=x$spp_key,by.x="sp1",by.y="num",all.x=T)
    sp2_name <- merge(x=data.frame(order=1:length(row.names(comat)),sp2=row.names(comat)),y=x$spp_key,by.x="sp2",by.y="num",all.x=T)
    
    colnames(comat) <- sp1_name[with(sp1_name,order(order)),"spp"]  
    row.names(comat) <- sp2_name[with(sp2_name,order(order)),"spp"]
      
  }  
  
  #ind <- apply(comat, 1, function(x) all(is.na(x)))
  #comat <- comat[!ind,]
  #ind <- apply(comat, 2, function(x) all(is.na(x)))
  #comat <- comat[,!ind]

  comat[is.na(comat)] <- 0
  
  origN <- nrow(comat)
  
  # SECTION TO REMOVE SPECIES INTERACTION WITH NO OTHERS

  #rmrandomspp <- function(orimat,plotrand = FALSE,randsummary = FALSE){
  if(plotrand == FALSE){
    ind <- apply(comat, 1, function(x) all(x==0))
    comat <- comat[!ind,]    
    ind <- apply(comat, 2, function(x) all(x==0))
    comat <- comat[,!ind]
    #ind <- apply(orimat, 1, function(x) all(x==0))
    #orimat <- orimat[!ind,]    
    #ind <- apply(orimat, 2, function(x) all(x==0))
    #orimat <- orimat[,!ind]
  }
  #return(orimat)
  #}
  
  #comat <- rmrandomspp(orimat = comat, dots)
  ####################################################### 

  postN <- nrow(comat)


  comat <- comat[order(rowSums(comat)),]
  comat <- comat[,order(colSums(comat))]
  
  #comat <- rmrandomspp(orimat = comat, ...)
  
  #ind <- apply(comat, 1, function(x) all(x==0))
  #comat <- comat[!ind,]
  #ind <- apply(comat, 2, function(x) all(x==0))
  #comat <- comat[,!ind]
  
  ind <- apply(comat, 1, function(x) all(x==0))
  comat <- comat[names(sort(ind)),]
  ind <- apply(comat, 2, function(x) all(x==0))
  comat <- comat[,names(sort(ind))]
  
  #comat
  data.m = melt(comat)
  colnames(data.m) <- c("X1","X2","value")
  data.m$X1 <- as.character(data.m$X1)
  data.m$X2 <- as.character(data.m$X2)
 
  meas <- as.character(unique(data.m$X2))

  dfids <- subset(data.m, X1 == X2)
  
  X1 <- data.m$X1
  X2 <- data.m$X2
  
  df.lower = subset(data.m[lower.tri(comat),],X1 != X2)
  
  ##### testing the rand summary
    if(randsummary == FALSE){  
      }else{
        dim <- nrow(comat)
        ext.dim <- round(dim*0.2,digits = 0)
          if(ext.dim<0){ext.dim<-1}
        placehold <- paste("ext_", rep(c(1:ext.dim),each = dim), sep="")
      
        randcol.df <- data.frame(
          X1 = placehold,
          X2 = rep(meas,times = ext.dim),
          value = rep(x = c(-2), times = dim*ext.dim))
        
        df.lower <- rbind(df.lower,randcol.df)
        meas <- c(meas,unique(placehold))
      }

  


  #####

  X1 <- df.lower$X1
  X2 <- df.lower$X2
  value <- df.lower$value
 


  ####
 if(randsummary == FALSE){  
    p <- ggplot(df.lower, aes(X1, X2)) + geom_tile(aes(fill = factor(value,levels=c(-1,0,1))), colour ="white") 
 p <- p + scale_fill_manual(values = c("#FFCC66","dark gray","light blue"), name = "", labels = c("negative","random","positive"),drop=FALSE) + 
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),plot.title = element_text(vjust=-4,size=20, face="bold"),panel.background = element_rect(fill='white', colour='white'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = c(0.9, 0.5),legend.text=element_text(size=18)) + 
    ggtitle("Species Co-occurrence Matrix") + 
    xlab("") + ylab("") + 
    scale_x_discrete(limits=meas, expand = c(0.3, 0),drop=FALSE) + 
    scale_y_discrete(limits=meas, expand = c(0.3, 0),drop=FALSE) 
 p <- p + geom_text(data=dfids,aes(label=X1),hjust=1,vjust=0,angle = -22.5)#, color="dark gray")
 
   
      }else{
        
         p <- ggplot(df.lower, aes(X1, X2)) + geom_tile(aes(fill = factor(value,levels=c(-1,0,1,-2))), colour ="white") 
 p <- p + scale_fill_manual(values = c("#FFCC66","dark gray","light blue","light gray"), name = "", labels = c("negative","random","positive","random"),drop=FALSE) + 
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),plot.title = element_text(vjust=-4,size=20, face="bold"),panel.background = element_rect(fill='white', colour='white'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = c(0.9, 0.5),legend.text=element_text(size=18)) + 
    ggtitle("Species Co-occurrence Matrix") + 
    xlab("") + ylab("") + 
    scale_x_discrete(limits=meas, expand = c(0.3, 0),drop=FALSE) + 
    scale_y_discrete(limits=meas, expand = c(0.3, 0),drop=FALSE) 
 p <- p + geom_text(data=dfids,aes(label=X1),hjust=1,vjust=0,angle = -22.5)#, color="dark gray")
 
 dim <- nrow(comat)
 ext_x <- dim + 0.5 #(ext.dim/2)
 ext_y <- dim + 1
 nrem <- origN - postN
 randtext <- paste(nrem, " completely\nrandom species")
 ext_dat <- data.frame(ext_x=ext_x,ext_y=ext_y,randtext=randtext)

  p <- p + geom_text(data=ext_dat,aes(x = ext_x,y = ext_y,label=randtext),hjust=0,vjust=0, color="dark gray")
}
 ####
 
 p
 
}
