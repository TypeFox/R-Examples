combfunc <-
function(kateg.zahl, item.zahl, eps.mat, patmat.o){

last <- unique(patmat.o[,kateg.zahl])

#create gammat-matrix

gammat <- vector("list",length(last))
gammat[[length(last)]] <- 1 


#pattern

fsum.s <- list()

lpart <- matrix(1.00,nrow=1,ncol=item.zahl+1)

hilfmat <- diag(rep(1,kateg.zahl-1))

patsplit <- split(patmat.o, patmat.o[,kateg.zahl])

patsplitmat <- lapply(patsplit, function(la) matrix(la, ncol=ncol(patmat.o)))

patmat.last <- as.matrix(patsplitmat[[item.zahl+1]][,-kateg.zahl])


  for (x in last[-1]){
    patmat.p <- patsplitmat[[x+1]][,-kateg.zahl]

    actpart <- matrix(NA,nrow=nrow(patmat.p),ncol=item.zahl+1)
    patmat.pl <- split(patmat.p, row(patmat.p))    
  
    for (i in seq(along=patmat.pl)){
      avo.r <- patmat.pl[[i]]-hilfmat
      ind.r <- which(colSums(patmat.pl[[i]]-hilfmat >= 0) == length(patmat.pl[[i]]))
      avo.rb <- matrix(avo.r[,ind.r], ncol=sum(colSums(patmat.pl[[i]]-hilfmat >= 0) == length(patmat.pl[[i]])))
      
      ind <- which(colSums(mapply(function(a,b) avo.rb[,a] == patmat.last[,b], a=rep(1:ncol(avo.rb), each=ncol(patmat.last)), b=1:ncol(patmat.last)))==nrow(avo.rb))
      ind2 <- ifelse(ind <= ncol(patmat.last), ind, ind %% ncol(patmat.last))
    
        for (z in (sum(!is.na(lpart[1,]))-1):1){
          fsum <- lpart[ind2,z+1]         
          actpart[i,z] <- sum(actpart[i,z+1],  eps.mat[z,ind.r]*fsum, na.rm=T)
        }
     
  }  
  gammat[[x+1]] <- actpart[,1]

  act <- matrix(NA, nrow=nrow(patmat.p), ncol=ncol(patmat.p))
  for (j in seq_len(ncol(patmat.p))){
    act[(patmat.p[,j]-1>=0),j] <- lpart[,2]/actpart[(patmat.p[,j]-1 >=0),1]  
  }

  fsum.s <- c(fsum.s, list(act))
      
  lpart <- actpart
  patmat.last <- t(patmat.p)
}

gun <- unlist(rev(gammat))

quot <- rbind(rep(0,(kateg.zahl-1)),do.call(rbind, fsum.s))

return(list(gammat = gun, gam.quot=quot,patmat.o = patmat.o))
}
