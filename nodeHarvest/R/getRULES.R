getRULES <-
function(tree,RULES,maxinter=2,levelvec=levelvec){


  norules <- sum( ind <- ((tree[,"level"]<= maxinter + 1) & (tree[,"level"] > 1)))
  RULES <- vector("list",norules)
  loop <- which( ind )
  
  for (rc in 1:norules){

    ruletmp <- matrix( nrow= tree[ loop[rc], "level"]-1, ncol=3)
    colnames(ruletmp) <- c("variable","lower","upper")
    tree[,"split point"] <- signif(tree[,"split point"],3)

    node <- loop[rc]
    depth <- tree[ loop[rc],"level"]
    while( depth >1){
      nodenew <- which( apply( tree[1:(node-1),c("left daughter","right daughter"),drop=FALSE]==node,1,sum) >=1 )
      left <- tree[nodenew,"left daughter"]==node
      depth <- depth-1
      nfact <- length(levelvec[[ tree[nodenew,"split var"]]]) 
      if(nfact==0){
        ruletmp[ depth,] <- c(tree[nodenew,"split var"], if(left) -Inf else tree[nodenew,"split point"], if(left) tree[nodenew,"split point"] else Inf)
      }else{
        tmp <- tree[nodenew,"split point"]
        ruletmp[depth,] <-  c(tree[nodenew,"split var"], if(left) tmp else round(dectobin(!(dectobin(tmp,nl=nfact)==1),forward=FALSE)) , Inf)
      }
      node <- nodenew
    }

    while(any( (tab <- table( ruletmp[,"variable"]))>1)){
      sel <- as.numeric( names( which(tab>1)[1]))[1]
      ind <- which( ruletmp[,"variable"]==sel)
      nfact <- length(levelvec[[sel]])
      if(nfact==0){
        ruletmp[ ind[1], ] <- c( sel, max(ruletmp[ind,"lower"]),min(ruletmp[ind,"upper"]))
      }else{
        tmp <- rep(1, nfact)
        for (kcc in ind) tmp <- tmp * dectobin(ruletmp[kcc,"lower"],nl=nfact)
        ruletmp[ind[1],] <- c(sel, round(dectobin(tmp==1,forward=FALSE)) ,Inf)
      }
      ruletmp <- ruletmp[ -ind[2:length(ind)], ,drop=FALSE]
    }
    

    RULES[[rc]] <- ruletmp
    
    vari <- ruletmp[,c("lower","upper")]
    attr(RULES[[rc]],"depth") <- sum( abs(vari) < Inf   )
    
  }
  return(RULES)
    
  
}

