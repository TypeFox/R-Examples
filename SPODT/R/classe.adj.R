classe.adj <-
function(bord)
{
    bord <-bord[order(bord$id),]
    
    bordX <- bord$X
    bordY <- bord$Y

    tab <- table(bord$id)
    nb.cl <- length(tab)
    d.cl <- c(0,cumsum(as.vector(tab)))
  
    adj <- c(rep(c(1,rep(0,nb.cl)),nb.cl-1),1)
    
    res <- .C("classeAdj",
              as.integer(nb.cl), as.integer(d.cl), as.double(bordX),  as.double(bordY),
              as.integer(adj)
              
             )#ajouté ,PACKAGE="SPODT"
             
    id <- as.character(names(tab))
                 
    adj <- matrix(res[[5]], nrow=nb.cl, ncol=nb.cl, dimnames=list(id, id))

    return(list(adj=adj, bord=bord))
}
