project1 <-
function(v0, mat, matsd=NULL, estamb=FALSE,
         estdem=FALSE, equalsign=TRUE, stmat=NULL, fecundity1=TRUE){
  # makes one generation projection
  if(estamb==FALSE & estdem == FALSE) v1 <- mat%*%v0  
  if(estamb==FALSE & estdem == TRUE & is.null(stmat) & fecundity1 ==TRUE) v1 <- estdemo(v0, mat=mat)
  if(estamb==FALSE & estdem == TRUE & is.null(stmat) & fecundity1 ==FALSE) v1 <- estdemo(v0, mat=mat, fecundity1=FALSE)
  if(estamb==FALSE & estdem == TRUE & !is.null(stmat)) v1 <- estdemo(v0, mat=mat, stmat=stmat)

  if(estamb==TRUE & estdem ==FALSE){
     if (is.null(matsd)) stop ("there is not SD matrix provided
                                (argument matsd=NULL)")
     v1<- estambi(mat=mat, matsd=matsd, equalsign=equalsign)%*% v0
  }
  if (estamb==TRUE & estdem == TRUE){
    if (is.null(matsd)) stop ("there is not SD matrix provided
                               (argument matsd=NULL)")
    if(is.null(stmat) & fecundity1 ==TRUE) v1 <- estdemo(v0, mat= estambi(mat,matsd, equalsign=equalsign))
    if(is.null(stmat) & fecundity1 ==FALSE) v1 <- estdemo(v0, mat= estambi(mat,matsd, equalsign=equalsign), fecundity1=FALSE)
    if(!is.null(stmat)) v1 <- estdemo(v0, mat= estambi(mat,matsd, equalsign=equalsign), stmat=stmat)
  }
  return (v1)
}

