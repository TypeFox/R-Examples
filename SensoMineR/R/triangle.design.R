triangle.design <- function(nbprod , nbpanelist, bypanelist = nbprod*(nbprod-1)/2, labprod = 1:nbprod, labpanelist = 1:nbpanelist){

  aux = as.data.frame(matrix(0,nbprod*(nbprod-1)/2,2))
  numligne = 0
  for (i in 1:(nbprod-1)){
    for (j in (i+1):nbprod){
      numligne = numligne + 1
      aux[numligne,1] = labprod[i]
      aux[numligne,2] = labprod[j]
    }
  }
  aux = aux[sample(nrow(aux)),]
  plan = as.data.frame(matrix(0,nbpanelist*bypanelist,4))
 
  num.test = 0
  ligne.aux = 0
  for (i in 1:nbpanelist){
    for (j in 1: bypanelist){
      ligne.aux = ligne.aux + 1
      if (ligne.aux == nrow(aux)+1) ligne.aux = 1
      num.test = num.test + 1
      plan[num.test,1] = labpanelist[i]
      plan[num.test,2:3] = aux[ligne.aux,]
    }
  }
  plan <- plan[order(plan[,2],plan[,3]),]
  for (i in 1:(nrow(plan))){
    if (i%%6 == 1) plan[i,4] = plan[i,2]
    if (i%%6 == 2) plan[i,4] = plan[i,3]
    if (i%%6 == 3){
      plan[i,4] = plan[i,3]
      plan[i,3] = plan[i,2]
    }
    if (i%%6 == 4){
      plan[i,4] = plan[i,2]
      plan[i,2] = plan[i,3]
      plan[i,3] = plan[i,4]
    }
    if (i%%6 == 5){
      plan[i,4] = plan[i,3]
      plan[i,3] = plan[i,2]
      plan[i,2] = plan[i,4]
    }
    if (i%%6 == 0){
      plan[i,4] = plan[i,2]
      plan[i,2] = plan[i,3]
    }
  }
  plan = plan[sample(nrow(plan)),]
#  plan = plan[order(plan[,1]),]
  
#  row.names(plan) = paste("Panelist",plan[,1],".Test",1:bypanelist,sep="")
  row.names(plan) = paste("Panelist",rep(1:nbpanelist,each=bypanelist),".Test",1:bypanelist,sep="")
  plan = plan[,-1]
  colnames(plan) = c("Product X","Product Y","Product Z")
  plan = as.data.frame(plan)
  return(plan)
}
