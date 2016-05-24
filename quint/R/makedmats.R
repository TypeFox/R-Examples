makedmats <-
function(dmat){
  #check of boundary condition: partition class cardinality condition:P1 and P2 may not be empty
  #creates D'(dmats): matrix D with admissible assignments(K'); dmats= K' * I matrix; 
  sel1<-numeric(dim(dmat)[1])
  sel2<-numeric(dim(dmat)[1])
  #count the assignments to p1 for each row of dmat
  sel1<-apply(dmat==1,1,sum)
  #count the assignments to p2 for each row of dmat
  sel2<-apply(dmat==2,1,sum)
  #select the rows for which sel1 & sel2 not equals 0
  dmats<-dmat[sel1&sel2!=0,] 
  return(dmats)}
