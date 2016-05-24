amari_distance <-
function(Q1,Q2){

 m=nrow(Q1)
 Per=Q1%*%ginv(Q2)
 Perf=(sum(apply(abs(Per),1,sum)/apply(abs(Per),1,max)-1)/(m-1)+sum(apply(abs(Per),2,sum)/apply(abs(Per),2,max)-1)/(m-1))/(2*m) 
 return(Perf)

}
