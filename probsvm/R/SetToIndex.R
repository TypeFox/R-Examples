"SetToIndex" <- 
function(n, Elbow, Left){
   set.id <- rep(1,n)
   set.id[Elbow] <- 0
   set.id[Left] <- -1
set.id
}
