get.counts <-
function( x, state.num){

count.mat<-matrix(0, nrow=state.num, ncol=state.num)

   row.len<-length(x)

   for(i in 2:row.len){

      count.mat[x[i-1], x[i]] <- count.mat[x[i-1], x[i]] + 1

   }
return(count.mat)
}

