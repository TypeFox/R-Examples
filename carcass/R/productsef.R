productsef <-
function(f, k, n, j){
res<-1
for(i in 0:(n-j-1)) res<-res*(1-f*k^i)
res
}

