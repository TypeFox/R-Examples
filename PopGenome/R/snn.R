
snn <- function(bial,populations){

npops  <- length(populations)
groups <- numeric(npops)
  
for(xx in 1:npops){

 what         <- populations[[xx]]
 groups[what] <- xx 

}

whole_pop_size <- length(unlist(populations)) 

if(whole_pop_size<2){return(NaN)}

perm           <- combn(whole_pop_size,2)

MAT  <<- matrix(,whole_pop_size,whole_pop_size)

apply(perm,2,function(x){
        seq1 <- bial[x[1],]
        seq2 <- bial[x[2],]
        diff <- sum(seq1!=seq2,na.rm=TRUE)
	MAT[x[1],x[2]] <<- diff
        MAT[x[2],x[1]] <<- diff 
})

MAT <- MAT 

# calculate nearest neighbours
GLB <- new.env()

GLB$seqq <- 1
res <- apply(MAT,1,function(x){
       nearest        <- which(x == min(x,na.rm=TRUE))
       member         <- groups[GLB$seqq]
       GLB$seqq       <- GLB$seqq + 1 
       Wk             <- sum(groups[nearest]==member,na.rm=TRUE)
       Tk             <- length(nearest)
return(cbind(Wk,Tk))

})

Xk  <- res[1,]/res[2,]
Snn <- sum(Xk)/whole_pop_size

return(Snn)

}
