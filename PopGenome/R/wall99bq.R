#Cai JJ (2008) PGEToolbox: A Matlab toolbox for population genetics and evolution
#Journal of Heredity Jul-Aug;99(4):438-40. doi:10.1093/jhered/esm127
#modified


wall99bq   <- function(matrix_pol,populations){

npops      <- length(populations)
popnames   <- paste("pop",1:npops)

segsites   <- get_segsites(matrix_pol,populations) # positions of the segsites of each population

init       <- rep(0,npops)
S          <- init
Q          <- init
B          <- init

names(Q) <- popnames
names(B) <- popnames

for(xx in 1:npops){

    segmatrix    <- matrix_pol[populations[[xx]],segsites[[xx]],drop=FALSE]
    if(dim(segmatrix)[2]<=1){Q[xx]<- NaN;B[xx]<-NaN;next;}
    con          <- .Call("count_congruent", segmatrix)
    #con          <- count_congruent(segmatrix) # neighbouring sites are identical !
    consum       <- sum(con)
    B[xx]        <- consum/(length(segsites[[xx]])-1)

    # Count the types of different congruent pairs
    checkmatrix  <- segmatrix[,con]
    checkmatrix  <- unique(t(checkmatrix))
    Scheck       <- dim(checkmatrix)[1]
    # --------------------------------------------
    Q[xx]        <- (consum + Scheck)/length(segsites[[xx]])
    
}

return(list(Q=Q,B=B))
}# End of function

######### SUBFUNCTIONS #############################################
# ----------------------------------------------------
# How many sitepairs are identical

count_congruent <- function(segmatrix){

back <- vector(,dim(segmatrix)[2]-1)
for(xx in 1:(dim(segmatrix)[2]-1)){
  site1 <- segmatrix[,xx]
  site2 <- segmatrix[,xx+1]
  if(dim((unique(cbind(site1,site2))))[1]==2){back[xx]<- TRUE }else{back[xx] <- FALSE}
  
}
return(back)
}

