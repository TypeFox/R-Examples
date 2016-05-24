
#I hardly know'er!

#code to get the spatially minimally spanning (convex hull/bounding box) box for two boxes:

#(note that it might not be the minimally spanning on a points basis)

#input:  two boxes (not box lists - though a box is a list)

unioner <- function(box1, box2){

#find which dimensions are in common:  Those are the only ones that can be restricted

commondims <- intersect(box1[[1]],box2[[1]])

b1ci <- c(1:length(box1[[1]]))[box1[[1]]%in%commondims]
b2ci <- c(1:length(box2[[1]]))[box2[[1]]%in%commondims]

#make the new bounds matrix
newbmat <- matrix(ncol = 2, nrow = length(commondims))

cd <- commondims

#fill it in with the extreme values along each commonly restricted dimension
for (i in 1:length(commondims)){

newbmat[i,1] <-  min(box1[[2]][b1ci[i],1],box2[[2]][b2ci[i],1])
newbmat[i,2] <-  max(box1[[2]][b1ci[i],2],box2[[2]][b2ci[i],2])
  
  
}

return(list(commondims,newbmat))

}


  
  