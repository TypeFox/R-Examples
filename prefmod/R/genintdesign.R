`genintdesign` <-
function(pattern,sinclair=TRUE)
################################################################
# interactions design matrix
################################################################
#
# generates design matrix for interaction effects between two comparisons
# with one common object
#
{
  nobj  <-get("nobj",get("ENV", environment(patt.design)))
  ncomp<-nobj*(nobj-1)/2          # number of comparisons
  compvec<-matrix(c(0:0),ncomp,2) # matrix to describe which objects in which comparison
  row<-1
  if (sinclair) {                 # (12),(13),(23),(14),(24),...
     for (j in 2:nobj) {
         for (i in 1:(j-1) ){
             compvec[row,1] <- i
             compvec[row,2] <- j
             row <- row + 1
         }
     }
  } else {                        # (12),(13),(14),...,(23),(24),...
     for (i in 1:(nobj-1)) {
         for (j in (i+1):nobj){
             compvec[row,1] <- i
             compvec[row,2] <- j
             row <- row + 1
         }
     }
  }
  vec1<-NULL
  vec2<-NULL
  vec3<-NULL
  desc.names<-NULL
  for (i in 1:(ncomp-1)) {
     for (j in (i+1):ncomp) {
        if (compvec[i,1]==compvec[j,1]) {
           vec1<-c(vec1,i)
           vec2<-c(vec2,j)
           vec3<-c(vec3,+1)  # sign + because common object preferred in both comps
           desc.names<-c(desc.names,paste("I",compvec[i,1],".",compvec[i,2],compvec[j,2],sep="",collapse=""))
        }
        if (compvec[i,2]==compvec[j,2]) {
           vec1<-c(vec1,i)
           vec2<-c(vec2,j)
           vec3<-c(vec3,+1)  # sign + because common object not preferred in both comps
           desc.names<-c(desc.names,paste("I",compvec[i,2],".",compvec[i,1],compvec[j,1],sep="",collapse=""))
        }
        if (compvec[i,1]==compvec[j,2]) {
           vec1<-c(vec1,i)
           vec2<-c(vec2,j)
           vec3<-c(vec3,-1)  # sign - because common object preferred/not preferred
           desc.names<-c(desc.names,paste("I",compvec[i,1],".",compvec[i,2],compvec[j,1],sep="",collapse=""))
        }
        if (compvec[i,2]==compvec[j,1]) {
           vec1<-c(vec1,i)
           vec2<-c(vec2,j)
           vec3<-c(vec3,-1)  # sign - because common object preferred/not preferred
           desc.names<-c(desc.names,paste("I",compvec[i,2],".",compvec[i,1],compvec[j,2],sep="",collapse=""))
        }
     }
  }
  multvec<-cbind(vec1,vec2)        # each row contains column numbers of comparison patterns
                                   #    to be multiplied for interactioneffect


  nintpars<-length(desc.names)    # number of interaction parameters
  assign("nintpars",nintpars,get("ENV", environment(patt.design)))

#pattern<-pattern-median(range(pattern))
#pattern<-pattern*2
#print(pattern)

  designmatrix<-NULL
  for (i in 1:nintpars) {
      work<-pattern[,multvec[i,1]]*pattern[,multvec[i,2]]*vec3[i]
      designmatrix<-cbind(designmatrix,work)
  }

  colnames(designmatrix)<-desc.names    # names for interaction parameters,
                                        #  e.g., I3.45, comparisons 3-4 and 3-5
  designmatrix
}
