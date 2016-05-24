# matrix of all rank patterns
Rpatternmat<-function(nobj)
{
      Rpatt<-perm(nobj,nobj) # permutations from gtools
      # ? RpattStr<-apply(Rpatt,1,function(x) paste(x,sep="",collapse=""))
      # ? ENV$RpattStr<-sort(RpattStr)

      ncomp<-nobj*(nobj-1)/2
      diffs<-matrix(,ncol=ncomp,nrow=nrow(Rpatt))
      c<-0
      for (j in 2:nobj)
          for (i in 1:(j-1) ){
              c<-c+1
              diffs[,c]<-Rpatt[,i]-Rpatt[,j]
          }
      diffs
}
