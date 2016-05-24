#Signature_2D

Signature_2D=function(out.BUUHWE_2D){

## out.BUUHWT_2D
# "d"  "decomp.hist" "vec.im"      "vec.weights"

edges_select= out.BUUHWE_2D$decomp.hist[1,,]
N1 <- out.BUUHWE_2D$d[1]
N2 <- out.BUUHWE_2D$d[2]
grid= matrix(1:(N1*N2),nrow=N1,byrow=FALSE)

 
signature = matrix(NA,nrow=ncol(edges_select),ncol=6)
colnames(signature)=c('x1','y1','x2','y2','coeff','abs')



for(i in 1:ncol(edges_select)){
BreakpointStart = c(row(grid==edges_select[1,i])[grid==edges_select[1,i]],col(grid==edges_select[1,i])[grid==edges_select[1,i]])
BreakpointEnd = c(row(grid==edges_select[2,i])[grid==edges_select[2,i]],col(grid==edges_select[2,i])[grid==edges_select[2,i]])

## Breakpoint is now the coordinates of the upper part of the edge 
signature[ncol(edges_select)-i+1,]=c(BreakpointStart[1],BreakpointStart[2],BreakpointEnd[1],BreakpointEnd[2],out.BUUHWE_2D$decomp.hist[3,1,i], abs(out.BUUHWE_2D$decomp.hist[3,1,i])) # first row=first rank



}

#signature=rbind(c(1,1,1,1, out.BUUHWT_2D$vec.im[1], abs(out.BUUHWT_2D$vec.im[1])), signature)

return(as.data.frame(signature))

} # e.o.f

#=========================================================================
