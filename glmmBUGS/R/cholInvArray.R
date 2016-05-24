cholInvArray = function(x, prefix="T", chol=FALSE) {


thenames = grep(paste("^", prefix, sep=""), dimnames(x)[[3]], value=TRUE)

 rowInd = gsub(paste("^", prefix, "\\[", sep=""), "", thenames)
 rowInd = gsub(",[[:digit:]+]\\]$", "", rowInd)
 rowInd = as.integer(rowInd)

 colInd = gsub(paste("^", prefix, "\\[[[:digit:]+],", sep=""), "", thenames)
 colInd = gsub("\\]$", "", colInd)
 colInd = as.integer(colInd)

  
 N = max(rowInd)
 Nseq = 1:N

 if( (!all(Nseq %in% colInd) | (!all(Nseq %in% rowInd) ) ) )
  warning("looks like the matrix isn't square")


 for(j in 1:N) {
    print(j)
           
  # Compute L(J,J) and test for non-positive-definiteness.
           
  jm1 = j-1

   #   v = ap::vdotproduct(a.getrow(j, 1, jm1), a.getrow(j, 1, jm1));
   # v = sum(A[k,1:jm1]^2)
   if(jm1>0) {

     v =apply(x[,,paste(prefix, "[", j, ",", 1:jm1, "]", sep=""),drop=FALSE]^2, 
           c(1,2), sum)
   } else {
     v = 0
   }  
         
   ajjString =paste(prefix, "[", j, ",", j, "]", sep="") 
   ajj = x[,,ajjString]-v;
   if( any(ajj<=0) ) {
     warning(paste(j, "row of matrix not pos def"))
    return(ajj)  
 }
             
   x[,,ajjString] = sqrt(ajj)
            
  #          // Compute elements J+1:N of column J.
  if( j<N ) {
    for(i in  seq(j+1,N) ) {
    #                    jm1 = j-1;
    #        v = ap::vdotproduct(a.getrow(i, 1, jm1), a.getrow(j, 1, jm1));\
    #        v = sum(A[i, 1:jm1] * A[j, 1:jm1])
      if(jm1>0) {
            v = apply(x[,,paste(prefix, "[", i, ",", 1:jm1, "]", sep=""),drop=FALSE]  *
                x[,,paste(prefix, "[", j, ",", 1:jm1, "]", sep=""),drop=FALSE], 
           c(1,2), sum)
      } else {
        v = 0
      }  

     #               a(i,j) = a(i,j)-v;
      aijString = paste(prefix, "[", i, ",", j, "]", sep="") 
      x[,,aijString] = (x[,,aijString] - v) / x[,,ajjString]
    }
  }
 } # for loop

 # cholesky finished.  now invert
 if(chol) return(x)
 
 # j is column
 for(j in Nseq) {

    ajjString =paste(prefix, "[", j, ",", j, "]", sep="")
      x[,, ajjString] = 1 / x[,,ajjString]
   #   print(j)

   if(j < N) {
    jp1 = j +1
    # i is row
    for(i in jp1:N) {
      x[,,paste(prefix, "[", i, ",", j, "]", sep="")] =
      (-1 / x[,,paste(prefix, "[", i, ",", i, "]", sep="")]) *
        apply(
         x[,,paste(prefix, "[", j:(i-1), ",",j , "]", sep=""),drop=FALSE] *
          x[,,paste(prefix, "[", i, ",", j:(i-1), "]", sep=""),drop=FALSE],
        c(1,2), sum)
    
    }
   
   }
 }
 

# multpily to create variance matrix

for(j in seq(1, N)) {
  for(i in seq(j, N)) {
     x[,,paste(prefix, "[", i, ",", j, "]", sep="")]  =
      apply(
         x[,,paste(prefix, "[", i:N, ",",i , "]", sep=""),drop=FALSE] *
          x[,,paste(prefix, "[", i:N, ",", j, "]", sep=""),drop=FALSE],
        c(1,2), sum)
  }
} 
 
 
# compute correlations, put them in the upper triangle
# change names of upper triangle
for(i in seq(1, N-1)) {

suffixii = paste(prefix, "[", i, ",", i, "]", sep="")
  for(j in seq(i+1, N)) {
  suffix =  paste(prefix, "[", i, ",", j, "]", sep="")
  suffixji =    paste(prefix, "[", j, ",", i, "]", sep="")
    x[,,suffix]  =
     x[,,paste(prefix, "[", j, ",", i, "]", sep="")] /
       sqrt(   x[,,suffixii] *
         x[,,paste(prefix, "[", j, ",", j, "]", sep="")])
  dimnames(x)[[3]][dimnames(x)[[3]]==suffix] =
    paste("corr", suffixji, sep="")    
  dimnames(x)[[3]][dimnames(x)[[3]]==suffixji] =
    paste("cov", suffixji, sep="")    
  }
  # variance to standard deviations
  x[,,suffixii] = sqrt(x[,,suffixii])
  dimnames(x)[[3]][dimnames(x)[[3]]== suffixii] = 
   paste("sd", suffixii, sep="")
} 
i=N
  suffixii = paste(prefix, "[", i, ",", i, "]", sep="")
  x[,,suffixii] = sqrt(x[,,suffixii])
  dimnames(x)[[3]][dimnames(x)[[3]]== suffixii] = 
   paste("sd", suffixii, sep="")
   
 
return(x)
        
} # function       

