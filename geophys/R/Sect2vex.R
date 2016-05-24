
Sect2vex<-function(V1, V2)
  {
    ###  find intersection point of two 2D vectors
    ####   vectors cannot be co-linear

    ##  inoput can be a list or a matrix of vectors

    if(is.list(V1))
      {
        x1 = V1$x[1]
        y1 = V1$y[1]
        x2 = V1$x[2]
        y2 = V1$y[2]
      }
    if(is.matrix(V1))
      {
        x1 = V1[1,1]
        y1 = V1[1,2]
        x2 = V1[2,1]
        y2 = V1[2,2]
      }

    if(is.list(V2))
      {
        
        x3 = V2$x[1] 
        y3 = V2$y[1]
        x4=  V2$x[2]   
        y4 = V2$y[2]
      }
    if(is.matrix(V2))
      {
        x3 = V2[1,1]
        y3 = V2[1,2]
        x4 = V2[2,1]
        y4 = V2[2,2]
      }




    

   ####   Tparallel = AXB.prod(c(x2-x1, y2-y1) , c(x4-x3, y4-y3) )

    Den = ( (y4-y3)*(x2-x1) - (x4-x3)*(y2-y1 ) )

    
    if(Den==0){
      cat("Parallel Vectors, No Intersection", sep="\n")
      return(NULL)
    }

    ua=(  (x4-x3)*(y1-y3) - (y4-y3)*(x1-x3) )/ Den

    xp = x1+ua*(x2-x1)
    yp = y1 + ua*(y2-y1)

    return(list(x=xp, y=yp))
    
}
