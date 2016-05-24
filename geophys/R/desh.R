desh<-function(M, add=TRUE, PTS=TRUE, colmesh=grey(.8), colpts=grey(.5), ... )
  {
#####  draw a meshgrid
    if(missing(add)) { add=TRUE }
    if(missing(PTS)) { PTS= TRUE }
    if(missing(colmesh)) { colmesh=grey(.8) }

   if(missing(colpts)) {  colpts=grey(.5) }
    
    DM = dim(M$x)

    nx =  DM[2]
    ny =  DM[1]
    if(!add)
      {
        plot(M, type='n', ann=FALSE, axes=FALSE, asp=1)
      }

    
    if(PTS)  points(M, col=colpts, ... )
    
    segments(M$x[1:(nx-1),], M$y[1:(nx-1),], M$x[2:nx,], M$y[2:nx,], col=colmesh, ...)

    segments(M$x[,1:(ny-1)], M$y[,1:(ny-1)], M$x[,2:ny], M$y[,2:ny], col=colmesh, ...)
  }
