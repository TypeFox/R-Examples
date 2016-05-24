interp2grid <-
  function(mat,xout,yout,xin=NULL,yin=NULL,type=2) {

    attrib = NULL
######check to ensure matrix
    mat = try(as.matrix(mat))
    if (!is.matrix(mat)) stop('objects must be a matrix')
    if (length(which(is.na(mat)))>0) warning('missing values in matrix can cause interpolation to stop if value used in interpolation')
#######ensure type is 1, 2 or 3
    if (!(type %in% 1:3)) stop('type must be a single numeric value of 1, 2 or 3. See help file')
    

   ###### mat.x = 1:dim(mat)[2]-1;
   ###### mat.y = 1:dim(mat)[1]-1

    mat.x = xin;
    mat.y = yin

    

#######check to ensure all point fall within the boundaries of the matrix
    extents.x = range(mat.x);
    extents.y = range(mat.y)


    
    if (!(min(xout,na.rm=TRUE)>=extents.x[1] & max(xout,na.rm=TRUE)<=extents.x[2]))
      stop('x interpolation data falls outside the input data boundaries')
    if (!(min(yout,na.rm=TRUE)>=extents.y[1] & max(yout,na.rm=TRUE)<=extents.y[2]))
      stop('y interpolation data falls outside the input data boundaries')
    

#######do the interpolation
    out = .Call('interp2grid',mat,mat.x,mat.y,xout,yout,as.integer(type))
    out = matrix(out,nrow=length(yout),byrow=FALSE)
    rownames(out)=yout;colnames(out)=xout

#######return the value
    return(out)
  }

