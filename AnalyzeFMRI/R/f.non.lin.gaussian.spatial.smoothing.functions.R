#functions to apply gaussian spatial smoothing to an array


NonLinearSmoothArray<-function(x,voxdim=c(1,1,1),radius=2,sm=3,mask=NULL)
{
    
    if(!is.array(x)) return("x should be an array")
    if(length(dim(x))!=3 && length(dim(x))!=4) return("array x should be 3D or 4D")
    if(is.null(mask)) mask <- array(1,dim=dim(x)[1:3])

    
    if(length(dim(x))==3)
    {
        d<-dim(x)
        a<-.C("non_lin_gauss_smooth",
              as.single(aperm(x,c(3,2,1))),
              as.integer(d),
              as.single(aperm(mask,c(3,2,1))),
              as.single(radius),
              as.single(sm),
              as.single(voxdim),
              single(length(x)),
              PACKAGE = "AnalyzeFMRI")
        
        a<-array(a[[7]],dim=d[3:1])
        a<-aperm(a,c(3,2,1))

        return(a)
    }
    else
    {
        d<-dim(x)
        a<-.C("temporal_non_lin_gauss_smooth",
              as.single(aperm(x,c(4,3,2,1))),
              as.integer(d),
              as.single(aperm(mask,c(3,2,1))),
              as.single(radius),
              as.single(sm),
              as.single(voxdim),
              res=single(length(x)),
              PACKAGE = "AnalyzeFMRI")
        
        a<-array(a$res,dim=d[4:1])
        a<-aperm(a,c(4,3,2,1))

        return(a)
    }
    
}    
