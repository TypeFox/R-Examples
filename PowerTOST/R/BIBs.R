# Incomplete block designs (period balanced) 
# from Chow, Liu
# chapter 2.6

bib.CL <- function(trt,p)
{
  if (trt<3 | trt>5) stop("BIB's only for trt=3,4,5 implemented.")
  if (p>=trt | p<=1) stop("p has to be >1 and <",trt," for trt=",trt,".")
  
  if (trt==3){
    # p=2
    # not given in Chow, Liu
    # in the same spirit as for trt>3
    dat <- c(1,2,3,
             2,3,1)
    # variant with is also carry-over balanced
#     dat <- c(1,2,3,3,1,2
#              2,3,1,2,3,1)
    return(matrix(dat, ncol=2))
  }
  if (trt==4){
    if(p==2){
      dat <- c(1,2,3,4,1,2,4,3,1,4,3,2,
               2,3,4,1,3,4,2,1,4,3,2,1)
    }
    if (p==3){
      dat <- c(2,3,4,1,
               3,4,1,2,
               4,1,2,3)
    }
    return(matrix(dat,ncol=p))
  }
  if (trt==5){
    if(p==2){
      dat <- c(1,2,3,4,5,1,3,5,2,4,
               2,3,4,5,1,3,5,2,4,1)
    }
    if(p==3){
      dat <- c(3,4,5,1,2,2,4,1,3,5,
               4,5,1,2,3,4,1,3,5,2,
               5,1,2,3,4,5,2,4,1,3)
    }
    if(p==4){
      dat <- c(2,3,4,5,1,
               3,4,5,1,2,
               4,5,1,2,3,
               5,1,2,3,4)
    }
    return(matrix(dat,ncol=p))
  }
}