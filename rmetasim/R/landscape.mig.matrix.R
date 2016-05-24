#
# this function produces a 
#

landscape.mig.matrix <- function(h=3,               #habitats
                                 s=2,               # stages
                                 mig.model="island",#also steppingstone.linear,
                                                    #steppingstone.circular, custom,
                                                    #twoD, twoDwDiagonal
                                 first.rep.s=s,     #first reproductive stage
                                 h.dim=NULL,        #habitat map ??
                                 distance.fun=NULL, #gives a density for a distance
                                                    ##and additional params
                                 distance.factor=1, #
                                 R.custom=NULL,...)
  {
    if(s<first.rep.s)
      stop("First reproductive stage must be smaller than total number of stages")
    if(mig.model=="island")
      {
        R.int <- matrix(0,nrow=h,ncol=h)
        R.comb <- combinations(h,s,1:h)[h:1,] #from gtools.  Could try a hand-written one to reduce depend
        for(a in 1:h)
          {
            R.int[R.comb[a,],a] <- 1          
          }
      }

     if(mig.model=="stepping.stone.linear")
       {
        R.int <- matrix(0,nrow=h+2,ncol=h+2)
        for(a in 2:(h+1))
          {
          R.int[c(a-1,a+1),a] <- 1
          }
        R.int <- R.int[-c(1,h+2),-c(1,h+2)]
#        h.dim=NULL
        distance.fun=NA
       }

       if(mig.model=="stepping.stone.circular")
       {
        R.int <- matrix(0,nrow=h+2,ncol=h+2)
        for(a in 2:(h+1))
          {
          R.int[c(a-1,a+1),a] <- 1
          }
        R.int <- R.int[-c(1,h+2),-c(1,h+2)]
        R.int[1,h] <- 1
        R.int[h,1] <- 1
       }


       if(mig.model=="custom")
       {
        if(class(R.custom)!="matrix"|dim(R.custom)[1]!=dim(R.custom)[2])
          {
          stop("there is an error in the matrix provided")
          }
        h <- dim(R.custom)[1]
        R.int <- R.custom
       }


       if(mig.model=="twoD"|mig.model=="twoDwDiagonal")
         {
          if(length(h.dim)!=2)
            {
            stop("the habitat map should have length 2")
            }
          if(h!=h.dim[1]*h.dim[2])
            {
            stop("the total number of habitats should be equal to the number of habitat defined in the map")
            }          
          R.int <- matrix(0,nrow=h,ncol=h)
          R.int[-(1:h.dim[1]),1:(dim(R.int)[1]-h.dim[1])] <-
            R.int[-(1:h.dim[1]),1:(dim(R.int)[1]-h.dim[1])] + diag(1,h.dim[1]*(h.dim[2]-1))
          R.int[1:(dim(R.int)[1]-(h.dim[1])),-(1:h.dim[1])] <-
            R.int[1:(dim(R.int)[1]-(h.dim[1])),-(1:h.dim[1])] + diag(1,h.dim[1]*(h.dim[2]-1))
          R.int[-1,-dim(R.int)[1]] <-
            R.int[-1,-dim(R.int)[1]] + diag(c(rep(c(rep(1,h.dim[1]-1),0),h.dim[2]-1),rep(1,h.dim[1]-1)),(h.dim[1]*h.dim[2])-1)
          R.int[-dim(R.int)[1],-1] <-
            R.int[-dim(R.int)[1],-1] + diag(c(rep(c(rep(1,h.dim[1]-1),0),h.dim[2]-1),rep(1,h.dim[1]-1)),(h.dim[1]*h.dim[2])-1)
        }



       if(mig.model=="twoDwDiagonal")
         {
         R.int[-(1:h.dim[1]-1),1:(dim(R.int)[1]-(h.dim[1]-1))] <-
           R.int[-(1:h.dim[1]-1),1:(dim(R.int)[1]-(h.dim[1]-1))] +
             diag(c(rep(c(0,rep(1,h.dim[1]-1)),h.dim[2]-1),0),h.dim[1]*(h.dim[2]-1)+1)

         R.int[1:(dim(R.int)[1]-(h.dim[1]-1)),-(1:h.dim[1]-1)] <-
           R.int[1:(dim(R.int)[1]-(h.dim[1]-1)),-(1:h.dim[1]-1)] +
             diag(c(rep(c(0,rep(1,h.dim[1]-1)),h.dim[2]-1),0),h.dim[1]*(h.dim[2]-1)+1)

         R.int[-(1:(h.dim[1]+1)),1:(dim(R.int)[1]-(h.dim[1]+1))] <-
           R.int[-(1:(h.dim[1]+1)),1:(dim(R.int)[1]-(h.dim[1]+1))] +
             diag(c(rep(c(rep(1,h.dim[1]-1),0),h.dim[2]-2),rep(1,h.dim[1]-1)),h.dim[1]*(h.dim[2]-1)-1)

         R.int[1:(dim(R.int)[1]-(h.dim[1]+1)),-(1:(h.dim[1]+1))] <-
           R.int[1:(dim(R.int)[1]-(h.dim[1]+1)),-(1:(h.dim[1]+1))] +
             diag(c(rep(c(rep(1,h.dim[1]-1),0),h.dim[2]-2),rep(1,h.dim[1]-1)),h.dim[1]*(h.dim[2]-1)-1)                 
         }


       if(mig.model=="distance")
         {
          if(length(h.dim)!=2)
            {
            stop("the habitat map should have length 2")
            }
          if(h!=h.dim[1]*h.dim[2])
            {
            stop("the total number of habitats should be equal to the number of habitat defined in the map")
             }
          if(class(distance.fun)!="function")
            {
            stop("distance.fun should be a function that takes one element (distance) and returns the migration rate")
            }

            h.position <- matrix(0,nrow=h,ncol=2)

            h.position[,1] <- rep(1:h.dim[2],each=h.dim[1])
            h.position[,2] <- rep(1:h.dim[1],h.dim[2])                 
            h.position <- h.position*distance.factor
          
            R.int <- dist(h.position,diag=F,upper=T)
            R.int  <- as.matrix(R.int)
            R.int <- distance.fun(R.int,...)
            R.int[as.logical(diag(1,h))] <- 0            
         }
    
    
    ##this part actually creates the big habitat*stage matrix,
    ##before that it only created the small, habitat-based matrix
  
      rep.s <- first.rep.s:s
      R <- matrix(0,nrow=h*s,ncol=h*s)
      for(a in 1:h)
       {
        R[seq(1,h*s,s),split(c(rep(rep.s,h)+rep(seq(0,(s*h)-1,s),each=length(rep.s))),rep(1:h,each=length(rep.s)))[[a]]] <- R.int[,a]
       }
  to.return <- list()
  to.return$R <- R
  to.return$h <- h
  to.return$s <- s
  to.return$first.rep.s <- first.rep.s
  to.return$h.dim <- h.dim
  to.return$mig.model <- mig.model
  to.return$distance.fun <- distance.fun
  to.return$R.int <- R.int
  
  return(to.return)
  }
