#######################################################################
##
## Function: summary.anchors.rank.type()
## Author  : Jonathan Wand <wand(at)stanford.edu>
## Created :  2006-10-02
##
## 
##
## MODIFIED:
##    2007-09-01 : JW
##    - added B
## 
##    2008-04-20 : JW
##    - was summary.anchors()
##    - option to look at B and C simultaneously
##   
#######################################################################
summary.anchors.rank.type <- function(A) {
  N <- NROW(A$span)
  scalar <- apply(A$span,1,function(x) { x[1] == x[2] })
  scalar.freq <- as.data.frame(matrix( 0, nrow=1, ncol=A$max))

  ## make sure that even scalar A that don't occur get included in summary, so loop:
  for (i in 1:A$max) 
    scalar.freq[1,i] <- sum( A$span[scalar,1] == i)
  names(scalar.freq) <- 1:A$max

  ## count:
  n.scalar   <- sum( scalar)
  n.interval <- sum(!scalar)
  avg.span   <- mean( A$span[,2] - A$span[,1] + 1 )
  
  ## summary for scalar cases
  nn <- as.numeric(names(scalar.freq))
  scalar.matrix <- as.data.frame(list(N=as.numeric(scalar.freq),
                                      Prop=as.numeric(scalar.freq/n.scalar),
                                      from=nn,to=nn))
  nn <- names(scalar.freq)
  nn2<- ifelse( nchar(nn)<2, paste(" ",nn,sep=""), nn)
  rownames(scalar.matrix) <- paste(nn2,"to",nn2)

  ## by using Cwmat, we get Cmax length vector
  ## without assuming that every C case is found in scalar or interval span
  ## i.e., this is not same as simply using table()
  ##       but does not require looping
  unif.matrix <- NULL
  if (!is.null(A$weight)) {
    unif.matrix <- as.data.frame(list(N=as.numeric( apply(A$weight,2,sum) ),
                                      Prop = as.numeric( apply(A$weight,2,mean) ),
                                      from= 1:A$max,
                                      to  = 1:A$max))
    rownames(unif.matrix) <- 1:A$max
  }


  #### make full C matrix:
  ## take advantage of the fact that we already calculated the scalar cases
  ## (INCLUDING scalars with zero frequency!)
  nn <- as.numeric(names(scalar.freq))
  interval.matrix <- as.data.frame(list(                           
                           N=as.numeric(scalar.freq),
                           Prop=as.numeric(scalar.freq/N),
                           from=nn,to=nn ))
  nn <- names(scalar.freq)
  nn2<- ifelse( nchar(nn)<2, paste(" ",nn,sep=""), nn)
  rownames(interval.matrix) <- paste(nn2,"to",nn2)
  
  if ( n.interval > 0 ) {
    ## find intervals cases...
    fn <- function(x) {
      paste(
            ifelse(nchar(x[1])<2,paste(" ",x[1],sep=""),x[1]),
            "to",
            ifelse(nchar(x[2])<2,paste(" ",x[2],sep=""),x[2]))
    }
    span <- table(apply( A$span[ A$span[,1] != A$span[,2], ], 1, fn))
    name.span <- names(span)
    n.name.span <- length( name.span )
    ## break apart names produced by table function
    sn <- function(x) {
      as.numeric(unlist(strsplit(x," to ",fixed=TRUE)))
    }
  
    nn.mat <- matrix(NA,nrow=n.name.span,ncol=2)
    for (i in 1:n.name.span ) {
      nn.mat[i,] <- sn(name.span[i])
    }
    zt <- as.data.frame(list(
                             N = as.numeric(span),
                             Prop= as.numeric(span/N),
                             from=nn.mat[,1],to=nn.mat[,2]))
    rownames(zt) <- names(span)
    interval.matrix <- rbind(interval.matrix,zt)
  }

    rval <- list(interval   = interval.matrix,
                 scalar     = scalar.matrix,
                 uniform    = unif.matrix,
                 n.interval = n.interval,
                 n.scalar   = n.scalar,
                 avg.span   = avg.span,
                 max        = A$max)

  return(rval)
}
