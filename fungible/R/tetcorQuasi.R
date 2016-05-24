tetcorQuasi<-function (x, y=NULL)
{

#------------------------------------------------------------#
# Function: quasi tetrachoric                                #
#                                                            #
# Niels Waller, December 2006                                #
# Ulrich, R., & Wirtz, M. (2004). On the correlation of a    #
# naturally and an artificially dichotomized variable.       #
# British Journal of Mathematical and Statistical Psychology,#
# 57, 235-252.                                               #
#                                                            #
#                                                            #
#ARGUMENTS:                                                  #
#  x:  An N x 2 matrix or an N x 1 vector of binary          #
#  responses coded 0/1                                       #
#  y : a vector (if x is N x 1) of binary responses          #
#coded 0/1                                                   #
#----------------------------------------------------------- #

# initialize warning flag
warning.msg<-list()


  if (is.data.frame(x))
      x <- as.matrix(x)
  if (is.data.frame(y))
       y <- as.matrix(y)
  if (!is.matrix(x) && is.null(y))
        stop("supply both x and y or an N x 2 matrix x")
  if(!is.null(y))
     x<-cbind(x,y)
  if(sum(is.na(x))>0)
        stop("Missing values are not allowed")

  nitems<-ncol(x)
  #--Do any items have zero variance?
  zeroVar<-apply(x,2,var)==0
  if(sum(zeroVar)>0){
   labs<-paste("Item",1:nitems,sep=" ")
   cat("SERIOUS ERROR\n")
   baditems<-paste("\nThe following item has no variance: ",labs[zeroVar],sep="")
   cat(baditems,"\n")
   stop()
  }



warning.k<-0
########################
#--MAIN loop begins here
########################

    v1 <- x[,1]
    v2 <- x[,2]

    N <- length(v1)
    v1 <- factor(v1,levels=c("0","1"))
    v2 <- factor(v2,levels=c("0","1"))
    abcd <- table(v1,v2)
  # cat(c(abcd,"\n"))
    
    p00 <- abcd[1,1]/N
    p10 <- abcd[1,2]/N
    p01 <- abcd[2,1]/N
    p11 <- abcd[2,2]/N
    

    

    p0. <- p00 + p01
    p1. <- p10 + p11
    if(p0.==0)p0.<-1
    if(p1.==0)p1.<-1
    
    delta <- qnorm(p00/p0.) - qnorm(p10/p1.)
    
    denominator<-sqrt(delta^2 + (1/(p1.*(1-p1.))))
    delta/denominator

}


 # N<-1000
 #  x<-rep(c(0,1),c(N/2,N/2))
 #  b<-c(rnorm(N/2),rnorm(N/2,3))
 #  cor(x,b)
 #  y<-rep(0,N)
 #  y[b>0]<-1
 #  cor(x,y)
 #  nu(x,y)

