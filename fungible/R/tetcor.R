tetcor<-function (X, y=NULL, BiasCorrect=TRUE,stderror=FALSE,Smooth=TRUE, max.iter=5000) 
{
# add: cparam=NULL to argument list
#------------------------------------------------------------#
# Function: tetcor                                           #
#                                                            #
#  warning messages updated Sept 2 2004                      # 
#                                                            #
# Niels Waller, May 2002                                     #
# Requires library:  R2Cuba from CRAN                        #
#                                                            #
# Compute ML estimate of tetrachoric r                       #
#                                                            #
#ARGUMENTS:                                                  #
#  X :  a matrix of binary responses coded 0/1               #
#  y : if x is a vector then y must be a vector of 0/1       #
#     responses                                              #
#  BiasCorrect : a logical that determines whether to apply  #
#  the Brown & Benedetti bias correction                     #
#  stderror : a logical that determines whether to compute   #
#             the s.e. of the estimate                       # 
#                                                            #
#  Smooth :  a logical that determines whether to smooth the #
#            correlation matrix                              #
#                                                            #
#  max.iter :                                                #
#                                                            #
#                                                            #
# VALUE                                                      #
#   r :   matrix of tetrachoric correlations                 #
#   se : matrix of standard error of r is requested          #
#   Warnings : a list of correlations that did not converge  #
#                                                            #
#                                                            #
#                                                            #
# Initial estimate of r calculated by Divgi's method         #
#                                                            #
# Divgi, D. R. (1979) Calculation of the tetrachoric         #
#  correlation coefficient. Psychometrika, 44, 169--172.     #
#                                                            #
# Bias correction:                                           #
# Brown, M. B. & Benedetti, J. K. (1977). On the mean and    #
#  variance of the tetrachoric correlation coefficient.      #
#  Psychometrika, 42, 347--355.                              #
#                                                            #
# S.E. of r_{tet}                                            #
# Hamdan, M. A. (1970). The equivalence of tetrachoric and   #
#  maximum likelihood estimates of $\rho$ in $2 \times 2$    #
#  tables. Biometrika, 57, 212--215.                         #
#----------------------------------------------------------- #

# initialize warning flag
warning.msg<-list()
 
x <- X

# Update: September 8, 2014
# use 'cuhre' for bivariate integral 
# library(R2Cuba)

  if (is.data.frame(x)) 
      x <- as.matrix(x) 
  if (is.data.frame(y)) 
       y <- as.matrix(y)
  if (!is.matrix(x) && is.null(y)) 
        stop("supply both x and y or a matrix-like x")
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
  
  
  r<-matrix(0,nitems,nitems)
  se <- 99
  if(stderror==TRUE) se<-r
  
warning.k<-0  
########################  
#--MAIN loop begins here
########################
for(iter.row in 2:nitems){
   for(iter.col in 1:(iter.row-1)){
   
   #cat(c("     Working on correlation: ",iter.row, iter.col, "\n"))   
          
    v1 <- x[,iter.row]
    v2 <- x[,iter.col]

    N <- length(v1)
    
#---------Compute start values for rhat by Divgi's method------#
 
    
    v1 <- factor(v1,levels=c("0","1"))
    v2 <- factor(v2,levels=c("0","1"))
    abcd <- table(v1,v2)  
  
    cella <- abcd[2,2]
    cellb <- abcd[2,1]
    cellc <- abcd[1,2]
    celld <- abcd[1,1]
    
  
  
#----Brown and Benedetti bias correction----------#

  if(BiasCorrect==TRUE){
       if(cella==0 & celld==0){
         r[iter.row,iter.col]<--1
         se<-NA
         next
       }
       if(cellc==0 & cellb==0){
         r[iter.row,iter.col]<-1
         se<-NA
          next
       }     
       if(cella==0 & cellb > 0 & cellc > 0 & celld > 0){   # only cell a is 0.00
         cella<-.5;  cellb<- cellb - .5;  cellc<- cellc - .5;  celld<- celld +.5      
       }
       if(cellb==0 & cella > 0 & cellc > 0 & celld > 0){   # only cell b is 0.00
         cellb<-cellb+.5; cella<-cella - .5; cellc<-cellc+.5; celld<-celld-.5                
       }
      if(cellc==0 & cellb > 0 & cella > 0 & celld > 0){   # only cell c is 0.00
         cellc<-cellc+.5; celld<-celld-.5; cella<-cella-.5; cellb<-cellb+.5
      } 
      if(celld==0 & cellb > 0 & cellc > 0 & cella > 0){   # only cell d is 0.00 
         celld <- .5; cellb <- cellb - .5; cellc <- cellc - .5; cella <- cella + .5
      }
   }#-------------end of bias correction----------------------------#
 
    p1 <- (cella+cellb)/N
    p2 <- (cellc+cella)/N
        
    p11 <- cella/N
    p00 <- celld/N 
    p01 <- cellc/N
    p10 <- cellb/N   
    

   
  #-------Correction for non-zero lower asymptotes
  # Currently not working
  # Carroll
  # if(!is.null(cparam)){
  #   gi <- cparam[iter.row]
  #   wi <- 1-gi
  #   gj <- cparam[iter.col]
  #   wj <- 1-gj 
  # 
  #   t00 <- p00/(wi * wj)
  #   
  #   t01 <- (wj*p01 * gj*p00)/(wi*wj)
  #   if(t01 < 0) t01 <-0
  #   
  #   t10 <- (wi*p10 - gi*p00)/(wi*wj)
  #   if(t10 < 0) t10<-0
  #   
  #   t11 <- 1 - t00 - t01 - t10 
  #   
  #   p11 <- t11
  #   p1 <- 1 - ( (1-p1)/(1-gi) )     
  #   p2 <- 1 - ( (1-p2)/(1-gj) )  
  #   
  #   cella <- t11 * N
  #   cellb <- t10 * N
  #   cellc <- t01 * N
  #   celld <- t00 * N
  #     
  #}  # end correction for guessing  
        
    h.x <- qnorm(p1)
    k.y <- qnorm(p2)
       
    h.sign <- sign(h.x)
    k.sign <- sign(k.y)
    hstar <- max(abs(h.x),abs(k.y))
    kstar <- min(abs(h.x),abs(k.y))
    
   h2<-hstar^2
   k2<-kstar^2

   Aa<-.5/(1+(h2+k2)*(.12454 - .27102*( 1 - hstar/sqrt(h2+k2) ) ) )
   Bb<-.5/(1+(h2+k2)*(.82281 - 1.03514*( kstar/sqrt(h2+k2) ) ) )
   Cc<- .07557*hstar + (hstar-kstar)^2 * (.51141/(hstar + 2.05793) -.07557/hstar)
   R<-(cella*celld)/(cellb*cellc)     
   Dd<-h.sign*k.sign*kstar*(.79289 + 4.28981/(1+3.30231*hstar))
   alpha<-Aa+Bb*(-1 + 1/(1+Cc*(log10(R) - Dd)^2))
    
          
   rhat<-cos(pi/(1+R^alpha)) # rhat is start value
   if(is.na(rhat>0))rhat<-0  # if x and y have means of exactly .5 rhat is NaN
     
   if(abs(rhat)==1) {        # integration breaks down if |rhat| = 1
     rhat<-.95*rhat
   }
 
#-----End Divgi's method for start values -----------------#

#-----compute bivariate normal density----------#
bvn<-function(z){
    1/(2*pi*sqrt(1-rhat^2)) * exp(-(z[1]^2+z[2]^2-2*rhat*z[1]*z[2])
     /(2*(1-rhat^2)) )
}


#------------------------
# Divgi sends the readers to Pearson 1900 for 
# the formula for this derivative. This form is
## reported in Kirk 1973, p. 261 Psychometrika
 
#---1st derivative of likelihood (a/N) wrt r------------#
dLdr<-function(r){
      1/(2*pi*sqrt(1-r^2)) * exp( - (h.x^2 + k.y^2 - 2*h.x*k.y*r)/(2*(1-r^2)) )
}


eps<-99
iterations<-0
 
#-----Newton Raphson Loop to improve initial estimate

 while( eps > .00001){
#----integrate over bivariate normal surface to estimate (d/N)


### use 'cuhre' from the R2Cuba in place of adapt 
    adaptOut <- R2Cuba::cuhre(ndim = 2, ncomp = 1, integrand = bvn,
             lower = c(-6, -6), upper = c(h.x, k.y), flags = list(verbose = 0),
             rel.tol = .000001, abs.tol = 0, min.eval = 0, max.eval=max.iter)         
                            
     ep11<-adaptOut$value
 
     firstderiv<-dLdr(rhat)
     
     if(firstderiv <= 0.02) firstderiv <- firstderiv * 5.5 #pull deriv back if too low
     rnew<- rhat - (ep11-p11)/firstderiv
     
     if(abs(rnew) > 1) {
         rnew<-.99
         rhat<-99 # do not stop during this iteration
     }
   
  
   
    eps<-abs(rnew-rhat)
    iterations<-iterations+1 
  
    rhat<-rnew
    if(iterations > max.iter){
       warning.k<-warning.k+1
       warn.msg<-paste("WARNING: Correlation ", iter.row, ", ", iter.col, 
                      " failed to converge!", sep="")
      # cat(warn.msg,"\n")
       warning.msg[warning.k]<-warn.msg
       rnew<-51; 
       ## replace with Pearson
       rhat<-cor(x[,c(iter.row,iter.col)])[1,2]
       eps<-.00000001
      }   
    
}


#-------ML estimate of standard error-(Hamdon, 1970)----#

if(stderror){
  se[iter.row,iter.col]<-1/(N*bvn(c(h.x,k.y)))*  (1/cella + 1/cellb + 1/cellc + 1/celld)^-.5
}  
#-------------------------------------------------------#

  r[iter.row,iter.col]<-rhat 
   }
  }

r<-r+t(r)
diag(r)<-1

if(Smooth){
     ULU<-eigen(r)
     U<-ULU$vectors
     L<-ULU$values

     if(min(L)<=0){  # renorm to make matrix positive definite     
         L[L<=0]<-.0001
         Ltot<-sum(L)
         L<-nitems*L/Ltot
         Lsqrt<-diag(sqrt(L))
         Fload<-U%*%Lsqrt
         r<- Fload %*% t(Fload)
         Dmat<-diag(1/sqrt(diag(r)))
         r<-Dmat%*%r%*%Dmat
     }
}   

if(stderror==FALSE) {
     list(r=r,Warnings=warning.msg)
   }

else if (stderror==TRUE){
       se<-se+t(se)
       diag(se)<-1
       list(r=r, se=se, Warnings=warning.msg)
  } 


}

######### CHECK #########

# ### generate bivariate normal data
# library(MASS)
# rho <- .85
# xy<- mvrnorm(100000, mu = c(0,0), Sigma = matrix(c(1,rho,rho,1),ncol=2))
# 
# ## dichotomize
# p1 <- .7
# p2 <- .1
# xy[,1] <- xy[,1] > p1; xy[,2] <- xy[,2] > p2
# 
# tetcor(x=xy[,1], y=xy[,2], max.iter=5000)



