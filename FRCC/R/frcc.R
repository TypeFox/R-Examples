# input:  scaled data matrix  (standardization is NOT checked)
#         weights of each data point
#         target (ignored, assumed to be 0)
#
#
# note: the fast algorithm in this function is due to Miika Ahdesm\"aki
#
off.diagonal.lambda = function(xs, p,q)
#xs must be standarized data, which in our case is performed in frcc but just in case
{  xs<-as.matrix(xs) # standardize input matrix by standard deviations
  # xs = wt.scale(xs, w, center=TRUE, scale=TRUE) # standardize data matrix   #comes from wt.scale.R
  #at this point xs is cbind(X,Y)
  w=matrix(data=c(1),ncol=dim(xs)[1],nrow=1)#all experimental units will be weigthed equally
  # bias correction factors
  w2 = sum(w*w)           # for w=1/n this equals 1/n   where n=dim(xs)[1]
  h1w2 = w2/(1-w2)        # for w=1/n this equals 1/(n-1)

  sw = sqrt(w)
  Q1.squared = (crossprod(sweep(xs, MARGIN=1, STATS=sw, FUN="*")))^2
  Q1.squared_XY = Q1.squared[1:p,p+1:q]#extrating part of XY
  Q2 = crossprod(sweep(xs^2, MARGIN=1, STATS=sw, FUN="*")) - Q1.squared
  Q2_XY = Q2[1:p,p+1:q]#extrating part of XY
  #removing the part about substracting the diagonal because now the target is the null-matrix, not the identity
  denominator = sum(Q1.squared_XY)#-sum(diag(Q1.squared))
  numerator = sum(Q2_XY)#-sum(diag(Q2))

  if(denominator == 0)
    lambda = 1
  else
    lambda = min(1, numerator/denominator * h1w2)

  return (lambda)
}
#===================draws circles========================================================
#from http://www.r-bloggers.com/circle-packing-with-r/ by Michael Bedward
custom.draw.circle <- function(x, y, r, col) {
    lines( cos(seq(0, 2*pi, pi/180)) * r + x, sin(seq(0, 2*pi, pi/180)) * r + y , col=col )
  }
#=====================Plot Units=========================================================
plot_units<-function(X,Y,res.mrcc,i,text_size=.8,point_size=2)
{#libraries needed for the plots
 require("calibrate")
 #horizontal axis for the ith CCA units
 U<- as.matrix(X) %*% as.numeric(res.mrcc$canonical_weights_X[,i])
 U<-(U-mean(U))/as.numeric(sapply(as.data.frame(U),sd))
 #vertical axis for the ith CCA units
 V<- as.matrix(Y) %*% res.mrcc$canonical_weights_Y[,i]
 V<-(V-mean(V))/as.numeric(sapply(as.data.frame(V),sd))
 #------------------Dr. Whitmore's graphs------------------------------------------
 #creating names for the axis
  s1<- paste("U",i,sep="")
  s2<-paste("r=",round(res.mrcc$cor[i], digits = 3),sep="")
  s3<-paste("p-value=",round(res.mrcc$p_value[i], digits = 3))
  my_xlab=paste(s1,s2,s3, sep= "       ")
  my_ylab=paste("V",i)
  plot(U,V, ylim=c(min(V)-.5,max(V)+.5),xlim=c(min(U)-.5,max(U)+.5),col="white",xlab=my_xlab, ylab=my_ylab)
  #Plotting   X cex =size, pch= figure  19=circle, col= color
  points(U,V, cex=point_size,pch=19,col="black")
  textxy(U,V, rownames(X) ,cx=text_size,dcol="black")
}
#=============================Plots Variables============================================
plot_variables<-function(res.mrcc,i,j,inner_circle_radius=.5,text_size=.8)
{#libraries needed for the plots
#for textxy
 require("calibrate")
#Plots
 p <-  nrow(res.mrcc$canonical_factor_loadings_X)
 q <-  nrow(res.mrcc$canonical_factor_loadings_Y)
 first_dimension<-matrix(data=c(0),nrow=p+q,ncol=1)
 first_dimension[1:p,1] <- res.mrcc$canonical_factor_loadings_X[,i]
 first_dimension[(p+1):(p+q),1] <- res.mrcc$canonical_factor_loadings_Y[,i]

 second_dimension<-matrix(data=c(0),nrow=p+q,ncol=1)
 second_dimension[1:p,1] <- res.mrcc$canonical_factor_loadings_X[,j]
 second_dimension[(p+1):(p+q),1] <- res.mrcc$canonical_factor_loadings_Y[,j]
  #creating names for the axis
  my_xlab=paste("CC",i)
  my_ylab=paste("CC",j)
 #no need to find the limits of the plot, they are always 1, because they are correlations
   plot(first_dimension,second_dimension, ylim=c(-1.25,1.25),xlim=c(-1.25,1.25),pch=" ", xlab=my_xlab, ylab=my_ylab)
  #Plotting   X cex =size, pch= figure  17=triangle, col= color ="grey70"
  points(first_dimension[1:p,],second_dimension[1:p,], cex=2,pch=17,col="grey70")
  #Plotting   X cex =size, pch= figure
  points(first_dimension[(p+1):(p+q),],second_dimension[(p+1):(p+q),], cex=2,pch=19,col=1)
  #obtining the names of the variates
 points_names<-c(rownames(res.mrcc$canonical_factor_loadings_X),rownames(res.mrcc$canonical_factor_loadings_Y))
 #deleting the names of the non-significan variables
 for (i in 1:(p+q))
 { if(  sqrt((first_dimension[i,1]^2) + (second_dimension[i,1]^2 ) )< inner_circle_radius )
   {
    points_names[i] <- " "
   }
 }#end for i
  #print(points_names)
 textxy(first_dimension,second_dimension, points_names ,cx=text_size,dcol="blue")
 custom.draw.circle(0,0,inner_circle_radius,col="black")
 custom.draw.circle(0,0,1.0,col="black")
 #Lines across the Axis
 segments(-1.45,0 , 1.45,0 )
 segments(0,-1.45,0 , 1.45 )
}

#==================Generate multivariate Normal Sample==========================================
generate_multivariate_normal_sample<-function(p,q,n)
{
#library needed for the function rmnorm
 require("MASS")
#p = #number of variables in X
#q = #number of variables in Y
#n = #number of observations
#Creating Correlation matrix, this will be the objective of the Sigmas
if( (p > 6) & (q>6))
{
Sigma_Z<-matrix(data=c(0.0),nrow<-(p+q), ncol<-(p+q))   # From Gnz. Journal of Biological Systems17(2):173-199 (2009) page 178
for (i in 1:(p+q))
{  for (j in 1:(q+p))
   { if (i==j)  # Filling the diagonal wih 1's
     {
     #Sigma_Z[i,j]<-1.0
      Sigma_Z[i,j]<- 1.0#*sample(random_number,1)
     }
     if( ((i== (p+1)) & (j==1)) || ((j==(p+1)) & (i==1)) ) #establishing that X1 and Y1 have a correlation of .9
     {Sigma_Z[i,j]<- .9
     }
     if( ((i== (p+2)) & (j==2)) || ((j==(p+2)) & (i==2)) ) #establishing that X2 and Y2 have a correlation of .7
     {Sigma_Z[i,j]<- .7 # 1.4
     }
      if( ((i== (p+3)) & (j==3)) || ((j==(p+3)) & (i==3)) ) #establishing that X2 and Y2 have a correlation of .7
     {Sigma_Z[i,j]<- .5 # 1.4
     }
     if( ((i== (p+4)) & (j==4)) || ((j==(p+4)) & (i==4)) ) #establishing that X2 and Y2 have a correlation of .7
     {Sigma_Z[i,j]<- .3 # 1.4
     }
     if( ((i== (p+5)) & (j==5)) || ((j==(p+5)) & (i==5)) ) #establishing that X2 and Y2 have a correlation of .7
     {Sigma_Z[i,j]<- .1 # 1.4
     }
   }#end for j
} #end for i
#creating random sample from multivariate normal with covariance matrix created above
means<- matrix(data=c(0.0),nrow<-1, ncol<-(p+q))   #notice mean =0
ndata<- mvrnorm(n,Sigma=Sigma_Z,mu=means)
X<-ndata[,1:p]
Y<-ndata[,(p+1):(p+q)]
res<-list(X=matrix(data=c(0),nrow=n,ncol=p),Y=matrix(data=c(0),nrow=n,ncol=q),Sigma_Z=matrix(data=c(0),nrow=(p+q),ncol=(p+q)))
res$X<-X
res$Y<-Y
res$Sigma_Z<- Sigma_Z
res
}#end if
else{print("Not enough variables")}
}#end function
#====================function needed to rearrange the res.frcc structure if th correlations are out of order due to insufficient data ===========================
rearrange.frcc <-function(res.frcc)
{#  print(res.frcc)
 #Dimensions of the problem
p<-dim(res.frcc$canonical_weights_X)[1]
#p = #number of variables in X
q<-dim(res.frcc$canonical_weights_Y)[1]
#q = #number of variables in Y
 aux.res.frcc<-list(cor=0,canonical_weights_X=matrix(data=c(0),nrow=p,ncol=1),canonical_weights_Y=matrix(data=c(0),nrow=q,ncol=1)) 
 for (i in 1:(q-1))
 {
   for (j in (i+1):q)
   {  #we have found a correlation grater than correlation i
      if (res.frcc$cor[j] > res.frcc$cor[i])  
      { #print("#swaping the correlations and the corresponding coefficients (at this point we don't have p-values or canonical factor loading to worry about")
        aux.res.frcc$cor <-  res.frcc$cor[i]
        aux.res.frcc$canonical_weights_X <-  res.frcc$canonical_weights_X[,i]
        aux.res.frcc$canonical_weights_Y <-  res.frcc$canonical_weights_Y[,i]
        
        res.frcc$cor[i]  <-  res.frcc$cor[j]
        res.frcc$canonical_weights_X[,i] <- res.frcc$canonical_weights_X[,j]
        res.frcc$canonical_weights_Y[,i] <- res.frcc$canonical_weights_Y[,j]
        
        res.frcc$cor[j] <- aux.res.frcc$cor 
        res.frcc$canonical_weights_X[,j] <- aux.res.frcc$canonical_weights_X 
        res.frcc$canonical_weights_Y[,j] <- aux.res.frcc$canonical_weights_Y 
     }#end if
   }#end for j now the CCA i is in the correct place
 }#end for i
 res.frcc
} #end rearrange.frcc function

#=========================function which performs tha frccA ======================================================================
frcc<-function(X,Y)
{
#library with the shrink functions from [Schafer and Strimmer, 2005]
 require("corpcor")
#library to calculate p-values
 require("CCP")
#Dimensions of the problem
p<-dim(X)[2]
#p = #number of variables in X
q<-dim(Y)[2]
#q = #number of variables in Y
n<-dim(X)[1] #or =dim(Y)[1]
#n = #number of observations
#not everyting is declared but ut works just fine
 res.frcc<-list(cor=matrix(data=c(0),nrow=1,ncol=q),p_values=matrix(data=c(0),nrow=1,ncol=q),canonical_weights_X=matrix(data=c(0),nrow=p,ncol=q),canonical_weights_Y=matrix(data=c(0),nrow=q,ncol=q))
#just in case that we receive a datasets that has not been standarized:
for (i in 1:p)
{ X[,i]<-(X[,i]-mean(X[,i]))/as.numeric(sapply(as.data.frame(X[,i]),sd))
} #end for i
for (i in 1:q)
{ Y[,i]<-(Y[,i]-mean(Y[,i]))/as.numeric(sapply(as.data.frame(Y[,i]),sd))
} #end for i
#------finding regularized cross-correlation matrix---------------
S_star_XX <- cor.shrink (X)
S_star_YY <- cor.shrink (Y)
lambda_XY<-off.diagonal.lambda(cbind(X,Y), p,q)
#shrinking the off-diagonal matrices
S_star_XY <- cor(X,Y) *(1- lambda_XY)# +lambda_XY*T_XY, but T_XY=null-matrix 
S_star_YX<- t(S_star_XY)
print(paste("Off-Diagonal Shrinkage Coefficient:",lambda_XY))
#--------------------calculating the canonical weights    -------------------------------------------------------
latent_roots<-eigen(solve(S_star_YY)%*%S_star_YX%*%solve(S_star_XX)%*%S_star_XY, only.values = FALSE) #from 1st Eq. page 14 "Understanding..." book
res.frcc$canonical_weights_Y<- latent_roots$vectors
#calculating canonical weights for X
for (i in 1:q)
{
 res.frcc$canonical_weights_X[,i] <- solve(S_star_XX)%*%S_star_XY %*% res.frcc$canonical_weights_Y[,i] 
}#end for i
#Assigning names to the variables
rownames(res.frcc$canonical_weights_X)<-colnames(X)
rownames(res.frcc$canonical_weights_Y)<-colnames(Y)
#------Calculating Canonical Correlations Original Method, not as in Clark book (can you believe that the guy from the CCP packae called them canonical coefficients?)---------------------------   
 for (i in 1:q)
 {
 U<- as.matrix(X) %*% as.numeric(res.frcc$canonical_weights_X[,i])
  U<-(U-mean(U))/as.numeric(sapply(as.data.frame(U),sd))

 #vertical axis for the ith CCA units
 V<- as.matrix(Y) %*% res.frcc$canonical_weights_Y[,i]
 V<-(V-mean(V))/as.numeric(sapply(as.data.frame(V),sd))

 res.frcc$cor[i]<-cor(U,V)
 }#end for i
#obtaining p-values
#------Checking if the correlations are in order, if we have too little data =>they are not 
flag_anomaly_in_correlations<-0
for(i in 1:(q-1))
{ if (res.frcc$cor[i] < res.frcc$cor[i+1])
  {flag_anomaly_in_correlations <-1}
}
if(flag_anomaly_in_correlations > .5)
{
 print("Warning: Too little data, results might be unreliable")
 res.frcc <- rearrange.frcc(res.frcc)
}
#Calculating the p-values
p_values_struct=0
p_values_struct<-p.asym(res.frcc$cor, n, p, q, tstat = "Wilks") #we can use   p.perm if we are not sure about the normality of the data
res.frcc$p_values<- p_values_struct$p.value
#-----------Calculating Canonical Factor loadings (Statistica) or interset correlation coefficients (Gittins)-----------------
canonical_factor_loadings_X <- matrix(data=c(0),nrow=p,ncol=q)    #page  38 and 39 of Gittins' book  (just not the final formula but the2nd after the definition)
canonical_factor_loadings_Y <- matrix(data=c(0),nrow=q,ncol=q)
#S_star_XY<-cor(X,Y)
#S_star_YX<- t(S_star_XY)
for (i in 1:q)
{
  #canonical_factor_loadings_X[,i] <-  S_star_XY %*%  as.matrix(res.frcc$canonical_weights_Y[,i])
   # Our data is already standarized hence our X is z^(x) and Y is z^(y) in Gittins, b= canonical_weights_Y
   #print(dim(as.matrix( res.frcc$canonical_weights_Y[,i] )))
   #print(class(Y))
   #print(dim(as.matrix(Y)))
   canonical_factor_loadings_X[,i] <- cor(as.matrix(X), as.matrix(Y) %*% as.matrix( res.frcc$canonical_weights_Y[,i] ) )
  #canonical_factor_loadings_Y[,i] <-  S_star_YX %*%  as.matrix(res.frcc$canonical_weights_X[,i])
   canonical_factor_loadings_Y[,i] <- cor(as.matrix(Y), as.matrix(X) %*%  as.matrix(res.frcc$canonical_weights_X[,i])  )
}#end for i
res.frcc$canonical_factor_loadings_X <- matrix(data=c(0),nrow=p,ncol=q)    #page  38 and 39 of Gittins' book
res.frcc$canonical_factor_loadings_Y <- matrix(data=c(0),nrow=q,ncol=q)
#As in [Schafer and Strimmer] top of p.10 we need to take care of value s greater than 1
for (i in 1:q)
{
 res.frcc$canonical_factor_loadings_X[,i] <- canonical_factor_loadings_X[,i]
 res.frcc$canonical_factor_loadings_Y[,i] <- canonical_factor_loadings_Y[,i]
}#end for i
#Assigning names to the variables
rownames(res.frcc$canonical_factor_loadings_X)<-colnames(X)
rownames(res.frcc$canonical_factor_loadings_Y)<-colnames(Y)
#Priting the final results
print(res.frcc)
res.frcc
}#end mrc function
