ordermscore.boot <- function(xobs, yobs, xeval=xobs, yeval=yobs, m=30, B=200, m.move=FALSE)
{
# initialization
 n1<-nrow(xobs)    # number of observations
 p.input<-ncol(xobs)     # number of inputs
 q.output<-ncol(yobs)   # number of outputs
 n2<-nrow(xeval)     # number of evaluation points
 

 # verification of some conditions
if(n1!=nrow(yobs)) stop("xobs and yobs have not the same number of observations")
if(n2!=nrow(yeval)) stop("xeval and yeval have not the same number of observations")

if(m<=0) stop("m must be an integer and positive")

 # matrix of results
 res<-matrix(0,n2,6)
  
 #########################################
 # order-m efficiency scores calculated on each observation
  for (i in 1:n2)
  {
  # 1. Preparation of the lambda score
   # selection of the observations which verify: xobs[ind] <= Xref[i]
    coord.ix<-matrix(rep(xeval[i,],n1),n1,p.input,byrow=TRUE)
    compare.ix<-ifelse(as.matrix(xobs)<=coord.ix,1,0)
    ind<-which(apply(compare.ix,1,sum)==p.input)
    # creation of the vector of the minimum
    min.yj<-apply(yobs[ind,]/matrix(rep(yeval[i,],length(ind)),length(ind),q.output,byrow=TRUE),1,min,na.rm=TRUE)
    Nx<-length(min.yj)
    if(Nx==1) min.yj=rep(min.yj,2)

  # 2. Preparation of the theta score
   # selection of the observations with y[ind2]>=ind[i]
    coord.iy<-matrix(rep(yeval[i,],n1),n1,q.output,byrow=TRUE)
    compare.iy<-ifelse(as.matrix(yobs)>=coord.iy,1,0)
    ind2<-which(apply(compare.iy,1,sum)==q.output)
   # creation of a vector of the maximum
    max.xy<-apply(xobs[ind2,]/matrix(rep(xeval[i,],length(ind2)),length(ind2),p.input,byrow=TRUE),1,max,na.rm=TRUE)
    Ny<-length(max.xy)
    if(Ny==1) max.xy=rep(max.xy,2)

  # 3. Preparation of the gamma score
    vec.gamma<-apply(cbind(apply(matrix(rep(xeval[i,],n1),n1,p.input,byrow=TRUE)/xobs,1,min),
    apply(yobs/matrix(rep(yeval[i,],n1),n1,q.output,byrow=TRUE),1,min)),1,min)
  
  ###########  Bootstrap Loop on m #######################

  # initialisation
  quant.lambda<-matrix(0,B)
  quant.theta<-matrix(0,B)
  quant.gamma<-matrix(0,B)
  
  # If m.move=TRUE, computation of the optimum value of m 
   m.x<-ceiling(ifelse(m.move,Nx^(1/2),m))  
   m.y<-ceiling(ifelse(m.move,Ny^(1/2),m))
   m.z<-ceiling(ifelse(m.move,length(vec.gamma)^(1/2),m))

    if(Nx!=0 & Ny!=0) # conditions on the boundings ...
     {
      for(b in 1:B)
      {# lambda efficiency score
       quant.lambda[b]<-max(sample(min.yj, m.x, replace=TRUE))
       # theta efficiency score
       quant.theta[b]<-min(sample(max.xy, m.y, replace=TRUE))
       # gamma efficiency score
       quant.gamma[b]<-max(sample(vec.gamma, m.z, replace=TRUE))
      }
     }
    else
     {if(Nx!=0 & Ny==0)  # conditions on the boundings ...
      {
       for(b in 1:B)
        {# lambda score
         quant.lambda[b]<-max(sample(min.yj,  m.x, replace=TRUE))
         # gamma efficiency score
         quant.gamma[b]<-min(sample(vec.gamma, m.z, replace=TRUE))
        }
       quant.theta<-NA
      }
      else
       {if(Nx==0 & Ny!=0)   # conditions on the boundings ...
        {
         for(b in 1:B)
         {#theta score
          quant.theta[b]<-min(sample(max.xy, m.y, replace=TRUE))
          # gamma efficiency score
          quant.gamma[b]<-min(sample(vec.gamma, m.z, replace=TRUE))
         }
        quant.lambda<-NA
       }
       else
       {quant.lambda<-NA
        quant.theta<-NA
       }
      }
    }
   res[i,1]<-mean(1/quant.lambda)
   res[i,2]<-sd(as.vector(1/quant.lambda))
   res[i,3]<-mean(quant.theta)
   res[i,4]<-sd(as.vector(quant.theta))   
   res[i,5]<-mean(1/quant.gamma)
   res[i,6]<-sd(as.vector(1/quant.gamma))  
 }
 
 res<-as.data.frame(res)
 colnames(res)<-c("output.mean","sd.output","input.mean","input.sd","hyper.mean","hyper.sd")

 return(res)
}


