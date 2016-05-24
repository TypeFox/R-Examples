ordermfrontier.2d <-
function(xobs, yobs, type="output", m=30, add=FALSE, confidence=FALSE, shade=FALSE, ...)
{
 # verification
 n1<-nrow(xobs)    # number of observations
 p.input<-ncol(xobs)     # number of inputs
 q.output<-ncol(yobs)   # number of outputs
 match.arg(type,c("output","input","hyper"))
 if(p.input>1|q.output>1) stop("input and output must be in one dimension only")

  if(type=="hyper")
  {confidence=FALSE
   shade=FALSE}
   
 # Representation of the orderm 
  # 1- Computation of the orderm score
  res1<-ordermscore(xobs, yobs, m=m)  # score computed on (xobs,yobs)

   # Initialisation of the graph and the bound
   if(add){  
    bound<-par()$usr}
    else{
     plot(xobs,yobs)
     bound<-par()$usr}
    
    # If the options type="input" and "confidence" a new graph is automatically
    # generated inversing the input and and the outpur on the y-axis and x-axis 
      if(type=="input" & confidence)
     {plot(yobs,xobs)}    
     
 # Computation of the frontier depending on the choosen direction  
  if(type=="output")
   {y.orderm<-yobs/res1[,1]
    nods.tri<-sort(xobs,index.return=TRUE)
    y.orderm<-y.orderm[nods.tri$ix]
    x.orderm<-c(0,rep(nods.tri$x[!is.na(y.orderm)],each=2), bound[2])
    y.orderm<-c(0,0,rep(y.orderm[!is.na(y.orderm)],each=2))

     if(confidence)
      {linspace.x<-as.matrix(c(seq(min(xobs),max(xobs),diff(range(xobs))/1000),bound[2]))
       linspace.y<-as.matrix(c(seq(min(yobs),max(yobs),diff(range(yobs))/1000),max(yobs))) 
       res2<-ordermscore(xobs, yobs, linspace.x, linspace.y, m=m)
      
       interval=1.96*(apply(linspace.x,1,sigma2hat.output,xobs,yobs,30)/length(xobs))^(1/2)
       x.Low=linspace.x
       x.Upp=linspace.x
       y.Low=linspace.y/res2$output-interval
       y.Upp=linspace.y/res2$output+interval
      } 
   }
   else
   {if(type=="input")
    {x.orderm<-xobs*res1[,2]
     ind.input<-order(yobs,round(x.orderm,4))
     x.orderm<-c(0,rep(x.orderm[ind.input],each=2),bound[2])
     y.orderm<-c(0,0,rep(yobs[ind.input],each=2))        
     
     if(confidence)
      {
       y.orderm=x.orderm
       x.orderm=c(0,0,rep(yobs[ind.input],each=2))  
       
       linspace.x<-as.matrix(c(seq(min(xobs),max(xobs),diff(range(xobs))/1000),bound[2]))
       linspace.y<-as.matrix(c(seq(min(yobs),max(yobs),diff(range(yobs))/1000),max(yobs))) 
       res2<-ordermscore(xobs, yobs, linspace.x, linspace.y, m=m)
      
       interval=1.96*(apply(linspace.y,1,sigma2hat.input,xobs,yobs,30)/length(xobs))^(1/2)
       
       y.Low=linspace.x*res2[,2]-interval
       y.Upp=linspace.x*res2[,2]+interval
       x.Low=linspace.y
       x.Upp=linspace.y
      } 
            
    }
    else
    {# the frontier in the hyperbolic direction is quite difficult to represent
     # by using only the reference points... That is why we create a sample of size 1000
     # equally distanced to compute the values of frontier for these points
    
     linspace.x<-as.matrix(seq(min(xobs),max(xobs),diff(range(xobs))/1000))
     linspace.y<-as.matrix(seq(min(yobs),max(yobs),diff(range(yobs))/1000))
     
     res2<-ordermscore(xobs, yobs, linspace.x, linspace.y, m=m)
     x.orderm<-linspace.x/res2[,3]
     y.orderm<-linspace.y*res2[,3]
     ind.hyper<-order(round(x.orderm,10),round(y.orderm,10))
     x.orderm<-c(0,rep(x.orderm[ind.hyper],each=2),bound[2])
     y.orderm<-c(0,0,rep(y.orderm[ind.hyper],each=2))
    }
   }
   
   
  # 3- Representation of the order-m frontier
   if(confidence)
     {  if(shade) polygon(c(x.Low,rev(x.Upp)), c(y.Upp,rev(y.Low)),col='lightgrey',border='lightgrey')
          else
          {
           cl <- match.call()
           mf <- match.call(expand.dots = TRUE)
           mtest <- match(c("lty"), names(mf), 0L)
           if(mtest==0)
           {lines(x.Upp,y.Upp,lty=2,...)
           lines(x.Low,y.Low,lty=2,...)}
           else
           {lines(x.Upp,y.Upp,...)
           lines(x.Low,y.Low,...)}           
          } 
      }
      
    lines(x.orderm,y.orderm,...)
   if(!add) 
   {if(type=="input"& confidence)
     {points(yobs,xobs)}
     else
     {points(xobs,yobs)}   
    }
 }

 # function which computes the confidence interval of the order-m frontier
 # in the case of output 
sigma2hat.output<-function(x,xobs,yobs,m)
{
   indx=which(xobs<=x)
   Nx=length(indx) 
   Yx=yobs[indx]
   Yx=sort(Yx)
   
    if(Nx<2)
     {
      if(Nx<=1) return(0)
      else 
      print(Nx);return(sum(0.5*(1:(Nx-1)/Nx)^((2*m)-1)*(1-(1:(Nx-1)/Nx))*((Yx[2:Nx]-Yx[1:(Nx-1)])^2)))
     }
     else
     {
      cvect=0.5*(1:(Nx-1)/Nx)^((2*m)-1)*(1-(1:(Nx-1)/Nx))*((Yx[2:Nx]-Yx[1:(Nx-1)])^2)
      a=c(0,((1:(Nx-2)/Nx)^m)*(Yx[2:(Nx-1)]-Yx[1:(Nx-2)]))
      b=cumsum(a)*((1:(Nx-1)/Nx)^(m-1)*(1-(1:(Nx-1)/Nx))*(Yx[2:Nx]-Yx[1:(Nx-1)]))
      return(2*(m^2)*length(xobs)*sum(b,cvect)/Nx) 
     }
     
}

 # function which computes the confidence interval of the order-m frontier
 # in the case of input
  
sigma2hat.input<-function(y,xobs,yobs,m)
{
   indy=which(yobs<=y)
   Ny=length(indy) 
   Xx=xobs[indy]
   Xx=sort(Xx)
   
    if(Ny<2)
     {
      if(Ny<=1) return(0)
      else 
      print(Ny);return(sum(0.5*(1:(Ny-1)/Ny)^((2*m)-1)*(1-(1:(Ny-1)/Ny))*((Xx[2:Ny]-Xx[1:(Ny-1)])^2)))
     }
     else
     {
      cvect=0.5*(1:(Ny-1)/Ny)^((2*m)-1)*(1-(1:(Ny-1)/Ny))*((Xx[2:Ny]-Xx[1:(Ny-1)])^2)
      a=c(0,((1:(Ny-2)/Ny)^m)*(Xx[2:(Ny-1)]-Xx[1:(Ny-2)]))
      b=cumsum(a)*((1:(Ny-1)/Ny)^(m-1)*(1-(1:(Ny-1)/Ny))*(Xx[2:Ny]-Xx[1:(Ny-1)]))
      return(2*(m^2)*length(yobs)*sum(b,cvect)/Ny) 
     }
     
}
