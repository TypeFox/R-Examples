alphafrontier.2d <-
function(xobs, yobs, type="output", alpha=0.95, add=FALSE, confidence=FALSE, shade=FALSE,...)
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
   
# Representation of the alpha quantile efficiency frontier 
  # 1- Computation of the alpha-quantile score
  res1<-alphascore(xobs,yobs, alpha=alpha)  # score computed on (xobs,yobs)
  
  # Initialisation of the graph and the bound
   if(add){  
    bound<-par()$usr}
    else{
     plot(xobs,yobs,type="n")
     bound<-par()$usr}

    # If the options type="input" and "confidence" a new graph is automatically
    # generated inversing the input and and the output on the y-axis and x-axis 
      if(type=="input" & confidence)
     {plot(yobs,xobs)}    
     
 # Computation of the frontier depending on the choosen direction  
  if(type=="output")
   { 
    y.alqf<-yobs/res1$output
    nods.tri<-sort(xobs,index.return=TRUE)
    y.alqf<-y.alqf[nods.tri$ix]
    x.alqf<-c(0,rep(nods.tri$x[!is.na(y.alqf)],each=2), bound[2])
    y.alqf<-c(0,0,rep(y.alqf[!is.na(y.alqf)],each=2))
    
     if(confidence)
     {
     linspace.x<-as.matrix(c(seq(min(xobs),max(xobs),diff(range(xobs))/1000),bound[2]))
     linspace.y<-as.matrix(c(seq(min(yobs),max(yobs),diff(range(yobs))/1000),max(yobs)))
     
      Nx0i=apply(linspace.x,1,function(x) length(which(xobs<=x)))  
      alphan1x0=alpha-(1.96*(alpha*(1-alpha)/Nx0i)^(0.5))
      alphan2x0=alpha+(1.96*(alpha*(1-alpha)/Nx0i)^(0.5))
      alphaLow=ifelse(alphan1x0<0,min(c(alphan1x0[alphan1x0>0],alpha)),alphan1x0)
      alphaUpp=ifelse(alphan2x0>1,1,alphan2x0)

      resLow<-alphascore(xobs, yobs, linspace.x, linspace.y, alpha=alphaLow)  # score computed on (xobs,yobs)
      resUpp<-alphascore(xobs, yobs, linspace.x, linspace.y, alpha=alphaUpp)  # score computed on (xobs,yobs)
    
      x.Low=linspace.x
      x.Upp=linspace.x
      y.Low<-linspace.y/resLow$output
      y.Upp<-linspace.y/resUpp$output
     } 
   }
   else
   {if(type=="input")
    {x.alqf<-xobs*res1$input
     ind.theta<-order(yobs,round(x.alqf,4))
     x.alqf<-c(0,rep(x.alqf[ind.theta],each=2),bound[2])
     y.alqf<-c(0,0,rep(yobs[ind.theta],each=2))
     
     if(confidence)
     {
     y.alqf=x.alqf
     x.alqf=c(0,0,rep(yobs[ind.theta],each=2)) 
       
     linspace.x<-as.matrix(c(seq(min(xobs),max(xobs),diff(range(xobs))/1000),bound[2]))
     linspace.y<-as.matrix(c(seq(min(yobs),max(yobs),diff(range(yobs))/1000),max(yobs)))
     
      Ny0i=apply(linspace.y,1,function(y) length(which(yobs>=y)))  
      alphan1y0=alpha-(1.96*(alpha*(1-alpha)/Ny0i)^(0.5))
      alphan2y0=alpha+(1.96*(alpha*(1-alpha)/Ny0i)^(0.5))
      alphaLow=ifelse(alphan1y0<0,min(c(alphan1y0[alphan1y0>0],alpha)),alphan1y0)
      alphaUpp=ifelse(alphan2y0>1,1,alphan2y0)


      resLow<-alphascore(xobs, yobs, linspace.x, linspace.y, alpha=alphaLow)  # score computed on (xobs,yobs)
      resUpp<-alphascore(xobs, yobs, linspace.x, linspace.y, alpha=alphaUpp)  # score computed on (xobs,yobs)
    
      y.Low<-linspace.x*resLow$input
      y.Upp<-linspace.x*resUpp$input
            
      x.Low<-linspace.y
      x.Upp<-linspace.y
     }
     
    }
    else
    {# In the hyperbolic direction, the frontier is difficult to draw exactly:
     # That is the reason why we took a sample of 1000 points where we computed
     # the score and then calculates the frontiers on these points
   
     linspace.x<-as.matrix(seq(min(xobs),max(xobs),diff(range(xobs))/1000))
     linspace.y<-as.matrix(seq(min(yobs),max(yobs),diff(range(yobs))/1000))
     
     res2<-alphascore(xobs, yobs, linspace.x, linspace.y, alpha=alpha)  
     x.alqf<-linspace.x*res2$hyper
     y.alqf<-linspace.y/res2$hyper
     ind.gamma<-order(round(x.alqf,10),round(y.alqf,10))
     x.alqf<-c(0,rep(x.alqf[ind.gamma],each=2),bound[2])
     y.alqf<-c(0,0,rep(y.alqf[ind.gamma],each=2))
    }
   }
   
   
  # 3- Representation of the alpha-quantile frontier
#  plot(xobs,yobs,type="n")
    if(confidence)
      {if(shade) polygon(c(x.Low,rev(x.Upp)), c(y.Upp,rev(y.Low)),col='lightgrey',border='lightgrey')
          else
          {
           cl <- match.call()
           mf <- match.call(expand.dots = TRUE)
           m <- match(c("lty"), names(mf), 0L)
           if(m==0)
           {lines(x.Upp,y.Upp,lty=2,...)
           lines(x.Low,y.Low,lty=2,...)}
           else
           {lines(x.Upp,y.Upp,...)
           lines(x.Low,y.Low,...)}           
          }  
    }
    lines(x.alqf,y.alqf,...)
  
    if(!add) 
        {if(type=="input"& confidence)
         {points(yobs,xobs)}
         else
         {points(xobs,yobs)}   
        }     
}

