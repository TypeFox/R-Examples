Qualitative=function(z,model,p,nrr,p2=NULL,alpha=NULL,t=NULL,k=NULL,mm=NULL,pm=NULL,z2=NULL){
  
  if(!is.vector(z)){stop("z must be a vector.")}
  if(any(is.na(z))){stop("There are missing values in z.")}  
  
  if((model<1)|(model>14)){stop("The value of model must be between 1 and 14.")}                         
 
   if((model!=7)&(model!=8)){ 
    if(any((z!=0)&(z!=1))){stop("There are invalid values in z.")} 
  } 
  
  if (model!=8){
    if((p<0)|(p>1)){stop("There are invalid values in p.")}
  }
  
  if((nrr!=1)&(nrr!=2)){stop("The value of nrr must be 1 or 2.")}
  
  if(nrr==1){
    if (model==1){
     r=(z-(1-p))/(2*p-1)
    }
    if ((model==2)|(model==3)){
     if((alpha<0)|(alpha>1)){stop("There are invalid values in alpha.")}
        
     r=(z-(1-p)*alpha)/p
    } 
    if (model==4){
     r=(z-(1-p))/p
    }
    if (model==5){
     if((t<0)|(t>1)){stop("There are invalid values in t.")}
        
     r=(z-(1-t)*(1-p))/(t+(1-t)*(2*p-1))
    } 
    if(model==6){
     if((p2<0)|(p2>1)){stop("There are invalid values in p2.")}
        
     p1=p
     r=(z-p1)/(1-p1-p2)   
    }
    if(model==7){
     if(k<1){stop("k must be a positive number.")}
     if(any((z<1)|(z>k))){stop("There are invalid values in z.")}
     
     if((p2<0)|(p2>1)){stop("There are invalid values in p2.")}
         
     p1=p
     fi=z
     r=((fi/k)-p2)/(p1-p2)
     b=(1-p1-p2)/(k*(p1-p2))
     c=(p2*(1-p2))/(k*(p1-p2)^2)
     vr=b*r+c
    }
    if(model==8){
     if(!is.vector(mm)){stop("mm must be a vector.")}
     if(any(is.na(mm))){stop("There are missing values in mm.")}  
     m=length(mm)
       
     if(!is.vector(pm)){stop("pm must be a vector.")}
     if(any(is.na(pm))){stop("There are missing values in pm.")}
     if(m!=length(pm)){stop("The lengths of mm and pm are different.")}
     if(sum(pm)!=1){stop("The sum of the probabilities pm must be 1.")}
     if(m<=2){stop("The length of the vectors mm and pm must be a positive number greater than 2.")}
     
     if(any((z<1)|(z>m))){stop("There are invalid values in z.")}
         
     mu=sum(mm*pm)
     var=sum(mm^2*pm)-mu^2
     r=(z-mu)/(m+1-2*mu)
     vr=var/(m+1-2*mu)^2
     vr=rep(vr,length(z))
    }
    if(model==9){
     r=(z-(1-p))/((2*p-1)+p*(1-p))
    }
    if(model==10){
     if((alpha<0)|(alpha>1)){stop("There are invalid values in alpha.")}
     if((t<0)|(t>1)){stop("There are invalid values in t.")}
        
     r=(z-((1-t)*(1-p)*alpha))/(t+(1-t)*p)
    }
    if(model==11){
     if((alpha<0)|(alpha>1)){stop("There are invalid values in alpha.")}
     
     r=(z-(1-p)*alpha)/(1-(1-p)*alpha)
    }
  }
  if(nrr==2){
    if(length(z)!=length(z2)){stop("The lengths of z and z2 are different.")} 
    
    if(!is.vector(z2)){stop("z2 must be a vector.")}
    if(any(is.na(z2))){stop("There are missing values in z2.")}
    if(any((z2!=0)&(z2!=1))){stop("There are invalid values in z2.")}
    
    if((p2<0)|(p2>1)){stop("There are invalid values in p2.")}
    
    if ((model==12)|(model==13)|(model==14)){
     p1=p
     r=((1-p2)*z-(1-p1)*z2)/(p1-p2)
    }
  }
  
  if((model!=7)&(model!=8)){
    vr=r*(r-1)
  }
  
  if(any(vr<0)){warning("The transformed variance estimation contains negative values.")}
   
   out=list(TransformedVariable=r,TransformedVariance=vr)
   return(out)
}