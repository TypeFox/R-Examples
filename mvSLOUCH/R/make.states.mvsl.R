## Function from the slouch package. Author : Jason Pienaar

`.make.states.mvsl` <-
function(pars.object){
x<-pars.object$Final.states
n<-length(x)
n.states<-length(pars.object$niche.code[,1])

#which nodes are ambiguous

count<-0

ambig<-NA
for(i in 1:n)
 {
  if(length(x[[i]])>=2)
   {
    count=count+1;
    ambig<-c(ambig,i)	
   }	
 }
ambig=ambig[-1]
n.ambig<-length(ambig)

# choose first character state

tmp<-matrix(data=NA, nrow=n, ncol=1) 
for(i in 1:n)
{
 tmp[i,1]<-x[[i]][1]
  for(j in 1:n.states)
  {
  if (tmp[i,1]==j) tmp[i,1]<-pars.object$niche.code[,1][[j]] 
  }
 }
 
# encode ambiguous states as character ambiguous

if(n.ambig!=0)
 {
  for(i in 1:n)
   {
    for(j in 1:n.ambig)
     {
      if(i==ambig[j]) tmp[i,1]="ambiguous"	
     }
   }
 }

return(as.factor(tmp))
}

