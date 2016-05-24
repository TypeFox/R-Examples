`plotRmean` <-
function(haplo.object,all=TRUE,name.var=" ",keep.rares.rmean=FALSE){

boa.init()
chain.import(haplo.object,keep.rares.chain=keep.rares.rmean)

x<-boa.chain("work")
boa.par()

vec<-as.vector(x[[1]])

if (all==TRUE){
       vars<-colnames(x[[1]])
       i<-1
       
            x11()
       par(mfrow=c(2,2))      
       for(i in 1:length(vars)){  
          if (i%%4==0){ 
             x11()
             par(mfrow=c(2,2))
          }        
              boa.plot.history("mcmc",vars[i])
       }
}else{
    par(mfrow=c(2,2))  
    for(i in 1:length(name.var)){
              if (i%%4==0){ 
             x11()
             par(mfrow=c(2,2))
          }        
              boa.plot.history("mcmc",name.var[i])
       }

                
}


}

