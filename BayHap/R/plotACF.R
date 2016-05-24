`plotACF` <-
function(x,all=TRUE,name.var=" ",keep.rares.acf=FALSE){

boa.init()
chain.import(x,keep.rares.chain=keep.rares.acf)
boa.par(acf.lags=1)
x<-boa.chain("work")


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
          boa.plot.acf(1,pname=vars[i])
       }
}else{
    par(mfrow=c(2,2))  
    for(i in 1:length(name.var)){
             if (i%%4==0){ 
                 x11()
                 par(mfrow=c(2,2))
             }        
             boa.plot.acf.ad(x,name.var[i])
    }

                
}
}

