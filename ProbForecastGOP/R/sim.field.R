"sim.field" <-
function(variog.model,param.est,x,y,n.sim){    
     model <- variog.model
     nugget <- param.est[1]
     variance <- param.est[2]
     scale <-  param.est[3]
     
     mean <- 0
     n.row.f <- length(x)
     n.col.f <- n.sim
     
     f.sim <- matrix(NA,nrow=n.row.f,ncol=n.sim)
     ifelse((variog.model=="gencauchy" | variog.model=="whittlematern"),
       req.par <-c(mean,variance,nugget,scale,param.est[-seq(1:3)]),req.par<-c(mean,variance,nugget,scale))

     if(variog.model=="exponential"){
      for(i in 1:n.sim){
        f.sim[,i] <- GaussRF(x=x,y=y,model=model,param=req.par,grid=FALSE)
      }
     }

     if(variog.model!="exponential"){
      for(i in 1:n.sim){
          f.sim[,i] <- GaussRF(x=x, y=y, model=model,param=req.par,grid=FALSE, method="TBM2")
      }
     }
     return(f.sim)      
}
