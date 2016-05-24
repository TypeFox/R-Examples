itemParEst<-function(data,model,c=NULL,engine="ltm",discr=1){
if (is.null(c)==TRUE) res<-switch(model,"1PL"=itemPar1PL(data,engine=engine,discr=discr),"2PL"=itemPar2PL(data),"3PL"=itemPar3PL(data))
else res<-itemPar3PLconst(data,c=c)
return(res)}