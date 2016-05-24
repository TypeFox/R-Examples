lambdaoo<-function(kmctime,delta,lambda,gtmat){
.Call('kmcomegalambda', PACKAGE = 'kmc',kmctime,delta,lambda,gtmat)
}

kmcdata_rcpp<-function(kmctime,delta,lambda,gtmat){
    .Call('kmcRCPP_KMCDATA', PACKAGE = 'kmc',kmctime,delta,lambda,gtmat)
}


kmc_routine4<-function( 
		delta, # status 
		lambda,# root finding 
		gtmat # g(t) Matrix p X n 
){ 
    np=as.integer(dim(gtmat)) 
    chk=numeric(np[1]) 
    delta=as.integer(delta) 
    lambda=as.vector(lambda); 
    gtmat=as.numeric(gtmat) 
    #.C extension will ignore variable name! 
    re=.C( 
		"nocopy_kmc_data", delta,gtmat,lambda,np,chk 
    ) 
    return(re[[5]]);
}

kmc_find0loc <- function(d){
 re=1L;
 d=as.integer(d);
 n=as.integer(length(d));
 return(.C("locLastZero",d,n,re)[[3]]);
}
