opt1D <-
    function(nsim=50,nprocessors=1,setpen="L1",cl=NULL,...){
        if (!(identical(setpen,"L1")|identical(setpen,"L2"))) stop("setpen must be L1 or L2")
        clusterIsSet <- "cluster" %in% class(cl)
        if(nprocessors>1 | clusterIsSet){
            if(!clusterIsSet){
                nprocessors <- as.integer(round(nprocessors))
                cl <- makeCluster(nprocessors, type="PSOCK")
            }
            clusterSetRNGStream(cl, iseed=NULL)
            if(identical(setpen,"L1")){
                thisopt <- parLapply(cl,1:nsim,function(n,...){
                    optL1(...)
                },...)
            }else{   ##if(identical(setpen,"L1")){
                thisopt <- parLapply(cl,1:nsim,function(n,...){
                    optL2(...)
                },...)
            }
            if(!clusterIsSet){
                stopCluster(cl)
            }
        }else{  ##if(nprocessors>1){
            if(identical(setpen,"L1")){
                thisopt <- lapply(1:nsim,function(n,...) optL1(...),...)
            }else{
                thisopt <- lapply(1:nsim,function(n,...) optL2(...),...)
            }
        }
        output <- sapply(thisopt,function(x){
            coefs <- coefficients(x$fullfit,"all")
            tmp <- c(x$lambda,x$cvl,coefs)
            names(tmp) <- c(setpen,"cvl",names(coefs))
            return(tmp)
        })
        return(t(output))
    }

