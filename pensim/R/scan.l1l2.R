scan.l1l2 <-
    function(L1range=c(0.1,100.1),L2range=c(0.1,100.1),L1.ngrid=50,L2.ngrid=50,nprocessors=1,polydegree=1,cl=NULL,...){ #... arguments for cvl
        ##a function for scanning L1 for a given value of L2
        scan.l1 <- function(L2,lambda1vals,...){
            sapply(lambda1vals,
                   function(thisL1,thisL2=L2,...){
                       cvl(lambda1=thisL1,lambda2=thisL2,...)$cvl
                   },...)
        }
        ##set up parallel processing and random number generation
        clusterIsSet <- "cluster" %in% class(cl)
        if(nprocessors>1 | clusterIsSet){
            if(!clusterIsSet){
                nprocessors <- as.integer(round(nprocessors))
                cl <- makeCluster(nprocessors, type="PSOCK")
            }
            clusterSetRNGStream(cl, iseed=NULL)
        }
        ##create the L1 and L2 sequences
        L1vals <- seq(L1range[1]^(1/polydegree),L1range[2]^(1/polydegree),length.out=L1.ngrid)^polydegree
        L2vals <- seq(L2range[1]^(1/polydegree),L2range[2]^(1/polydegree),length.out=L2.ngrid)^polydegree
        ##do the actual work
        if(clusterIsSet | nprocessors > 1){
            ##randomize the order of L1vals and L2vals, so that slow
            ##computations get more evenly distributed across the cluster.  In
            ##some situations where a low value of one of the penalties results
            ##in longer computation time, this should speed up the final
            ##result.
            L1.reorder <- sample(1:length(L1vals),length(L1vals))
            L2.reorder <- sample(1:length(L2vals),length(L2vals))
            L1vals <- L1vals[L1.reorder]
            L2vals <- L2vals[L2.reorder]
            cvl.matrix <- parSapply(cl,L2vals,function(thisL2,...){
                scan.l1(L2=thisL2,lambda1vals=L1vals,...)},
                                    ...)
            ##now re-rorder things
            cvl.matrix <- cvl.matrix[order(L1.reorder),order(L2.reorder)]
            L1vals <- L1vals[order(L1.reorder)]
            L2vals <- L2vals[order(L2.reorder)]
            ##shut down the cluster
            if(!clusterIsSet){
                stopCluster(cl)
            }
        }else{
            cvl.matrix <- sapply(L2vals,function(thisL2,...){
                scan.l1(L2=thisL2,lambda1vals=L1vals,...)},...)
        }
        ##rows correspond to values of L1, columns to values of L2
        rownames(cvl.matrix) <- L1vals
        colnames(cvl.matrix) <- L2vals
        ##return the results in a list
        x <- list(cvl=cvl.matrix,
                  L1range=L1range,L2range=L2range,
                  xlab=paste("L1 (",L1range[1]," to ",L1range[2],")",sep=""),
                  ylab=paste("L2 (",L2range[1]," to ",L2range[2],")",sep=""),
                  zlab=paste("Log-likelihood (",round(min(cvl.matrix),1)," to ",round(max(cvl.matrix),1),")",sep=""),
                  note="rows of cvl correspond to values of lambda1, columns to lambda2")
        return(x)
    }

