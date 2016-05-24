"kselect" <- function(dudi, factor, weight, scannf = TRUE, nf = 2, ewa = FALSE)
{

    ## 1. Verifications
    if (!inherits(dudi, "dudi")) stop("Object of class dudi expected")
    X<-dudi$tab
    f<-factor
    ab<-weight
    if (nrow(X) != length(f))
        stop("factor should have the same number of observations as dudi")
    if (nrow(X) != length(ab))
        stop("weight should have the same number of observations as dudi")
    if (!is.vector(weight))
        stop("weight should be placed in a vector")
    if (!is.factor(f)) f<-factor(f)


    ## Split the table into a list of tables (one per animal) giving
    ## the values of variables (columns) in each pixel (rows) of the
    ## home range of the animals. Idem for weight
    lo<-split(X,f)
    ab<-split(ab,f)

    ## The weight given to the animals in the analysis
    if (!ewa) {
        poco<-unlist(lapply(ab, function(x) sum(x)/sum(weight)))
    } else {
        poco<-rep(1/nlevels(f), nlevels(f))
    }
    ab<-lapply(ab, function(x) x/sum(x))

    ## Computation of the coordinates of the centroids "available"
    ## in the ecological space
    m<-data.frame(lapply(lo, function(x) apply(x,2,mean)))

    ## Computation of the coordinates of the centroids "used"
    ## in the ecological space
    n<-list()
    for (i in 1:length(lo)) {
        w<-ab[[i]]
        D<-lo[[i]]
        n[[names(lo)[i]]]<-apply(D,2,function(x) sum(w*x))
    }
    n<-data.frame(n)

    ## Note that the data frames are such that variables are in rows
    ## and animals in columns

    ## Analysis: a non centered PCA of the difference between use and
    ## available centroids
    z<-as.dudi(df=n-m, col.w=poco,
               row.w=dudi$cw, call=match.call(), type="kselect",
               scannf = scannf, nf = nf)

    ## The output
    z$initab<-dudi$tab
    z$initfac<-factor
    z$initwei<-weight

    ## Coordinates of the PCA axes on the K-select axes
    U <- as.matrix(z$l1) * unlist(z$lw)
    U <- data.frame(t(as.matrix(dudi$c1)) %*% U)
    row.names(U) <- names(dudi$li)
    names(U) <- names(z$li)
    z$as <- U

    ## Coordinates of the RUs on the K-select axes
    U <- as.matrix(z$l1 * z$lw)
    z$ls <- as.matrix(z$initab) %*% U

    ## Mean used and mean available
    liani <- split(as.data.frame(z$ls), z$initfac)
    liwei <- split(z$initwei, z$initfac)
    mav <- as.data.frame(t(as.matrix(data.frame(lapply(liani,
        function(x) apply(x, 2, mean))))))
    names(mav) <- names(z$li)
    row.names(mav) <- names(z$tab)
    mutemp <- list()
    for (i in 1:length(liwei)) mutemp[[i]] <- apply(liani[[i]],
        2, function(x) weighted.mean(x, liwei[[i]]))
    mut <- as.data.frame(t(as.matrix(data.frame(mutemp))))
    names(mut) <- names(z$li)
    row.names(mut) <- names(z$tab)
    z$mus <- mut
    z$mav <- mav

    ## Output
    return(z)
}

