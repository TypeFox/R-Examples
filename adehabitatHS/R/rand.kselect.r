"rand.kselect" <- function(dudi, factor, weight, nrep=200,
                           alpha=0.05, ewa = FALSE)
{
    ## Verifications
    if (!inherits(dudi, "dudi"))
        stop("Object of class dudi expected")
    if (nrow(dudi$tab) != length(factor))
        stop("The factor should have the same length as the dudi object")
    if (nrow(dudi$tab) != length(weight))
        stop("The vector of weights should have the same length as the dudi object")
    if (!is.vector(weight))
        stop("The weights should be placed in a vector")

    ## Data preparation before the external call to the
    ## C function "permutksel"
    tab<-as.matrix(dudi$tab)
    fac<-as.numeric(factor)
    poids<-as.numeric(weight)
    col1<-ncol(dudi$tab)
    lig1<-nrow(dudi$tab)
    lev1<-nlevels(factor)
    nombreani<-tapply(dudi$tab[,1], factor, length)

    ## The C function "permutksel" permutes the weight and computes
    ## the marginality vectors
    res<-.C("permutksel", as.double(t(tab)),
            as.integer(fac), as.double(poids),
            as.integer(col1), as.integer(lev1), as.integer(lig1),
            double(col1*lev1), double(col1*lev1), double(col1*lev1),
            as.integer(nombreani), as.integer(nrep), double(1),
            double(nrep), double(lev1), double(lev1*nrep), double(col1),
            double(nrep*col1*lev1), as.double(dudi$cw), as.integer(ewa),
            PACKAGE="adehabitatHS")

    names(res)<-c("tab", "fac", "weight", "nh", "na",
                  "nl", "use", "ava", "mar",
                  "nbani", "npermut", "obseig",
                  "simeig", "obsmarg", "simmarg", "eigenvp",
                  "simtout", "poco", "ewa")


    ## Marginality observed for each animal ad each variable
    obsval<-matrix(res$mar, ncol=lev1)
    margs<-matrix(res$simtout, nrow=nrep, byrow=FALSE)
    pval1<-obsval
    pval2<-obsval
    pval<-obsval
    sign<-obsval

    ### Marginality for each animal ad each variable
    for (i in 1:nrow(pval)) {
        for (j in 1:ncol(pval)) {
            k<-(i-1)*lev1+j
            pval1[i,j]<-as.randtest(margs[,k], obsval[i,j])$pvalue
            pval2[i,j]<-as.randtest(-margs[,k], -obsval[i,j])$pvalue
        }
    }
    for (i in 1:nrow(pval)) {
        for (j in 1:ncol(pval)) {
            pval[i,j]<-min(c(pval1[i,j], pval2[i,j]))
        }
    }
    ## Corresponding ranking matrix
    for (i in 1:nrow(pval)) {
        for (j in 1:ncol(pval)) {
            if (obsval[i,j]>=0) {
                sign[i,j]<-"+"
                if (pval[i,j]<(alpha/(2*nrow(pval)*ncol(pval))))
                    sign[i,j]<-"+++"
            } else {
                sign[i,j]<-"-"
                if (pval[i,j]<(alpha/(2*nrow(pval)*ncol(pval))))
                    sign[i,j]<-"---"
            }
        }
    }
    colnames(pval)<-levels(factor)
    row.names(pval)<-colnames(tab)
    colnames(obsval)<-levels(factor)
    row.names(obsval)<-colnames(tab)
    colnames(sign)<-levels(factor)
    row.names(sign)<-colnames(tab)


    ## Marginality per animal
    o<-matrix(res$simmarg, ncol=lev1, byrow=FALSE)
    mat<-matrix(0, nrow=lev1, ncol=2)
    colnames(mat)<-c("observed", "pvalue")
    mat<-as.data.frame(mat)
    for (i in 1:lev1) {
        mat[i,1]<-res$obsmarg[i]
        tmp<-as.randtest(o[,i], res$obsmarg[i])
        mat[i,2]<-tmp$pvalue
    }
    row.names(mat)<-levels(factor)

    ## Significance of the first axis of the K-select
    global<-c(0,0)
    names(global)<-c("observed","pvalue")
    global[1]<-res$obseig
    global[2]<-as.randtest(res$simeig, res$obseig)$pvalue


    ## output
    lili<-list(obsval=obsval, pvalue=pval, signification=sign)
    so<-list(global=global, marg=mat, per.ind=lili, alpha=alpha)
    class(so)<-"rand.kselect"
    return(so)
}

