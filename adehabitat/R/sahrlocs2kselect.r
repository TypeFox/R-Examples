"sahrlocs2kselect" <- function(sahr)
{
    ## Verifications
    if (!inherits(sahr, "sahrlocs"))
        stop("sahr should be of class sahrlocs")

    ## Converts the kasc of the study area as a matrix of numbers
    ## and stores the levels of factors in lev
    sa<-sahr$sa
    fac <- unlist(lapply(sa, is.factor))
    if (any(fac)) {
        lev<-list()
        for (i in names(fac[fac]))
            lev[[i]]<-levels(sa[[i]])
        sa<-lapply(sa, as.numeric)
        sa<-as.data.frame(sa)
    }
    sa<-as.matrix(sa)

    ## Converts the kasc of the home ranges as a matrix of numbers
    hr<-sahr$hr
    hr<-lapply(hr, as.numeric)
    hr<-as.matrix(as.data.frame(hr))

    ## Converts the kasc of the relocations as a matrix of numbers
    locs<-as.matrix(sahr$locs)

    ## Replace the missing values by -999
    sa[is.na(sa)]<-(-9999)
    hr[is.na(hr)]<-(-9999)

    ## Use of the C function nls2k to compute the number of lines of
    ## the table in the output
    nh<-ncol(sa)
    np<-nrow(sa)
    na<-ncol(hr)

    so1<-.C("nls2k", as.double(t(sa)), as.double(t(hr)), as.integer(nh),
            as.integer(np), as.integer(na), PACKAGE="adehabitat")[[5]]

    ## Use of this number of lines to reserve memory, and computes the
    ## tables and vectors to be analysed by the K-select
    so2<-.C("sahr2ksel", as.double(t(sa)), as.double(t(hr)),
            as.double(t(locs)), as.integer(nh), as.integer(np),
            as.integer(na), as.integer(so1), double (so1*nh),
            integer(so1), double(so1), PACKAGE="adehabitat")

    ## ta contains the table giving the concatenated tables giving
    ## the values of environmental variables (columns) in each pixel of
    ## the home ranges
    ta<-so2[[8]]
    ta<-as.data.frame(matrix(ta, ncol=nh, byrow=TRUE))
    names(ta)<-names(sahr$sa)

    ## retransform into factors when convenient
    if (any(fac>0)) {
        for (i in names(lev))
            ta[,i]<-factor(ta[[i]],
                           levels=c(1:length(lev[[i]])),
                           labels=lev[[i]])
    }

    ## factor and weight contain respectively a factor giving the ID of
    ## the animals for each row of ta, and the number of relocations in
    ## each pixel in the rows of ta
    factor<-so2[[9]]
    weight<-so2[[10]]
    factor<-factor(factor, labels=names(sahr$hr))

    ## output
    sorties<-list(tab=ta, factor=factor, weight=weight)
    return(sorties)
}

