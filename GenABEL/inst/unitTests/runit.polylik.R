### --- Test setup ---
#
# regression tests
#

if(FALSE) {
    ## Not really needed, but can be handy when writing tests
    library(RUnit)
    library(GenABEL)
}

### do not run
#stop("SKIP THIS TEST")
###

### ---- common functions and data -----

#source(paste("../inst/unitTests/shared_functions.R"))
#source(paste(path,"/shared_functions.R",sep=""))

### --- Test functions ---

test.polylik <- function(tol = 2*.Machine$double.eps^0.5)
{
  	require(GenABEL.data)
    data(ge03d2.clean)
    df <- ge03d2.clean[1:150,autosomal(ge03d2.clean)]
    gkin <- ibs(df,w="freq")
    set.seed(1)
    phdata(df)$height <- phdata(df)$height+rnorm(nids(df),sd=.1*sd(phdata(df)$height,na.rm=TRUE))
    formula <- height ~ sex + age
    mf <- model.frame(formula,phdata(df),na.action=na.omit,drop.unused.levels=TRUE)
    y <- model.response(mf)
    desmat <- model.matrix(formula,mf)
    phids <- rownames(phdata(df))[rownames(phdata(df)) %in% rownames(mf)]
    relmat <- gkin; relmat[upper.tri(relmat)] <- t(relmat)[upper.tri(relmat)]; 
    relmat <- 2.*relmat; relmat <- relmat[phids,phids]
    eigenRO <- eigen(ginv(relmat),symm=TRUE)
    eigenRes <- eigen(relmat,symm=TRUE)
    tO <- system.time(h2htOld <- polygenic(formula,kin=gkin,df,llfun="polylik"))
    tN <- system.time(h2htNew <- polygenic(formula,kin=gkin,df,llfun="polylik_eigen"))
    tO; tN; tO/tN;
    print(h2htNew$h2an)
    print(h2htOld$h2an)
    checkEquals(h2htNew$h2an$est,h2htOld$h2an$est, tolerance = tol)
    checkEquals(h2htNew$h2an$min,h2htOld$h2an$min, tolerance = tol)
    checkEquals(h2htNew$residualY,h2htOld$residualY, tolerance = tol)
    checkEquals(h2htNew$pgresidualY,h2htOld$pgresidualY, tolerance = tol)
    checkEquals(h2htNew$InvSigma,h2htOld$InvSigma, tolerance = tol)
    checkEquals(h2htNew$measuredIDs,h2htOld$measuredIDs, tolerance = tol)
}
