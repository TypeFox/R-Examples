# R port of GENECOUNTING/PREPARE
# 29-1-2004 start implementing
# 30-1-2004 in shape
# 31-1-2004 working

pgc <- function (data,handle.miss=1,is.genotype=0,with.id=0)
{
    nobs <- dim(data)[1]
    nloci2 <- dim(data)[2]
    if (is.genotype)
    {
       nloci <- nloci2
       data<-cbind(data,data)
       a1 <- a2 <- gid <- 0
       for (i in 1:nobs)
       {
           row.i <- data[i,]
           for (j in 1:nloci)
           {
               .C("g2a_",s=as.integer(row.i[j]),a1=as.integer(a1),a2=as.integer(a2),gid=as.integer(gid),PACKAGE="gap")
               data[i,2*j-1] <- a1
               data[i,2*j] <- a2
           }
       }
    }
    else nloci <- nloci2/2
    data <- as.matrix(data)
    stack <- rbind(data[,(2*1:nloci)-1],data[,(2*1:nloci)])
    alleles <- apply(stack,2,max)
    idsave <- wt <- array(0,nobs)
    obscom <- nobs
    data <- t(data)
    gret <- matrix(array(0,nobs*nloci2),nrow=nobs)
    z <- .C("pgc",gdata=as.integer(data),handlemiss=as.integer(handle.miss),nobs=as.integer(nobs),nloci=as.integer(nloci),
            alleles=as.integer(alleles), wt=as.integer(wt),gret=as.integer(gret),
            withid=as.integer(with.id),idsave=as.double(idsave),obscom=as.integer(obscom),PACKAGE="gap")
    subset <- 1:(z$obscom)
    gret <- matrix(z$gret,nrow=nloci2)[,subset]
    if (with.id) list(cdata=t(gret),obscom=z$obscom,idsave=z$idsave[subset],wt=z$wt[subset])
    else list(cdata=t(gret),obscom=z$obscom,wt=z$wt[subset])
}
