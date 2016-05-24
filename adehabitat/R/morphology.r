"morphology" <- function(x, operation = c("erode", "dilate"), nt=5)
{
    ## Verifications
    op<-match.arg(operation)
    if (nt<1)
      stop("nt should be > 0")
    if (op=="erode")
      ope<-0
    if (op=="dilate")
      ope<-1
    if (!inherits(x,"asc"))
      stop("should be of class asc")

    ## Bases
    nc<-ncol(x)
    nr<-nrow(x)
    tmpc<-rep(NA,nc)
    u<-rbind(tmpc,x,tmpc)
    tmpl<-rep(NA,nr+2)
    u<-cbind(tmpl,u,tmpl)
    o<-as.vector(t(u))
    o[!is.na(o)]<-1
    o[is.na(o)]<-0

    ## External call to the C function "erodil"
    res<-.C("erodil", as.double(o), as.integer(nr+2), as.integer(nc+2),
            as.integer(nt), as.integer(ope), PACKAGE="adehabitat")

    ## Output
    res[[1]][res[[1]]==0]<-NA
    gr<-matrix(res[[1]], nrow=(nr+2), byrow=TRUE)
    gr <- gr[-c(1, nrow(gr)), -c(1, ncol(gr))]
    if (all(is.na(gr)))
        stop("all the image has been erased\n Please consider a lower value for nt")
    gr<-getascattr(x,gr)
    return(gr)
}

