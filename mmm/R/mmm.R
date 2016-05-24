mmm <-
function (formula, id, data=NULL, correlation = NULL, initEstim = NULL, tol = 0.001, maxiter = 25, 
                  family = "gaussian", corStruct = "independence", Mv = 1, silent = TRUE, 
                  scale.fix = FALSE, scale.value = 1) 
{
    mf<-model.frame(formula=formula,data=data)
    x <- as.matrix(model.matrix(attr(mf, "terms"), data=mf)[,-1])
    if (ncol(x)==1) colnames(x)<-colnames(model.matrix(attr(mf, "terms"), data=mf))[2:(ncol(x)+1)]
    y<-model.response(mf)
    data<-as.matrix(cbind(id,y,x))
    nresp<-ncol(y)
    id2 <- unique(data[, 1])
    resp <- NULL
    for (i in id2) {
        resp1 <- data[data[, 1] == i, 2:(1 + nresp)]
        resp2 <- NULL
        for (j in 1:nresp) {
            resp2 <- c(resp2, resp1[, j])
        }
        resp <- c(resp, resp2)
    }
    resp <- matrix(resp)
    ones <- rep(1, dim(data)[1])
    covmat <- cbind(data[, 1], ones, data[, (1 + nresp + 1):dim(data)[2]])
    cov3 <- NULL
    for (k in id2) {
        cov1 <- covmat[covmat[, 1] == k, 2:dim(covmat)[2]]
        cov2 <- kronecker(diag(1, nresp), as.matrix(cov1))
        cov3 <- rbind(cov3, cov2)
    }
    cov3<-data.frame(cov3)
    if (length(colnames(y))==0){
    vn1<-gsub("^.*\\$", "", unlist(strsplit(colnames(mf)[1], "[,)[:blank:]]")))
    vn1<-vn1[!((1:length(vn1))%%2==0)]
    }else{
    vn1<-colnames(y)
    }   
    vn2<-c("Intercept",colnames(x))
    vn3<-rep(vn1,each=length(vn2))
    vnames<-paste(vn3,vn2,sep=".")
    colnames(cov3)<-vnames
    id5 <- NULL
    for (t in id2) {
        id3 <- data[data[, 1] == t, 1]
        id4 <- rep(id3, nresp)
        id5 <- c(id5, id4)
    }
   colnames(cov3)<-gsub('$','.',colnames(cov3),fixed=T) 
   formula2<-as.formula(paste("resp ~ -1+", paste(colnames(cov3), collapse= "+")))
    #library(gee)
    fit <- gee(formula=formula2, id = id5, data=cov3, R=correlation, b=initEstim, tol=tol, maxiter=maxiter, 
                   family = family, corstr = corStruct, Mv = Mv, silent = silent,
                   scale.fix=scale.fix, scale.value=scale.value)
    fit$title<-"Multivariate Marginal Models"
    fit$version<-"Version 1.4 (01/2014)"
    fit
}
