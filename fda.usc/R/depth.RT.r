 depth.RT  <-function (fdataobj,fdataori=fdataobj, trim = 0.25, nproj = 10, proj = 1, xeps = 1e-07, 
    draw = FALSE, ...) 
{
    if (!is.fdata(fdataobj)) fdataobj = fdata(fdataobj)
    if (!is.fdata(fdataori)) fdataobj=fdata(fdataori)    
#    nas <- apply(fdataobj$data, 1, count.na)
#    if (any(nas)) {
#        fdataobj$data <- fdataobj$data[!nas, ]
#        cat("Warning: ", sum(nas), " curves with NA are not used in the calculations \n")
#    }
 if (is.null(rownames(fdataobj$data)))  rownames(fdataobj$data)<-1:nrow(fdataobj$data)
 nms<-rownames(fdataobj$data)
 m0<-nrow(fdataobj)
 fdataobj2<-fdataobj
 fdataobj<-na.omit.fdata(fdataobj)
 nas<-na.action(fdataobj)
 nullans<-!is.null(nas) 
 if (nullans)  fdataori<-fdataori[-nas]  
    data  <- fdataobj[["data"]]
    data2 <- fdataori[["data"]]    
    n <- nrow(data)
    m <- ncol(data)
    m2<-ncol(data2)
    n2<-nrow(data2)
    if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS OF fdataobj")
    if (is.null(n2) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS fdataori")
    t = fdataobj[["argvals"]]
    rtt <- fdataobj[["rangeval"]]
    names1 <- names2 <- names <- fdataobj[["names"]]
    names1$main <- "depth.RT median"
    names2$main <- paste("RT trim", trim * 100, "%", sep = "")
   
    if (is.fdata(proj)) {
        nproj <- nrow(proj)
        if (fdataobj$argvals != proj$argvals || ncol(fdataobj) != 
            ncol(proj)) 
            stop("Error en proj dimension")
        z <- proj
    }
    else {
        z <- rproc2fdata(nproj, t, sigma = proj, norm = TRUE,...)
    }
    Fn <- list()
    Prod=data%*%t(z$data)
    Prod2=data2%*%t(z$data)    
#    dep = array(NA, dim = c(nproj, n))
    dep2 = array(NA, dim = c(nproj, n2))    
    
 #   dep2=rep(0.0,n) 
    for (j in 1:nproj) {
        Fn[[j]] = ecdf(Prod2[,j])  
#        dep[j, ] = pmin(Fn[[j]](Prod[,j]) ,(1 - Fn[[j]](Prod[,j]-xeps)))       
        dep2[j, ] = pmin(Fn[[j]](Prod[,j]) ,(1 - Fn[[j]](Prod[,j]-xeps)))               
#        dep2 =dep2+ pmin(Fn[[j]](Prod[,j]) ,(1 - Fn[[j]](Prod[,j]))) #idem que dep
    }
#    print(dep)
#    dep = apply(dep, 2, min)
#    dep = colMeans(dep)
    dep2 = apply(dep2,2,mean)
    if (nullans) {
        ans1<-rep(NA,len=m0)
        ans1[-nas] <-dep2
        dep2<-ans1      
        }
     names(dep2)<-nms      
    k = which.max(dep2)
    med = fdataobj2$data[k, ]
    nl = length(trim)
    mtrim = matrix(NA, nrow = nl, ncol = m)
    for (j in 1:length(trim)) {
        lista = which(dep2 >= quantile(dep2, probs = trim[j], na.rm = TRUE))
        if (length(lista)==1) {
          mtrim[j, ]<-fdataobj2$data[lista,]
          if (draw) {draw=FALSE;warning("The plot is not shown")}
        }
        else mtrim[j, ]=apply(fdataobj2$data[lista,],2,mean,na.rm=TRUE)       
    }
    tr <- paste("RT.tr", trim * 100, "%", sep = "")
    med <- fdata(med, t, rtt, names1)
    mtrim <- fdata(mtrim, t, rtt, names2)
    rownames(med$data) <- "RT.med"
    rownames(mtrim$data) <- tr
    if (draw) {
        ans <- dep2
        ind1 <- !is.na(ans)
        cgray = 1 - (ans - min(ans, na.rm = TRUE))/(max(ans, 
            na.rm = TRUE) - min(ans, na.rm = TRUE))
        plot(fdataobj2[ind1, ], col = gray(cgray[ind1]), main = "RT Depth")
        lines(mtrim, lwd = 2, col = "yellow")
        lines(med, col = "red", lwd = 2)
        legend("topleft", legend = c(tr, "Median"), lwd = 2,box.col=0, 
            col = c("yellow", "red"))
    }
    return(invisible(list(median = med, lmed = k, mtrim = mtrim, 
        ltrim = lista, dep = dep2, proj = z, Fn = Fn)))
}


