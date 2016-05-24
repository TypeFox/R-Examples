
 depth.RPD<-function (fdataobj,fdataori=fdataobj, nproj = 50, proj=1,deriv = c(0, 1), trim = 0.25,
    dfunc2 = depth.mode, method = "fmm", draw = FALSE, ...)
{
    if (!is.fdata(fdataobj))         fdataobj = fdata(fdataobj)
    if (!is.fdata(fdataori)) fdataobj=fdata(fdataori)    
#     nas <- apply(fdataobj$data, 1, count.na)
#     if (any(nas)) {
#         fdataobj$data <- fdataobj$data[!nas, ]
#         cat("Warning: ", sum(nas), " curves with NA are not used in the calculations \n")
#     }
#     data <- fdataobj[["data"]]
 if (is.null(rownames(fdataobj$data)))  rownames(fdataobj$data)<-1:nrow(fdataobj$data)
 nms<-rownames(fdataobj$data)
 m0<-nrow(fdataobj)
  fdataobj2<- fdataobj
 fdataobj<-na.omit.fdata(fdataobj)
 fdataori<-na.omit.fdata(fdataori) 
 nas<-na.action(fdataobj)
 nullans<-!is.null(nas) 
    data2 <- fdataori[["data"]]    
    data <- fdataobj[["data"]]        
    names1 <- names2 <- names <- fdataobj[["names"]]
    names1$main <- "depth.RPD median"
    names2$main <- paste("RPD trim ", trim * 100, "%", sep = "")
    n <- nrow(data)
    m <- ncol(data)
    m2<-ncol(data2)
    n2<-nrow(data2)    
    if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
    modulo = function(z) {        sqrt(sum(z^2))    }
    if (is.null(n) || is.null(m))       stop("Input must be a matrix")
    tt = fdataobj[["argvals"]]
    rtt <- fdataobj[["rangeval"]]
    newfunc = array(NA, dim = c(n, m, length(deriv)))
    newfunc2 = array(NA, dim = c(n2, m2, length(deriv)))    
    for (ider in 1:length(deriv)) {
        if (deriv[ider] == 0) {   
          newfunc[, , ider] = data     
          newfunc2[, , ider] = data2     
          } 
        else {
            newfunc[, , ider] = fdata.deriv(fdataobj, nderiv = (ider -
                1), method = method, ...)[["data"]]
            newfunc2[, , ider] = fdata.deriv(fdataori, nderiv = (ider -
                1), method = method, ...)[["data"]]               
        }
    }
    dep = rep(0, n)
    dep2 = rep(0, n2)    
    vproject = matrix(0, nrow = n, ncol = length(deriv))
    vproject2 = matrix(0, nrow = n2, ncol = length(deriv))    
    z = matrix(rnorm(m * nproj), nrow = nproj, ncol = m)
    modu = apply(z, 1, modulo)
    z = z/modu
    if (is.fdata(proj)) {
     if (fdataobj$argvals!=proj$argvals || m!=ncol(proj)) stop("Error en proj dimension")
     z<-proj
     nproj<-nrow(z)
    }
  else {	 z<-rproc2fdata(nproj,tt,sigma=proj,norm=TRUE,...)	}
  pb = txtProgressBar(min = 0, max = nproj, style = 3)
  for (j in 1:nproj) {
        setTxtProgressBar(pb, j - 0.5)
        for (ider in 1:length(deriv)) {
            matriz = newfunc[, , ider]
            vproject[, ider] = matriz %*%z$data[j, ]
            vproject2[, ider] = newfunc2[, , ider] %*%z$data[j, ]              
        }
#        resul = dfunc2(vproject, ...)       
        par.dfunc = list()
        par.dfunc$fdataobj <- fdata(vproject)
        par.dfunc$fdataori <- fdata(vproject2)
        par.dfunc$trim <- trim              
        par.dfunc$scale<-TRUE
        resul = do.call(dfunc2, par.dfunc)
        dep = dep + resul$dep
        setTxtProgressBar(pb, j)
    }                                                               
    close(pb)
if (nullans) {
        ans1<-rep(NA,len=m0)
        ans1[-nas] <-dep
        dep<-ans1      
        }
    names(dep)<-nms       
    dep = dep/nproj
    k = which.max(dep)    
    med = data[k, ]
    lista = which(dep >= quantile(dep, probs = trim,na.rm=TRUE))
    mtrim = apply(fdataobj2$data[lista, ,drop=FALSE], 2, mean, na.rm = TRUE)
    tr <- paste("RPD.tr", trim * 100, "%", sep = "")
    med <- fdata(med, tt, rtt, names1)
    mtrim <- fdata(mtrim, tt, rtt, names2)
    rownames(med$data) <- "RPD.med"
    rownames(mtrim$data) <- tr
    if (draw) {
        ans <- dep
        ind1 <- !is.nan(ans)
        ans[is.nan(ans)] = NA
        cgray = 1 - (ans - min(ans, na.rm = TRUE))/(max(ans,
            na.rm = TRUE) - min(ans, na.rm = TRUE))
        plot(fdataobj[ind1, ], col = gray(cgray[ind1]), main = "RPD Depth")
        lines(mtrim, lwd = 2, col = "yellow")
        lines(med, col = "red", lwd = 2)
        legend("topleft", legend = c(tr, "Median"), lwd = 2,box.col=0,
            col = c("yellow", "red"))
    }
    return(invisible(list(median = med, lmed = k, mtrim = mtrim,
        ltrim = lista, dep = dep,proj = z)))
}  
