EPRM_red <-
function(daten, score_par=NULL){
  
  call <- match.call()
    
  if(is.data.frame(daten)) {daten <- as.matrix(daten)}
  
  kateg.zahl <- length(table(daten))
  item.zahl <- ncol(daten)
  
  if(kateg.zahl <= 2 ){stop("there are only 2 categories for testing reduction!")}
  # margin vector groups of persons
  
  row.table <- apply(daten+1,1,function(x) sprintf("%04d",(tabulate(x,nbins=kateg.zahl))))
  
  pat   <- apply(row.table,2, function(n) paste0(n, collapse=""))
  patt  <- table(pat)
  
  #first term (last category left out because these item parameters are 0)
  
  col.table <- apply(daten+1, 2, function(s) tabulate(s,nbins=kateg.zahl))
  
  #starting values
  if(!is.null(score_par)){
    startval <- rep(0,(item.zahl-1))
  } else {
    startval <- rep(0,(item.zahl-1)+(kateg.zahl-2)) 
  }
    
  #pattern
  
  patmat <- t(xsimplex(kateg.zahl,item.zahl))
  patmat.o <- patmat[order(patmat[,kateg.zahl], decreasing=T),]
  
  patt.c <- apply(patmat.o,1,function(p) paste0(sprintf("%04d",p),collapse=""))
  
  cL <- function(para=startval,kateg.zahl=kateg.zahl, item.zahl=item.zahl,col.table=col.table, patmat.o=patmat.o, patt.c=patt.c, patt=patt, score_par=score_par){
    
    if(!is.null(score_par)){
      para_ext <- outer(c(1,score_par,0),c(para,0))
    } else {
      para_ext <- outer(c(1,para[1:(kateg.zahl-2)],0),c(para[(kateg.zahl-1):(length(para))],0))
    }
    
    fir <- sum(col.table*para_ext)
    
    cf.g <- combfunc(kateg.zahl, item.zahl, eps.mat=t(exp(para_ext)), patmat.o)
    
    ord <- sapply(names(patt), function(hpat) which(patt.c %in% hpat))
    cfg.o <- log(cf.g$gammat[ord])
    sec <- sum(patt*cfg.o,na.rm=T)
    
    fir-sec
  }
    
  
  der1 <- function(para=startval, col.table=col.table, kateg.zahl=kateg.zahl, item.zahl=item.zahl, patmat.o=patmat.o, patt=patt, patt.c=patt.c, score_par=score_par){
    
    if(!is.null(score_par)){
      
      para_ext <- outer(c(1,score_par,0),c(para,0))
      
    } else {
      
      para_ext <- outer(c(1,para[1:(kateg.zahl-2)],0),c(para[(kateg.zahl-1):(length(para))],0))
    }
        
    eps.f <- list(exp(para_ext))
    
    for (e in seq_len(item.zahl-1)){
      eps.f[[e+1]] <- cbind(eps.f[[e]][,2:item.zahl],eps.f[[e]][,1])
    }
    cf.all <- lapply(eps.f, function(gg) {combfunc(kateg.zahl, item.zahl, eps.mat=t(gg), patmat.o)})
    
    #match pattern in patt.c
    
    ord2 <- sapply(names(patt), function(hpat) which(patt.c %in% hpat))
    cf.o <- lapply(cf.all, function(lo) lo$gam.quot[ord2,])
    
    cf.oNR <- lapply(cf.o, function(l3) {colSums(as.vector(patt)*l3, na.rm=TRUE)})
    
    cf.oM <- do.call(rbind,cf.oNR)
    
    cf.oMp <- t(exp(para_ext)[-kateg.zahl,]*t(cf.oM))
    
    if(!is.null(score_par)){
      firD  <- NULL
      secD  <- NULL
      
      firD2 <- colSums(c(1,score_par)*col.table[-kateg.zahl,-item.zahl])
      secD2 <- colSums(t(cf.oMp[-item.zahl,])*c(1,score_par))
      
    } else {
      firD <- colSums(as.matrix(c(para[(kateg.zahl-1):length(para)],0)*t(col.table[-c(kateg.zahl),])))[-1]
      
      secD <- colSums(as.matrix((cf.oMp*c(para[(kateg.zahl-1):length(para)],0))[,-1]))

      firD2 <- colSums(c(1,para[1:(kateg.zahl-2)])*col.table[-kateg.zahl,-item.zahl])
      
      secD2 <- colSums(t(cf.oMp[-item.zahl,])*(c(1,para[1:(kateg.zahl-2)])))
    }
    
    c(firD,firD2) - c(secD,secD2)
  }
                                                                                                                                                                                               
  res <- optim(startval, cL, gr=der1,kateg.zahl=kateg.zahl, item.zahl=item.zahl,col.table=col.table, patmat.o=patmat.o,patt.c=patt.c, patt=patt, score_par=score_par, method="BFGS", control=list(maxit=500, fnscale=-1), hessian=TRUE)   
  
  estpar_se <- sqrt(diag(solve(res$hessian*(-1))))
    
  if(!is.null(score_par)){
    item_par  <- c(res$par,0)
    score_par <- c(1, score_par,0)
  } else {
    item_par  <- c(res$par[(kateg.zahl-1):length(res$par)],0)
    score_par <- c(1,res$par[1:(kateg.zahl-2)],0)
  } 
    
  itmat_v <- outer(score_par,item_par)
  itmat <- t(apply(itmat_v,1, function(ro) scale(ro, TRUE, FALSE)))
 
  if(!is.null(colnames(daten))){
    colnames(itmat) <- paste("beta", colnames(daten))
  } else {
    colnames(itmat) <- paste("beta item", 1:ncol(itmat))
  }
    rownames(itmat) <- paste("cat", 1:nrow(itmat))
  
  
  res_all <- list(logLikelihood=res$value, estpar=res$par, estpar_se=estpar_se, score_par=score_par, item_par=item_par,itempar=itmat, data=daten, hessian=res$hessian, convergence=res$convergence, fun_calls=res$counts, call=call)
  class(res_all) <- "EPRM_red"
  res_all
  
}
