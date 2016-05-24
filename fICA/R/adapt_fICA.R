adapt_fICA <- function(X, gs=gf, dgs=dgf, name=gnames, kj=0, inR=TRUE, eps=1e-06, maxiter=100)
{
   n <- nrow(X)
   p <- ncol(X)
   eps <- p*eps 
 
   init_est <- "k-JADE"   

   if(!(kj %in% 1:p)){
     W0 <- FOBI(X)$W  
     init_est <- "FOBI"  
   }else{
     W0 <- k_JADE(X,k=kj,eps=eps,maxiter=maxiter)$W 
     init_est <- paste(kj,"-JADE",sep="")
   }   
   Z <- tcrossprod(X,W0)
   Z <- sweep(Z,2,colMeans(Z))

   if(inR){
     if(length(name)!=length(gs)){
       name<-NULL
       for(i in 1:length(gs)){
        name[i] <- paste("g",i) 
       }
     } 
     res <- adapt_fICA.R(Z,gs=gs,dgs=dgs,name=name,kj=kj,eps=eps,maxiter=maxiter)
   }else{
     name <- gnames
     res <- .Call("adfica",Z,eps,maxiter,PACKAGE="fICA")
   }  
  
   cnam <- NULL
   for(i in 1:p){
     cnam[i] <- paste("comp",i)
   }

   V <- res$W 
   alphas <- res$alphas
   ord <- res$ord+1
   usedg <- res$usedg+1   
   used_gs <- NULL
   for(i in 1:(p-1)){
    used_gs[i] <- name[usedg[i]] 
   }
   
   if(length(ord)==(p-1)){
     ord[p] <- sum(1:p)-sum(ord)
   }else ord <- 1:p

   W <- crossprod(V,W0)
   W <- crossprod(diag(sign(rowMeans(W))),W)
   S <- tcrossprod(sweep(X,2,colMeans(X)),W)

   alphas <- matrix(alphas[,ord],ncol=p)
   dimnames(alphas) <- list(name, cnam)
   res <- list(W=W, gs=name, used_gs=used_gs, alphas=alphas, init_est=init_est, S=S)
   class(res) <- "bss"
   res
}
 
adapt_fICA.R <- function(Z, gs, dgs, name, kj, eps, maxiter)
{
   n <- nrow(Z)
   p <- ncol(Z)
   ng <- length(gs)
   usedg <- NULL 

   alphas <- compute_alphas(Z,gs,dgs,name)
   ca <- ifelse(alphas>0,alphas,Inf)
   
   ord <- NULL
   VN <- diag(p)
   V <- matrix(0,p,p)
 
   for(i in 1:(p-1)){
     mina <- which.min(ca)
     comp <- ceiling(mina/ng)
     gc <- mina-(comp-1)*ng
     vn <- VN[,comp]
     iter <- 0
     a <- 0
     if(min(ca)==Inf){ stop("no convergence")
     }else{
     while(TRUE){
       iter <- iter+1
       v <- vn
       Zv <- crossprod(t(Z),v)
       vn <- colMeans(sweep(Z,1,gs[[gc]](Zv),"*"))-mean(dgs[[gc]](Zv))*v
       vn <- vn/sqrt(sum(vn^2))
       if(sqrt(sum((vn-v)^2))>1) vn <- -vn
       if((a>0)&&(floor(iter/a)==iter/a)){
        vn <- (1-1/5)*vn+v/5
       }else vn <- (1-1/(iter+20))*vn+v/(iter+20) 

       vn <- vn-crossprod(tcrossprod(V,V),vn) 
       vn <- vn/sqrt(sum(vn^2))
    
       if(sqrt(sum((v-vn)^2))<eps || sqrt(sum((v+vn)^2))<eps){
         usedg[i] <- mina-(comp-1)*ng
         ord[i] <- comp
         ca[,comp] <- Inf 
         break
       }
       if(iter==maxiter){
        a <- a+1
        iter <- 0
        if(a>10){ 
         ca[mina] <- Inf
         mina <- which.min(ca)[1]
         comp <- ceiling(mina/ng)   
         vn <- VN[,comp]
         if(min(ca)==Inf) stop("no convergence")
        }
       }
     }
     V[,i] <- t(vn)
    }
   }
   vn <- VN[,1]   
   vn <- vn-crossprod(tcrossprod(V,V),vn)   
   vn <- vn/sqrt(sum(vn^2))
   V[,p] <- t(vn)

   res <- list(W=V, usedg=usedg-1, alphas=alphas, ord=ord-1)
   res
}

