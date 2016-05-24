reloaded_fICA <- function(X, g="tanh", dg=NULL, G=NULL, kj=0, method="alpha", inR=TRUE, eps=1e-06, maxiter=100)
{
   n <- nrow(X)
   p <- ncol(X)
   eps <- p*eps
   name <- c("pow3","tanh","gaus")
   init_est <- "k-JADE"   
   
   method <- match.arg(method,c("alpha","G"))

   if(!(kj %in% 1:p)){
     W0 <- FOBI(X)$W  
     init_est <- "FOBI"  
   }else{
     W0 <- k_JADE(X,k=kj,eps=eps,maxiter=maxiter)$W 
     init_est <- paste(kj,"-JADE",sep="")
   }   
   Z <- tcrossprod(X,W0)
   Z <- sweep(Z,2,colMeans(Z))
  
   res<-switch(method,
        "alpha"={
           if(!inR){
            gi <- which(name==g[1])
            res <- .Call("relfica",Z,gi,eps,maxiter,PACKAGE="fICA")
           }else{ 
            if(!(is.function(g)&&is.function(dg))){
             gi <- which(name==g[1])
             g1 <- gf[[gi]]
             dg1 <- dgf[[gi]]    
            }else{
             g1 <- g
             dg1 <- dg
            }
           res <- reloaded_fICA.R(Z,g=g1,dg=dg1,eps=eps,maxiter=maxiter)
          }
         }
        ,
        "G"={ 
             res <- G_fICA.R(Z,g=g,dg=dg,G=G,eps=eps,maxiter=maxiter)       
            }
    )
   
   V <- res$W 
   alphas <- res$alphas
   ord <- res$ord+1

   if(length(ord)==(p-1)){
     ord[p] <- sum(1:p)-sum(ord)
   }else ord <- 1:p


   W <- crossprod(V,W0)
   W <- crossprod(diag(sign(rowMeans(W))),W)
   S <- tcrossprod(sweep(X,2,colMeans(X)),W)

   res <- list(W=W, g=g, method=method, alphas=alphas[ord], init_est=init_est, S=S)
   class(res) <- "bss"
   res
}


reloaded_fICA.R <- function(Z, g=NULL, dg=NULL, eps=1e-04, maxiter=100)
{
   n <- nrow(Z)
   p <- ncol(Z)
   eps <- p*eps 

   alphas <- compute_alphas(Z,c(g),c(dg),"g")
       
   ca <- ifelse(alphas>0,alphas,Inf)
   
   ord <- NULL
   VN <- diag(p)
   V <- matrix(0,p,p)
 
   for(i in 1:(p-1)){
     comp <- which.min(ca)
     vn <- VN[,comp]
     iter <- 0
     a <- 0
     while(TRUE){
       iter <- iter+1
       v <- vn
       Zv <- crossprod(t(Z),v)
       vn <- colMeans(sweep(Z,1,g(Zv),"*"))-mean(dg(Zv))*v
       vn <- vn/sqrt(sum(vn^2))
       if(sqrt(sum((vn-v)^2))>1) vn <- -vn
       if((a>0)&&(floor(iter/a)==iter/a)){ 
        vn <- (1-1/5)*vn+v/5
       }else vn <- (1-1/(iter+20))*vn+v/(iter+20) 

       vn <- vn-crossprod(tcrossprod(V,V),vn) 
       vn <- vn/sqrt(sum(vn^2))
    
       if(sqrt(sum((v-vn)^2))<eps || sqrt(sum((v+vn)^2))<eps){
         ord[i] <- comp
         ca[,comp] <- Inf 
         break
       }
       if(iter==maxiter){
        a <- a+1
        iter <- 0
        if(a>10) stop("no convergence")
       }
     }
     V[,i] <- t(vn)
   }
   vn <- VN[,1]   
   vn <- vn-crossprod(tcrossprod(V,V),vn)   
   vn <- vn/sqrt(sum(vn^2))
   V[,p] <- t(vn)

   res <- list(W=V, alphas=alphas, ord=ord-1)
   class(res) <- "bss"
   res
}


G_fICA.R <- function(Z,g,dg,G,eps,maxiter)
{
   n <- nrow(Z)
   p <- ncol(Z)
   name <- c("pow3","tanh","gaus")
   EGn <- c(0.75,0.3745672,-0.7071068)
   if((is.function(g)&&is.function(dg)&&is.function(G))){
     Gn <- integrate(Vectorize(function(x){G(x)*dnorm(x)}),-10,10)$value  
   }else{
     gi <- which(name==g[1])
     g <- gf[[gi]]
     dg <- dgf[[gi]]                
     G <- Gf[[gi]]
     Gn <- EGn[gi] 
   }
 
   EG <- NULL
   for(i in 1:p){
     EG[i] <- abs(mean(G(Z[,i]))-Gn)
   }  
 
   ord <- order(EG,decreasing=TRUE)

   VN <- diag(p)
   V <- matrix(0,p,p)

   for(i in 1:p){ 
     vn <- VN[,ord[i]]
     a <- 0
     iter <- 0 
     while(TRUE){
       iter <- iter+1 
       v <- vn
       Zv <- crossprod(t(Z),v)
       vn <- colMeans(sweep(Z,1,g(Zv),"*"))-mean(dg(Zv))*v
       vn <- vn/sqrt(sum(vn^2))
       if(sqrt(sum((vn-v)^2))>1) vn <- -vn
       if((a>0)&&(floor(iter/a)==iter/a)){
        vn <- (1-1/5)*vn+v/5
       }else vn <- (1-1/(iter+20))*vn+v/(iter+20) 
       vn <- vn-crossprod(tcrossprod(V,V),vn)  
       vn <- vn/sqrt(sum(vn^2))
       if(sqrt(sum((v-vn)^2))<eps || sqrt(sum((v+vn)^2))<eps) break
       if(iter==maxiter){
        a <- a+1
        iter <- 0
        if(a>10) stop("maxiter reached without convergence") 
       } 
      }
    V[,i] <- t(vn)
  }   
  list(W=V, alphas=EG, ord=ord[-p]-1) 
}

