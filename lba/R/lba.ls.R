lba.ls <- function(obj        ,             
                   A          ,      
                   B          ,      
                   K          ,
                   row.weights, 
                   col.weights, 
                   tolA       , 
                   tolB       ,
                   itmax.unide,
                   itmax.ide  ,
                   trace.lba  ,
                   what,
                   ...)
{

 I  <- nrow(obj)     # row numbers of data matrix
 J  <- ncol(obj)       # column numbers of data matrix

 P  <- obj/rowSums(obj) # observed components I x J

 # Generating the initial mixing parameters if they aren't informed 

 if(is.null(A)){
  A <- rdirich(I, 
               runif(K))
 }
 # Generating the initial latent components if they are not informed
 if(is.null(B)){
  B <- t(rdirich(K, 
                 runif(J)))
 }
 # Generating the identity matrix of row weights if they aren't informed
 if(is.null(row.weights)){
  vI <- rep(1,I)
  V  <- vI * diag(I)
 } else {
  vI <- row.weights
  V <- vI * diag(I)
 }

 # Generating the identity matrix of column weights if they aren't informed 
 if(is.null(col.weights)){
  wi <- rep(1,J)
  W  <- wi * diag(J)
 } else {
  wi <- col.weights
  W <- wi * diag(J)
 }

 P_ast_trans <- W %*% t(P) %*% V
 P_ast       <- t(P_ast_trans)

 iter_unide <-  0 # interactions counter

 vec <- function(X){
  matrix(X)
 }

 ls_func <- function(V,P,pij,W){
  val_function <- sum((V %*% (P - pij) %*% W)^2)   
 }

 repeat{
  Aa <- A  # start values
  Ba <- B  # start values

  pij <- A %*% t(B)

  iter_unide <- iter_unide + 1

  if(trace.lba){
   cat('Unidentified iteration: ',iter_unide,'\n')
   cat('fval = ', signif(ls_func(V,P,pij,W),4), '\n')
  }
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Step 1. Estimanting the betas (latent components) - A constant
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
  C1 <- diag(K) %x% t(rep(1,J))  # is a K x JK matrix 
  d1 <- matrix(1,
               nrow = K, 
               ncol = 1) # is a K x 1 matrix 

  Qb = (V %*% A) %x% W # is a IJ x JK matrix

  rb = vec(P_ast_trans) # is a IJ x 1 matrix

  QQinvb <- ginv(t(Qb) %*% Qb)

  x0b = QQinvb %*% t(Qb) %*% rb 

  lambdab <- ginv(C1 %*% QQinvb %*% t(C1)) %*% (d1 - C1 %*% x0b)

  xB  <- vec(x0b + QQinvb %*% t(C1) %*% lambdab)

  if(any(as.vector(xB) < 0)){
   C2 <- NULL
   d2 <- NULL
   while(any(as.vector(xB) < -1e-10)){
    C2 <- C2
    d2 <- d2

    Bw <- which(xB < 0, 
                arr.ind=TRUE) 

    b0  <- matrix(0, 
                  nrow=nrow(xB), 
                  ncol=dim(Bw)[1])

    for(i in 1:dim(Bw)[1]) b0[sapply(Bw[i,1],function(x)x),i] <- 1

    C2 <- rbind(C2,t(b0))

    d2 <- matrix(c(d2,rep(0,dim(Bw)[1])))

    C <- rbind(C1,C2)

    d <- rbind(d1,d2)

    lambdab <- ginv(C %*% QQinvb %*% t(C)) %*% (d - C %*% x0b)

    xB = matrix(x0b + QQinvb %*% t(C) %*% lambdab)

    if(any(dim(lambdab)[1] > J*K)) break
   }
   xB[xB < 0] <- 0
  }
  B <- matrix(xB, ncol = K)

  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Step 2. Estimanting the alfas (mixing parameters) - B constant
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  C1 <- matrix(1,nrow = 1, ncol = K) # uma matrix 1 x K 
  d1 <- matrix(1,1,1) 

  Qa <- lapply(vI, function(x) x * W %*% B)

  tQa <- lapply(Qa,
                t)

  ra <- lapply(apply(P_ast, 
                     1, 
                     function(x) list(matrix(x,ncol=1))),'[[',1)

  QQinva <- lapply(mapply('%*%', 
                          tQa,
                          Qa, 
                          SIMPLIFY=F),
                   function(x) ginv(x))

  x0a <- mapply('%*%',
                mapply('%*%',
                       QQinva,
                       tQa,
                       SIMPLIFY=F),
                ra,
                SIMPLIFY=F)

  auxlambda1 <-  lapply(QQinva,
                        function(x) ginv(C1 %*% x %*%t(C1)))

  auxlambda2 <-  lapply(x0a,
                        function(x) d1 - C1 %*% x)

  lambdaa <- mapply('%*%',
                    auxlambda1,
                    auxlambda2,
                    SIMPLIFY=F)

  xA <- mapply('+',
               x0a,
               mapply('%*%',                         
                      lapply(QQinva,
                             function(x) x %*% t(C1)),
                      lambdaa,
                      SIMPLIFY=F),
               SIMPLIFY=F)

  if(any(unlist(xA) < 0)){
   C2 <- NULL
   d2 <- NULL

   while(any(unlist(xA) < -1e-10)){
    C2 <- C2
    d2 <- d2

    Aw <- lapply(xA, function(x) which(x <0, arr.ind  = TRUE))

    a0 <- lapply(Aw, function(x) matrix(0,nrow=dim(x)[1],ncol=K))

    for(i in 1:I){
     for(j in 0:dim(Aw[[i]])[1]){
      a0[[i]][j , Aw[[i]][j,1]]  <- 1
     }
    }

    a0 <- lapply(a0, function(x) if(dim(x)[1] == 0) NULL else x)

    if(is.null(C2)){
     C2 <-  a0
    } else {
     C2 <- mapply('rbind',
                  C2,
                  a0,
                  SIMPLIFY=F)
    }

    C <- list()

    for(i in 1:length(C2)){

     C[[i]] <- rbind(C1,C2[[i]])

    }

    a1 <- lapply(Aw,
                 function(x) matrix(0,nrow = dim(x)[1],ncol=1))

    if(is.null(d2)){
     d2 <- a1 
    } else {
     d2 <- mapply('rbind',
                  d2,
                  a1,
                  SIMPLIFY=F)
    }

    d <- mapply('rbind',
                d1,
                d2,
                SIMPLIFY=F)

    CT <-   lapply(C,
                   function(x)t(x))

    QTQ <- mapply('%*%',  
                  lapply(Qa,
                         function(x)t(x)),
                  Qa,
                  SIMPLIFY=F)# transposta de Q x Q

    invQTQ <-  lapply(QTQ,
                      function(x)ginv(x))# inversa de QTQ

    CinvQTQCT <-  mapply('%*%',             
                         mapply('%*%',
                                C,
                                invQTQ,
                                SIMPLIFY=F),                       
                         CT,
                         SIMPLIFY=F)

    invCinvQTQCT <-   lapply(CinvQTQCT,
                             function(x)ginv(x))

    d_Cx0a <-  mapply('-',
                      d,
                      mapply('%*%',
                             C,
                             x0a,
                             SIMPLIFY=F),
                      SIMPLIFY=F)


    lambdaa <- mapply('%*%',
                      invCinvQTQCT,                                     
                      d_Cx0a,
                      SIMPLIFY=F)

    CTlambdaa <- mapply('%*%',
                        CT,
                        lambdaa,
                        SIMPLIFY=F)

    invQTQCTlambdaa <- mapply('%*%',
                              invQTQ,
                              CTlambdaa,
                              SIMPLIFY=F)

    xA <- mapply('+',
                 x0a,
                 invQTQCTlambdaa,
                 SIMPLIFY=F)

    if(any(unlist(lapply(lambdaa, function(x) dim(x)[1] > I*K))) == TRUE) break
   }
   for(i in 1:I) xA[[i]][xA[[i]] < 0]  <- 0
  }

  ifelse(K==1,
         A <- matrix(unlist(xA),
                     ncol = 1),

         A <- t(sapply(xA, 
                       rbind, 
                       simplify=T))) 

  at <- max(abs(A - Aa))
  bt <- max(abs(B - Ba))

  if(iter_unide > itmax.unide){
   warning("maximum number of the unidentified iteractions exceeded" )
   break 
  }

  if(at < tolA & bt < tolB) break

 }

 A <- abs(A)

 B <- abs(B)

 pimais <- rowSums(obj)/sum(obj)

 if(K==1){

  Aoi <- A

  Boi <- B

  iter_ide  <- 'not applicable'

 } else {

  AB_outerAB_inner <- identif(A,
                               B,
                               P,
                               K,
                               trace.lba,
                               itmax.ide,
                               what)

  Aoi <- AB_outerAB_inner[[1]]

  Boi <- AB_outerAB_inner[[2]]

  iter_ide <- AB_outerAB_inner[[3]]  

 }

 pij <- Aoi %*% t(Boi) # expected budget

 rownames(pij) <- rownames(A) <- rownames(P)
 colnames(pij) <- rownames(B) <- colnames(P)

 residual <- P - pij

 pk <- pimais %*% Aoi # budget proportions
 colnames(pk) <- colnames(A) <- colnames(B) <- paste('LB',1:K,sep='') 
 val_func <- ls_func(V,P,pij,W) 

 rescB <- rescaleB(obj,
                   Aoi,
                   Boi)

 colnames(rescB) <- colnames(B)
 rownames(rescB) <- rownames(B)

 res <- list(P,       
             pij,
             residual,
             A, 
             B,
             Aoi,
             Boi,
             rescB,
             pk,
             val_func,
             iter_unide,
             iter_ide)

 names(res) <- c('P',
                 'pij',
                 'residual',
                 'A',
                 'B',
                 'Aoi',
                 'Boi',
                 'rescB',
                 'pk',
                 'val_func',
                 'iter_unide',
                 'iter_ide')

# 
 #  names(res) <- c('Composition data matrix',
 #                  'Expected budget',
 #                  'Residual matrix',
 #                  'Unidentified mixing parameters',
 #                  'Unidentified latent budgets',
 #                  ifelse(what=='outer',
 #                         'Outer extreme mixing parameters',
 #                         'Inner extreme mixing parameters'),
 #                  ifelse(what=='outer',
 #                         'Outer extreme latent budgets',
 #                         'Inner extreme latent budgets'),
 #                  'Rescaled latent budgets',
 #                  'Budget proportions',
 #                  'Value of the ls function',
 #                  'Number of unidentified iteractions',
 #                  'Number of identified iteractions')
 # 

 class(res) <- 'lba.ls'

 invisible(res)
} 
