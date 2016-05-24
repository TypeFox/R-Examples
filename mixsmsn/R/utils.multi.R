######################################################
#     algoritmo para graficar os contornos

mix.contour <- function(y, model, slice = 100, ncontour = 10, x.min=1, x.max=1, y.min=1,y.max=1, ...){
   # Essa funcao so serve para graficar os contornos de misturas finitas de SMSN BIVARIADO!!!
   # dat: o cosliceunto de dat a ser plotado no R^2
   # model: deve ser um objeto resultante da funcao EMmulti.MIXSNI
   # slice e ncountor sao parametros passados para a funcao countor
   # ?contour para detalhes
   dat <- y
   y <- NULL
   n <- nrow(dat)
   p <- ncol(dat)

   if(p != 2) stop("The mix.contour function is only appropriate for the bivariate analysis.\n")
   if((class(model) != "t") && (class(model) != "Skew.t") && (class(model) != "Skew.cn") && (class(model) != "Skew.slash") && (class(model) != "Skew.normal") && (class(model) != "Normal")) stop(paste("Family",class(model),"not recognized.",sep=" "))
   if(length(model$group) == 0) stop("The groups were not save in the model.\n")
   if ((x.min < 0) || (x.max < 0) || (y.min < 0) || (y.max < 0)) stop("All limits must be non negative.\n")
   g <- length(model$pii)
   if (class(model) == "Normal"){
      mixed.Normal <- function(x, y, pii, mu, Sigma) {
        dens <- 0
        for (j in 1:g) dens <- dens + pii[j]*dmvnorm(cbind(x, y), mu[[j]], Sigma[[j]])
        dens
      }
      #x <- seq(min(dat[,1])-1,max(dat[,1])+1, length = slice)    #faithfull
      #y <- seq(min(dat[,2])-15,max(dat[,2])+10, length = slice)  #faithfull
      #x <- seq(min(dat[,1])-1,max(dat[,1])+1, length = slice)     # Swiss
      #y <- seq(min(dat[,2])-1,max(dat[,2])+1, length = slice)     # Swiss
      x <- seq(min(dat[,1])-x.min,max(dat[,1])+x.max, length = slice)    
      y <- seq(min(dat[,2])-y.min,max(dat[,2])+y.max, length = slice)  
      f <- matrix(0,slice,slice)
      for (i in 1:slice) for (j in 1:slice) f[i,j] <- mixed.Normal(x[i], y[j], model$pii, model$mu, model$Sigma)
      contour(x, y, f, , nlevels = ncontour, ...) 
      #contour(x, y, f, , levels = c(0.3,0.15, 0.075,0.0375,0.01875)) # Swiss
      points(dat[,1], dat[,2], col = (model$group+1))
   }

   if (class(model) == "Skew.normal"){
      mixed.SN <- function(x, y, pii, mu, Sigma, lambda) {
        dens <- 0
        for (j in 1:g) dens <- dens + pii[j]*2*dmvnorm(cbind(x, y), mu[[j]], Sigma[[j]])*pnorm(t(lambda[[j]])%*%solve(matrix.sqrt(Sigma[[j]]))%*%t(cbind(x, y) - mu[[j]]))
        dens
      }
      #x <- seq(min(dat[,1])-1,max(dat[,1])+1, length = slice)
      #y <- seq(min(dat[,2])-15,max(dat[,2])+10, length = slice)
      #x <- seq(min(dat[,1])-1,max(dat[,1])+1, length = slice)     # Swiss
      #y <- seq(min(dat[,2])-1,max(dat[,2])+1, length = slice)     # Swiss
      x <- seq(min(dat[,1])-x.min,max(dat[,1])+x.max, length = slice)    
      y <- seq(min(dat[,2])-y.min,max(dat[,2])+y.max, length = slice)  
      f <- matrix(0,slice,slice)
      for (i in 1:slice) for (j in 1:slice) f[i,j] <-  mixed.SN(x[i], y[j], model$pii, model$mu, model$Sigma, model$shape)
      #contour(x, y, f, nlevels = 5, levels = c(0.1, 0.015, 0.005, 0.0009, 0.00015))
      contour(x, y, f, , nlevels = ncontour, ...)
      #contour(x, y, f, , levels = c(0.3,0.15, 0.075,0.0375,0.01875)) # Swiss
      points(dat[,1], dat[,2], col = (model$group+1))
   }

   if (class(model) == "t"){
      mixed.ST <- function(x, y, pii, mu, Sigma, lambda, nu) {
#        n <- nrow(dat)
#        p <- ncol(dat)
        dens <- 0
        for (j in 1:g) {
          lambda[[j]] <- rep(0,length(lambda[[j]]))
          denst <- (gamma((p+nu)/2)/(gamma(nu/2)*pi^(p/2)))*nu^(-p/2)*det(Sigma[[j]])^(-1/2)*(1 + mahalanobis(cbind(x,y), mu[[j]], Sigma[[j]])/nu)^(-(p+nu)/2)
          dens <- dens + pii[j] * 2*(denst)*pt(sqrt((p + nu)/(mahalanobis(cbind(x,y), mu[[j]], Sigma[[j]]) + nu))*t(lambda[[j]])%*%solve(matrix.sqrt(Sigma[[j]]))%*%(c(x,y) - mu[[j]]), df = nu + p)
        }
        dens
      }
      #x <- seq(min(dat[,1])-1,max(dat[,1])+1, length = slice)
      #y <- seq(min(dat[,2])-15,max(dat[,2])+10, length = slice)
      #x <- seq(min(dat[,1])-1,max(dat[,1])+1, length = slice)     # Swiss
      #y <- seq(min(dat[,2])-1,max(dat[,2])+1, length = slice)     # Swiss
      x <- seq(min(dat[,1])-x.min,max(dat[,1])+x.max, length = slice)    
      y <- seq(min(dat[,2])-y.min,max(dat[,2])+y.max, length = slice)  
      f <- matrix(0,slice,slice)
      for (i in 1:slice) for (j in 1:slice) f[i,j] <-  mixed.ST(x[i], y[j], model$pii, model$mu, model$Sigma, model$shape, model$nu)
      #contour(x, y, f, nlevels = 5, levels = c(0.1, 0.015, 0.005, 0.0009, 0.00015))
      contour(x, y, f, , nlevels = ncontour, ...)
      #contour(x, y, f, , levels = c(0.3,0.15, 0.075,0.0375,0.01875)) # Swiss
      points(dat[,1], dat[,2], col = (model$group+1))
   }

   if (class(model) == "Skew.t"){
      mixed.ST <- function(x, y, pii, mu, Sigma, lambda, nu) {
#        n <- nrow(dat)
#        p <- ncol(dat)
        dens <- 0
        for (j in 1:g) {
          denst <- (gamma((p+nu)/2)/(gamma(nu/2)*pi^(p/2)))*nu^(-p/2)*det(Sigma[[j]])^(-1/2)*(1 + mahalanobis(cbind(x,y), mu[[j]], Sigma[[j]])/nu)^(-(p+nu)/2)
          dens <- dens + pii[j] * 2*(denst)*pt(sqrt((p + nu)/(mahalanobis(cbind(x,y), mu[[j]], Sigma[[j]]) + nu))*t(lambda[[j]])%*%solve(matrix.sqrt(Sigma[[j]]))%*%(c(x,y) - mu[[j]]), df = nu + p)
        }
        dens
      }
      #x <- seq(min(dat[,1])-1,max(dat[,1])+1, length = slice)
      #y <- seq(min(dat[,2])-15,max(dat[,2])+10, length = slice)
      #x <- seq(min(dat[,1])-1,max(dat[,1])+1, length = slice)     # Swiss
      #y <- seq(min(dat[,2])-1,max(dat[,2])+1, length = slice)     # Swiss
      x <- seq(min(dat[,1])-x.min,max(dat[,1])+x.max, length = slice)    
      y <- seq(min(dat[,2])-y.min,max(dat[,2])+y.max, length = slice)  
      f <- matrix(0,slice,slice)
      for (i in 1:slice) for (j in 1:slice) f[i,j] <-  mixed.ST(x[i], y[j], model$pii, model$mu, model$Sigma, model$shape, model$nu)
      #contour(x, y, f, nlevels = 5, levels = c(0.1, 0.015, 0.005, 0.0009, 0.00015))
      contour(x, y, f, , nlevels = ncontour, ...)
      #contour(x, y, f, , levels = c(0.3,0.15, 0.075,0.0375,0.01875)) # Swiss
      points(dat[,1], dat[,2], col = (model$group+1))
   }

   if (class(model) == "Skew.cn"){
      mixed.SNC <- function(x, y, pii, mu, Sigma, lambda, nu) {
#        n <- nrow(y)
#        p <- ncol(y)
        dens <- 0
        for (j in 1:g) dens <- dens + pii[j]* 2*(nu[1]*dmvnorm(cbind(x,y), mu[[j]], Sigma[[j]]/nu[2])*pnorm(sqrt(nu[2])*t(lambda[[j]])%*%solve(matrix.sqrt(Sigma[[j]]))%*%(c(x,y) - mu[[j]]) ) + (1 - nu[1])*dmvnorm(cbind(x,y), mu[[j]], Sigma[[j]])*pnorm(t(lambda[[j]])%*%solve(matrix.sqrt(Sigma[[j]]))%*%(c(x,y) - mu[[j]])) )
        dens
      }
#      xc <- yc <- 0
#      for (j in 1:g) {
#        xc <- xc + model$mu[[j]][1]
#        yc <- yc + model$mu[[j]][2]
#      }
#      xc <- xc / g
#      yc <- yc / g
#      x <- seq(xc - x.min, xc + x.max, length = slice)
#      y <- seq(yc - y.min, yc + y.max, length = slice)
#      f <- matrix(0,slice,slice)
      x <- seq(min(dat[,1])-x.min,max(dat[,1])+x.max, length = slice)    
      y <- seq(min(dat[,2])-y.min,max(dat[,2])+y.max, length = slice)  
      f <- matrix(0,slice,slice)
      for (i in 1:slice) for (j in 1:slice) f[i,j] <-  mixed.SNC(x[i], y[j], model$pii, model$mu, model$Sigma, model$shape, model$nu)
      contour(x, y, f, nlevels = ncontour, ...)
      points(dat[,1], dat[,2], col = (model$group+1))
   }

   if (class(model) == "Skew.slash"){
      #x <- seq(min(dat[,1])-1,max(dat[,1])+1, length = slice)
      #y <- seq(min(dat[,2])-15,max(dat[,2])+10, length = slice)
      #x <- seq(min(dat[,1])-1,max(dat[,1])+1, length = slice)     # Swiss
      #y <- seq(min(dat[,2])-1,max(dat[,2])+1, length = slice)     # Swiss
      x <- seq(min(dat[,1])-x.min,max(dat[,1])+x.max, length = slice)    
      y <- seq(min(dat[,2])-y.min,max(dat[,2])+y.max, length = slice)  
      f <- matrix(0,slice,slice)
      for (i in 1:slice) for (j in 1:slice) f[i,j] <-  d.mixedmvSS(cbind(x[i],y[j]), model$pii, model$mu, model$Sigma, model$shape, model$nu)
      #contour(x, y, f, nlevels = 5, levels = c(0.1, 0.015, 0.005, 0.0009, 0.00015))
      contour(x, y, f, , nlevels = ncontour, ...)
      #contour(x, y, f, , levels = c(0.3,0.15, 0.075,0.0375,0.01875)) # Swiss
      points(dat[,1], dat[,2], col = (model$group+1))
   }
   title(main=paste("Contour plot for",class(model)), font = 2)


}

###########################################################
####### Matriz de informacao caso multivariado  ###########
adjoint <- function(A) det(A)*solve(A)
deriv.der <- function(A,B,C) det(A)*sum(B * t(C))

imm.smsn <- function(y, model){
  if((class(model) != "t") && (class(model) != "Skew.t") && (class(model) != "Skew.cn") && (class(model) != "Skew.slash") && (class(model) != "Skew.normal") && (class(model) != "Normal")) stop(paste("Family",class(model),"not recognized.",sep=" "))
  if (ncol(y) <= 1) stop(paste("The dimension of y (p) is: ", ncol(y),". We need p >= 2.",sep=" "))
  if (model$uni.Gama) stop("Sorry. The information matrix cannot be calculated when the uni.Gama was used!")
  y <- as.matrix(y)
  dimnames(y) <- NULL
  n <- nrow(y)
  p <- ncol(y)
  g <- length(model$pii)
  Sipi <- Simu <- Silambda <- c()
  Ssigma <- c()
  
  for(i in 1:length(model$Sigma)) model$Sigma[[i]] <- model$Sigma[[i]] %*% model$Sigma[[i]]    

  mu <- model$mu
  Sigma <- model$Sigma
  lambda <- model$shape
  pii <- model$pii
  nu <- model$nu

  if (class(model) == "t"){
    soma <- soma2 <- 0
    
    I.Phi <- function(w=0,Ai=NULL,di,nu=0) as.numeric((( 2^w*nu^(nu/2)*gamma(w + nu/2))/(gamma(nu/2)*(nu + di)^(nu/2 + w)))*pt( ((Ai)/(di + nu)^(0.5))*sqrt(2*w + nu), 2*w + nu))
    I.phi <- function(w=0,Ai=NULL,di,nu=0) as.numeric(((2^w*nu^(nu/2))/(sqrt(2*pi)*gamma(nu/2)))*(1/(di + Ai^2 + nu))^((nu + 2*w)/2)*gamma((nu + 2*w)/2))
    
    for (i in 1:n){
      S <- c() # vetor com todas as derivadas em relacao a cada parametro desconhecido do modelo
      dPsi.dnu <- 0
      yi <- matrix(y[i,], 1, p)#2)
      for (j in 1:g){
        Dr <- matrix.sqrt(Sigma[[j]])
        Dr.inv <- solve(Dr)
        d.sig <- det(Sigma[[j]])        
        
        Ai <- as.numeric(t(lambda[[j]])%*%Dr.inv%*%(y[i,] - mu[[j]]))
        di <- as.numeric(mahalanobis(yi, mu[[j]], Sigma[[j]]))
        

        #derivadinhas       
        dir.dmu <- -2*(Dr.inv%*%Dr.inv)%*%(y[i,] - mu[[j]])
        dAir.dmu <- -Dr.inv%*%lambda[[j]]
        #dAir.dlambda <- solve(Dr)%*%(y[i,] - mu[[j]])
      
        dPsi.dmu <- ((2*d.sig^(-1/2))/(2*pi)^(p/2))*( dAir.dmu * I.phi((p+1)/2, Ai, di, nu) - (1/2)*dir.dmu*I.Phi((p/2)+1, Ai, di, nu) )
        #dPsi.dlambda <- ((2*det(Sigma[[j]])^(-1/2))/(2*pi)^(p/2))*dAir.dlambda*I.phi((p+1)/2, Ai, di, nu)

        #para os elementos de sigma                      
        l <- m <- 1
        for(k in 1:((p+1)*p/2)) {
           Vis <- FALSE
           D <- matrix(rep(0,p*p),p,p)
           D[l,m] <- D[m,l] <- 1
           
           #ddet.ds <- -(1/det(Dr)^2)*Dadj[l,m]
           ddet.ds <- -(1/det(Dr)^2)*deriv.der(Dr,Dr.inv,D)
           dir.ds <- - t(y[i,] - mu[[j]])%*%Dr.inv%*%(D%*%Dr.inv + Dr.inv%*%D)%*%Dr.inv%*%(y[i,] - mu[[j]])
           dAir.ds <- - t(lambda[[j]])%*%Dr.inv%*%D%*%Dr.inv%*%(y[i,] - mu[[j]])

           dPsi.dsigma <- (2/(2*pi)^(p/2))*( ddet.ds*I.Phi(p/2, Ai, di, nu) - (1/2)*dir.ds*d.sig^(-1/2)*I.Phi(p/2+1, Ai, di, nu) + d.sig^(-1/2)*dAir.ds*I.phi((p+1)/2, Ai, di, nu) )           
           Ssigma[k] <- (pii[j]/ d.mixedmvST(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dsigma

           if(((l*m - p*floor((l*m)/p)) == 0) && (l != m)) {
              l <- l+1
              m <- l
              Vis <- TRUE
           }
           if(!Vis) m <- m+1
        }


        ui <- rgamma(10000, shape = nu/2, rate = nu/2)
        resto <- mean(ui^(p/2)*log(ui)*exp(-ui*di/2)*pnorm(ui^(1/2)*Ai))
        dPsi.dnu <- dPsi.dnu + pii[j]*((d.sig^(-1/2))/(2*pi)^(p/2))*((log(nu/2)+1-digamma(nu/2))*I.Phi(p/2, Ai, di, nu) - I.Phi((p+2)/2, Ai, di, nu) + resto)

        Simu <- as.vector((pii[j]/ d.mixedmvST(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dmu)
        #Silambda <- as.vector((pii[j]/ d.mixedmvST(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dlambda )

      
        S <- c(S, Simu, Ssigma)
      }
      Sinu <- (1/d.mixedmvST(yi, pii, mu, Sigma, lambda, nu))*dPsi.dnu
      if(g>1){
        for(j in 1:(g-1)) Sipi[j] <- (1/d.mixedmvST(yi, pii, mu, Sigma, lambda, nu))*( dmvt.ls(yi, mu[[j]], Sigma[[j]], lambda[[j]], nu) - dmvt.ls(yi, mu[[g]], Sigma[[g]], lambda[[g]], nu))        
        S <- c(S, Sipi, Sinu)
      }
      if(g == 1) S <- c(S, Sinu)
      soma <- soma + S%*%t(S)
 #     soma2 <- soma2 + S
    }
    NAME <- piiN <- c()
    for(i in 1:g){
        SigmaN <- muN <- shapeN <- c()
        for (k in 1:p){
           muN <- c(muN,paste("mu",i,"_",k,sep=""))
           #shapeN <- c(shapeN,"RETIRA");#c(shapeN,paste("shape",i,"_",k,sep=""))
        }
        l <- m <- 1
        for(k in 1:((p+1)*p/2)) {
           Vis <- FALSE
           SigmaN <- c(SigmaN,paste("Sigma",i,"_",l,m,sep=""))
           if(((l*m - p*floor((l*m)/p)) == 0) && (l != m)) {
              l <- l+1
              m <- l
              Vis <- TRUE
           }
           if(!Vis) m <- m+1
        }
       NAME <- c(NAME,muN,shapeN,SigmaN)
    }
    
    for(i in 1:(g-1)){
       piiN <- c(piiN, paste("pii",i,sep=""))
       NAME <- c(NAME,piiN,"nu")
    }
    if( g==1) NAME <- c(NAME,"nu")
    dimnames(soma)[[1]] <- NAME
    dimnames(soma)[[2]] <- NAME   

  }
  
  if (class(model) == "Skew.t"){
    soma <- soma2 <- 0
    
    I.Phi <- function(w=0,Ai=NULL,di,nu=0) as.numeric((( 2^w*nu^(nu/2)*gamma(w + nu/2))/(gamma(nu/2)*(nu + di)^(nu/2 + w)))*pt( ((Ai)/(di + nu)^(0.5))*sqrt(2*w + nu), 2*w + nu))
    I.phi <- function(w=0,Ai=NULL,di,nu=0) as.numeric(((2^w*nu^(nu/2))/(sqrt(2*pi)*gamma(nu/2)))*(1/(di + Ai^2 + nu))^((nu + 2*w)/2)*gamma((nu + 2*w)/2))
    
    for (i in 1:n){
      S <- c() # vetor com todas as derivadas em relacao a cada parametro desconhecido do modelo
      dPsi.dnu <- 0
      yi <- matrix(y[i,], 1, p)#2)
      for (j in 1:g){
        Dr <- matrix.sqrt(Sigma[[j]])
        Dr.inv <- solve(Dr)
        d.sig <- det(Sigma[[j]])        
        
        Ai <- as.numeric(t(lambda[[j]])%*%Dr.inv%*%(y[i,] - mu[[j]]))
        di <- as.numeric(mahalanobis(yi, mu[[j]], Sigma[[j]]))
        
        #derivadinhas
        dir.dmu <- -2*(Dr.inv%*%Dr.inv)%*%(y[i,] - mu[[j]])
        dAir.dmu <- -Dr.inv%*%lambda[[j]]
        dAir.dlambda <- Dr.inv%*%(y[i,] - mu[[j]])
      
        dPsi.dmu <- ((2*d.sig^(-1/2))/(2*pi)^(p/2))*( dAir.dmu * I.phi((p+1)/2, Ai, di, nu) - (1/2)*dir.dmu*I.Phi((p/2)+1, Ai, di, nu) )
        dPsi.dlambda <- ((2*d.sig^(-1/2))/(2*pi)^(p/2))*dAir.dlambda*I.phi((p+1)/2, Ai, di, nu)

        #para os elementos de sigma                      
        l <- m <- 1
        for(k in 1:((p+1)*p/2)) {
           Vis <- FALSE
           D <- matrix(rep(0,p*p),p,p)
           D[l,m] <- D[m,l] <- 1

           ddet.ds <- -(1/det(Dr)^2)*deriv.der(Dr,Dr.inv,D)
           dir.ds <- - t(y[i,] - mu[[j]])%*%Dr.inv%*%(D%*%Dr.inv + Dr.inv%*%D)%*%Dr.inv%*%(y[i,] - mu[[j]])
           dAir.ds <- - t(lambda[[j]])%*%Dr.inv%*%D%*%Dr.inv%*%(y[i,] - mu[[j]])

           dPsi.dsigma <- (2/(2*pi)^(p/2))*( ddet.ds*I.Phi(p/2, Ai, di, nu) - (1/2)*dir.ds*d.sig^(-1/2)*I.Phi(p/2+1, Ai, di, nu) + d.sig^(-1/2)*dAir.ds*I.phi((p+1)/2, Ai, di, nu) )           
           Ssigma[k] <- (pii[j]/ d.mixedmvST(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dsigma

           if(((l*m - p*floor((l*m)/p)) == 0) && (l != m)) {
              l <- l+1
              m <- l
              Vis <- TRUE
           }
           if(!Vis) m <- m+1
        }


        ui <- rgamma(10000, shape = nu/2, rate = nu/2)
        resto <- mean(ui^(p/2)*log(ui)*exp(-ui*di/2)*pnorm(ui^(1/2)*Ai))
        dPsi.dnu <- dPsi.dnu + pii[j]*((d.sig^(-1/2))/(2*pi)^(p/2))*((log(nu/2)+1-digamma(nu/2))*I.Phi(p/2, Ai, di, nu) - I.Phi((p+2)/2, Ai, di, nu) + resto)

        Simu <- as.vector((pii[j]/ d.mixedmvST(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dmu)
        Silambda <- as.vector((pii[j]/ d.mixedmvST(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dlambda )
      
        S <- c(S, Simu, Silambda, Ssigma)
      }
      Sinu <- (1/d.mixedmvST(yi, pii, mu, Sigma, lambda, nu))*dPsi.dnu
      if(g>1){
        for(j in 1:(g-1)) Sipi[j] <- (1/d.mixedmvST(yi, pii, mu, Sigma, lambda, nu))*( dmvt.ls(yi, mu[[j]], Sigma[[j]], lambda[[j]], nu) - dmvt.ls(yi, mu[[g]], Sigma[[g]], lambda[[g]], nu))        
        S <- c(S, Sipi, Sinu)
      }
      if(g == 1) S <- c(S, Sinu)
      soma <- soma + S%*%t(S)
 #     soma2 <- soma2 + S
    }
    NAME <- piiN <- c()
    for(i in 1:g){
        SigmaN <- muN <- shapeN <- c()
        for (k in 1:p){
           muN <- c(muN,paste("mu",i,"_",k,sep=""))
           shapeN <- c(shapeN,paste("shape",i,"_",k,sep=""))
        }
        l <- m <- 1
        for(k in 1:((p+1)*p/2)) {
           Vis <- FALSE
           SigmaN <- c(SigmaN,paste("Sigma",i,"_",l,m,sep=""))
           if(((l*m - p*floor((l*m)/p)) == 0) && (l != m)) {
              l <- l+1
              m <- l
              Vis <- TRUE
           }
           if(!Vis) m <- m+1
        }
       NAME <- c(NAME,muN,shapeN,SigmaN)
    }
    
    for(i in 1:(g-1)){
       piiN <- c(piiN, paste("pii",i,sep=""))
       NAME <- c(NAME,piiN,"nu")
    }
    if( g==1) NAME <- c(NAME,"nu")
    dimnames(soma)[[1]] <- NAME
    dimnames(soma)[[2]] <- NAME   

  }

  if (class(model) == "Skew.cn"){
    soma <- soma2 <- 0

    I.Phi <- function(w=0,Ai=NULL,di,nu=0) as.numeric( sqrt(2*pi)*(nu[1]*nu[2]^(w -0.5)*dnorm(sqrt(di), 0, sqrt(1/nu[2]))*pnorm(nu[2]^(1/2)*Ai) + (1 - nu[1])*(dnorm(sqrt(di), 0,1)*pnorm(Ai)) )   )
    I.phi <- function(w=0,Ai=NULL,di,nu=0) as.numeric( nu[1]*nu[2]^(w - 0.5)*dnorm(sqrt(di + Ai^2), 0, sqrt(1/nu[2])) + (1 - nu[1])*dnorm(sqrt(di + Ai^2))   )

    for (i in 1:n){
      S <- c() # vetor com todas as derivadas em relacao a cada parametro desconhecido do modelo
      dPsi.dnu1 <- dPsi.dnu2 <- 0
      yi <- matrix(y[i,], 1, p)#2)
      for (j in 1:g){
        Dr <- matrix.sqrt(Sigma[[j]])
        Dr.inv <- solve(Dr)
        d.sig <- det(Sigma[[j]])        
        
        Ai <- as.numeric(t(lambda[[j]])%*%Dr.inv%*%(y[i,] - mu[[j]]))
        di <- as.numeric(mahalanobis(yi, mu[[j]], Sigma[[j]]))
        
        #derivadinhas
        dir.dmu <- -2*(Dr.inv%*%Dr.inv)%*%(y[i,] - mu[[j]])
        dAir.dmu <- -Dr.inv%*%lambda[[j]]
        dAir.dlambda <- Dr.inv%*%(y[i,] - mu[[j]])
      
        dPsi.dmu <- ((2*d.sig^(-1/2))/(2*pi)^(p/2))*( dAir.dmu * I.phi((p+1)/2, Ai, di, nu) - (1/2)*dir.dmu*I.Phi((p/2)+1, Ai, di, nu) )
        dPsi.dlambda <- ((2*d.sig^(-1/2))/(2*pi)^(p/2))*dAir.dlambda*I.phi((p+1)/2, Ai, di, nu)

        dPsi.dnu1 <- dPsi.dnu1 + pii[j]*2*(dmvnorm(yi, mu[[j]], nu[2]^(-1)*Sigma[[j]])*pnorm(nu[2]^(1/2)*Ai) - dmvnorm(yi, mu[[j]], Sigma[[j]])*pnorm(Ai) )
        dPsi.dnu2 <- dPsi.dnu1 + pii[j]*((nu[1]*d.sig^(-1/2)*nu[2]^(p/2))/(2*pi)^(p/2))*exp(-nu[2]*di/2)*(p*nu[2]^(-1)*pnorm(nu[2]^(1/2)*Ai) + dnorm(nu[2]^(1/2)*Ai)*Ai*nu[2]^(-1/2) - pnorm(nu[2]^(1/2)*Ai)*di )

        #para os elementos de sigma
        l <- m <- 1
        for(k in 1:((p+1)*p/2)) {
           Vis <- FALSE
           D <- matrix(rep(0,p*p),p,p)
           D[l,m] <- D[m,l] <- 1

           ddet.ds <- -(1/det(Dr)^2)*deriv.der(Dr,Dr.inv,D)
           dir.ds <- - t(y[i,] - mu[[j]])%*%Dr.inv%*%(D%*%Dr.inv + Dr.inv%*%D)%*%Dr.inv%*%(y[i,] - mu[[j]])
           dAir.ds <- - t(lambda[[j]])%*%Dr.inv%*%D%*%Dr.inv%*%(y[i,] - mu[[j]])

           dPsi.dsigma <- (2/(2*pi)^(p/2))*( ddet.ds*I.Phi(p/2, Ai, di, nu) - (1/2)*dir.ds*d.sig^(-1/2)*I.Phi(p/2+1, Ai, di, nu) + d.sig^(-1/2)*dAir.ds*I.phi((p+1)/2, Ai, di, nu) )           
           Ssigma[k] <- (pii[j]/ d.mixedmvSNC(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dsigma

           if(((l*m - p*floor((l*m)/p)) == 0) && (l != m)) {
              l <- l+1
              m <- l
              Vis <- TRUE
           }
           if(!Vis) m <- m+1
        }


        Simu <- as.vector((pii[j]/ d.mixedmvSNC(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dmu)
        Silambda <- as.vector((pii[j]/ d.mixedmvSNC(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dlambda )

        S <- c(S, Simu, Silambda, Ssigma)
      }
      Sinu1 <- (1/d.mixedmvSNC(yi, pii, mu, Sigma, lambda, nu))*dPsi.dnu1
      Sinu2 <- (1/d.mixedmvSNC(yi, pii, mu, Sigma, lambda, nu))*dPsi.dnu2
      if (g >1){ 
        for(j in 1:(g-1)) Sipi[j] <- (1/d.mixedmvSNC(yi, pii, mu, Sigma, lambda, nu))*( dmvSNC(yi, mu[[j]], Sigma[[j]], lambda[[j]], nu) - dmvSNC(yi, mu[[g]], Sigma[[g]], lambda[[g]], nu))
        S <- c(S, Sipi, Sinu1, Sinu2)
      }
      if( g==1) S <- c(S, Sinu1, Sinu2)
      soma <- soma + S%*%t(S)
#      soma2 <- soma2 + S
    }
    NAME <- piiN <- c()
    for(i in 1:g){
        SigmaN <- muN <- shapeN <- c()
        for (k in 1:p){
           muN <- c(muN,paste("mu",i,"_",k,sep=""))
           shapeN <- c(shapeN,paste("shape",i,"_",k,sep=""))
        }
        l <- m <- 1
        for(k in 1:((p+1)*p/2)) {
           Vis <- FALSE
           SigmaN <- c(SigmaN,paste("Sigma",i,"_",l,m,sep=""))
           if(((l*m - p*floor((l*m)/p)) == 0) && (l != m)) {
              l <- l+1
              m <- l
              Vis <- TRUE
           }
           if(!Vis) m <- m+1
        }
       NAME <- c(NAME,muN,shapeN,SigmaN)
    }
    
    if ( g>1){
      for(i in 1:(g-1)) piiN <- c(piiN, paste("pii",i,sep=""))
       NAME <- c(NAME,piiN,"nu1","nu2")
    }
    if( g==1) NAME <- c(NAME,"nu1","nu2")
    dimnames(soma)[[1]] <- NAME
    dimnames(soma)[[2]] <- NAME   

  }
  
    if (class(model) == "Skew.slash"){
    soma <- soma2 <- 0

    I.Phi <- function(w=0,Ai=NULL,di,nu=0) {
      Esper <- vector(mode = "numeric", length = length(di))
      for(i in 1:length(di)){
        U <- runif(2500)
        V <- pgamma(1,w + nu, di[i]/2)*U
        S <- qgamma(V,w + nu, di[i]/2)
        Esper[i] <- mean(pnorm(S^(1/2)*Ai[i]))
      }
      res1 <-  (nu*(2^(w + nu)*gamma(w + nu))/(di^(w + nu)))*pgamma(1, w + nu, di/2)*Esper
      return(res1)
    }

    I.phi <- function(w=0,Ai=NULL,di,nu=0){
      res2 <- ((nu*2^(w + nu)*gamma(w + nu))/(sqrt(2*pi)*(di + Ai^2)^(w + nu)))*pgamma(1, w + nu, (di + Ai^2)/2)
      return(res2)
    }

    for (i in 1:n){
      S <- c() # vetor com todas as derivadas em relacao a cada parametro desconhecido do modelo
      dPsi.dnu <- 0
      yi <- matrix(y[i,], 1, p)#2)
      for (j in 1:g){
        Dr <- matrix.sqrt(Sigma[[j]])
        Dr.inv <- solve(Dr)
        d.sig <- det(Sigma[[j]])        
        
        Ai <- as.numeric(t(lambda[[j]])%*%Dr.inv%*%(y[i,] - mu[[j]]))
        di <- as.numeric(mahalanobis(yi, mu[[j]], Sigma[[j]]))
        
        #derivadinhas
        dir.dmu <- -2*(Dr.inv%*%Dr.inv)%*%(y[i,] - mu[[j]])
        dAir.dmu <- -Dr.inv%*%lambda[[j]]
        dAir.dlambda <- Dr.inv%*%(y[i,] - mu[[j]])
      
        dPsi.dmu <- ((2*d.sig^(-1/2))/(2*pi)^(p/2))*( dAir.dmu * I.phi((p+1)/2, Ai, di, nu) - (1/2)*dir.dmu*I.Phi((p/2)+1, Ai, di, nu) )
        dPsi.dlambda <- ((2*d.sig^(-1/2))/(2*pi)^(p/2))*dAir.dlambda*I.phi((p+1)/2, Ai, di, nu)

        #para os elementos de sigma
        l <- m <- 1
        for(k in 1:((p+1)*p/2)) {
           Vis <- FALSE
           D <- matrix(rep(0,p*p),p,p)
           D[l,m] <- D[m,l] <- 1

           ddet.ds <- -(1/det(Dr)^2)*deriv.der(Dr,Dr.inv,D)
           dir.ds <- - t(y[i,] - mu[[j]])%*%Dr.inv%*%(D%*%Dr.inv + Dr.inv%*%D)%*%Dr.inv%*%(y[i,] - mu[[j]])
           dAir.ds <- - t(lambda[[j]])%*%Dr.inv%*%D%*%Dr.inv%*%(y[i,] - mu[[j]])

           dPsi.dsigma <- (2/(2*pi)^(p/2))*( ddet.ds*I.Phi(p/2, Ai, di, nu) - (1/2)*dir.ds*d.sig^(-1/2)*I.Phi(p/2+1, Ai, di, nu) + d.sig^(-1/2)*dAir.ds*I.phi((p+1)/2, Ai, di, nu) )           
           Ssigma[k] <- (pii[j]/ d.mixedmvSS(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dsigma

           if(((l*m - p*floor((l*m)/p)) == 0) && (l != m)) {
              l <- l+1
              m <- l
              Vis <- TRUE
           }
           if(!Vis) m <- m+1
        }
       
        #ui <- rgamma(10000, shape = nu/2, rate = nu/2)
        #resto <- mean(ui^(p/2)*log(ui)*exp(-ui*di/2)*pnorm(ui^(1/2)*Ai))
        #dPsi.dnu <- dPsi.dnu + pii[j]*((det(Sigma[[j]])^(-1/2))/(2*pi)^(p/2))*((log(nu/2)+1-4*digamma(nu/2))*I.Phi(p/2, Ai, di, nu) - I.Phi((p+2)/2, Ai, di, nu) + resto)

        u <- runif(8000)
        dPsi.dnu <- dPsi.dnu + pii[j]*mean(2*u^(nu - 1)*(1 + nu*log(u))*( (det(Sigma[[j]])^(-1/2)/(2*pi)^(p/2))*exp(-(u^(-1)/2)*di) )*pnorm(u^(1/2)*as.numeric(t(lambda[[j]])%*%solve(matrix.sqrt(Sigma[[j]]))%*%t(yi-mu[[j]]))))#integrate(f,0,1)$value

        Simu <- as.vector((pii[j]/ d.mixedmvSS(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dmu)
        Silambda <- as.vector((pii[j]/ d.mixedmvSS(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dlambda )

        S <- c(S, Simu, Silambda, Ssigma)
      }
      Sinu <- (1/d.mixedmvSS(yi, pii, mu, Sigma, lambda, nu))*dPsi.dnu
      if(g>1){
        for(j in 1:(g-1)) Sipi[j] <- (1/d.mixedmvSS(yi, pii, mu, Sigma, lambda, nu))*( dmvSS(yi, mu[[j]], Sigma[[j]], lambda[[j]], nu) - dmvSS(yi, mu[[g]], Sigma[[g]], lambda[[g]], nu))
        S <- c(S, Sipi, Sinu)
      }
      if( g==1) S <- c(S, Sinu)
      soma <- soma + S%*%t(S)
#      soma2 <- soma2 + S
    }
    NAME <- piiN <- c()
    for(i in 1:g){
        SigmaN <- muN <- shapeN <- c()
        for (k in 1:p){
           muN <- c(muN,paste("mu",i,"_",k,sep=""))
           shapeN <- c(shapeN,paste("shape",i,"_",k,sep=""))
        }
        l <- m <- 1
        for(k in 1:((p+1)*p/2)) {
           Vis <- FALSE
           SigmaN <- c(SigmaN,paste("Sigma",i,"_",l,m,sep=""))
           if(((l*m - p*floor((l*m)/p)) == 0) && (l != m)) {
              l <- l+1
              m <- l
              Vis <- TRUE
           }
           if(!Vis) m <- m+1
        }
       NAME <- c(NAME,muN,shapeN,SigmaN)
    }
    
    if( g > 1){
       for(i in 1:(g-1)) piiN <- c(piiN, paste("pii",i,sep=""))
       NAME <- c(NAME,piiN,"nu")
    }
    if( g==1) NAME <- c(NAME,"nu")
    dimnames(soma)[[1]] <- NAME
    dimnames(soma)[[2]] <- NAME   
  }

  if (class(model) == "Skew.normal"){
    soma <- soma2 <- 0

    I.Phi <- function(w=0,Ai=NULL,di,nu=0) as.numeric( exp(-di/2)*pnorm(Ai) )
    I.phi <- function(w=0,Ai=NULL,di,nu=0) as.numeric( exp(-di/2)*dnorm(Ai) )

    for (i in 1:n){
      S <- c() # vetor com todas as derivadas em relacao a cada parametro desconhecido do modelo
      yi <- matrix(y[i,], 1, p)#2)
      for (j in 1:g){
        Dr <- matrix.sqrt(Sigma[[j]])
        Dr.inv <- solve(Dr)
        d.sig <- det(Sigma[[j]])        
        
        Ai <- as.numeric(t(lambda[[j]])%*%Dr.inv%*%(y[i,] - mu[[j]]))
        di <- as.numeric(mahalanobis(yi, mu[[j]], Sigma[[j]]))
        
        #derivadinhas
        dir.dmu <- -2*(Dr.inv%*%Dr.inv)%*%(y[i,] - mu[[j]])
        dAir.dmu <- -Dr.inv%*%lambda[[j]]
        dAir.dlambda <- Dr.inv%*%(y[i,] - mu[[j]])
      
        dPsi.dmu <- ((2*d.sig^(-1/2))/(2*pi)^(p/2))*( dAir.dmu * I.phi(Ai=Ai, di=di) - (1/2)*dir.dmu*I.Phi(Ai=Ai, di=di) )
        dPsi.dlambda <- ((2*d.sig^(-1/2))/(2*pi)^(p/2))*dAir.dlambda*I.phi(Ai=Ai, di=di)

        #para os elementos de sigma
        l <- m <- 1
        for(k in 1:((p+1)*p/2)) {
           Vis <- FALSE
           D <- matrix(rep(0,p*p),p,p)
           D[l,m] <- D[m,l] <- 1

           ddet.ds <- -(1/det(Dr)^2)*deriv.der(Dr,Dr.inv,D)
           dir.ds <- - t(y[i,] - mu[[j]])%*%Dr.inv%*%(D%*%Dr.inv + Dr.inv%*%D)%*%Dr.inv%*%(y[i,] - mu[[j]])
           dAir.ds <- - t(lambda[[j]])%*%Dr.inv%*%D%*%Dr.inv%*%(y[i,] - mu[[j]])

           dPsi.dsigma <- (2/(2*pi)^(p/2))*( ddet.ds*I.Phi(Ai=Ai, di=di) - (1/2)*dir.ds*d.sig^(-1/2)*I.Phi(Ai=Ai, di=di) + d.sig^(-1/2)*dAir.ds*I.phi(Ai=Ai, di=di) )           
           Ssigma[k] <- (pii[j]/ d.mixedmvSN(yi, pii, mu, Sigma, lambda) )*dPsi.dsigma

           if(((l*m - p*floor((l*m)/p)) == 0) && (l != m)) {
              l <- l+1
              m <- l
              Vis <- TRUE
           }
           if(!Vis) m <- m+1
        }

        Simu <- as.vector((pii[j]/ d.mixedmvSN(yi, pii, mu, Sigma, lambda) )*dPsi.dmu)
        Silambda <- as.vector((pii[j]/ d.mixedmvSN(yi, pii, mu, Sigma, lambda) )*dPsi.dlambda )

        S <- c(S, Simu, Silambda, Ssigma)
      }      
      if(g>1){
        for(j in 1:(g-1)) Sipi[j] <- (1/d.mixedmvSN(yi, pii, mu, Sigma, lambda))*( dmvSN(yi, mu[[j]], Sigma[[j]], lambda[[j]]) - dmvSN(yi, mu[[g]], Sigma[[g]], lambda[[g]]))    
        S <- c(S, Sipi)
      }
      soma <- soma + S%*%t(S)
#      soma2 <- soma2 + S
    }
    NAME <- piiN <- c()
    for(i in 1:g){
        SigmaN <- muN <- shapeN <- c()
        for (k in 1:p){
           muN <- c(muN,paste("mu",i,"_",k,sep=""))
           shapeN <- c(shapeN,paste("shape",i,"_",k,sep=""))
        }
        l <- m <- 1
        for(k in 1:((p+1)*p/2)) {
           Vis <- FALSE
           SigmaN <- c(SigmaN,paste("Sigma",i,"_",l,m,sep=""))
           if(((l*m - p*floor((l*m)/p)) == 0) && (l != m)) {
              l <- l+1
              m <- l
              Vis <- TRUE
           }
           if(!Vis) m <- m+1
        }
       NAME <- c(NAME,muN,shapeN,SigmaN)
    }
    
    if(g >1){
       for(i in 1:(g-1)) piiN <- c(piiN, paste("pii",i,sep=""))
       NAME <- c(NAME,piiN)
    }
    dimnames(soma)[[1]] <- NAME
    dimnames(soma)[[2]] <- NAME   

#    dimnames(soma)[[1]] <- c("mu1_1","mu1_2","shape1_1","shape1_2","Sigma1_11","Sigma1_12","Sigma1_22","mu2_1","mu2_2","shape2_1","shape2_2","Sigma2_11","Sigma2_12","Sigma2_22","pii")
#    dimnames(soma)[[2]] <- c("mu1_1","mu1_2","shape1_1","shape1_2","Sigma1_11","Sigma1_12","Sigma1_22","mu2_1","mu2_2","shape2_1","shape2_2","Sigma2_11","Sigma2_12","Sigma2_22","pii")

  }
  
  if (class(model) == "Normal"){
    soma <- soma2 <- 0

    I.Phi <- function(w=0,Ai=NULL,di,nu=0) as.numeric( exp(-di/2)*1/2 )
    I.phi <- function(w=0,Ai=NULL,di,nu=0) as.numeric( exp(-di/2)*dnorm(0) )

    for (i in 1:n){
      S <- c() # vetor com todas as derivadas em relacao a cada parametro desconhecido do modelo
      yi <- matrix(y[i,], 1, p)#2)
      for (j in 1:g){
        Dr <- matrix.sqrt(Sigma[[j]])
        Dr.inv <- solve(Dr)
        d.sig <- det(Sigma[[j]])        
        
        Ai <- as.numeric(t(lambda[[j]])%*%Dr.inv%*%(y[i,] - mu[[j]]))
        di <- as.numeric(mahalanobis(yi, mu[[j]], Sigma[[j]]))
        
        #derivadinhas
        dir.dmu <- -2*(Dr.inv%*%Dr.inv)%*%(y[i,] - mu[[j]])
        dAir.dmu <- -Dr.inv%*%lambda[[j]]
      
        dPsi.dmu <- ((2*d.sig^(-1/2))/(2*pi)^(p/2))*( dAir.dmu * I.phi(di=di) - (1/2)*dir.dmu*I.Phi(di=di) )

        #para os elementos de sigma                      
        l <- m <- 1
        for(k in 1:((p+1)*p/2)) {
           Vis <- FALSE
           D <- matrix(rep(0,p*p),p,p)
           D[l,m] <- D[m,l] <- 1

           ddet.ds <- -(1/det(Dr)^2)*deriv.der(Dr,Dr.inv,D)
           dir.ds <- - t(y[i,] - mu[[j]])%*%Dr.inv%*%(D%*%Dr.inv + Dr.inv%*%D)%*%Dr.inv%*%(y[i,] - mu[[j]])
           dAir.ds <- - t(lambda[[j]])%*%Dr.inv%*%D%*%Dr.inv%*%(y[i,] - mu[[j]])

           dPsi.dsigma <- (2/(2*pi)^(p/2))*( ddet.ds*I.Phi(di=di) - (1/2)*dir.ds*d.sig^(-1/2)*I.Phi(di=di) + d.sig^(-1/2)*dAir.ds*I.phi(di=di) )           
           Ssigma[k] <- (pii[j]/ d.mixedmvSN(yi, pii, mu, Sigma, lambda) )*dPsi.dsigma

           if(((l*m - p*floor((l*m)/p)) == 0) && (l != m)) {
              l <- l+1
              m <- l
              Vis <- TRUE
           }
           if(!Vis) m <- m+1
        }
      
        Simu <- as.vector((pii[j]/ d.mixedmvSN(yi, pii, mu, Sigma, lambda) )*dPsi.dmu)

        S <- c(S, Simu, Ssigma)
      }
      if(g>1){
        for(j in 1:(g-1)) Sipi[j] <- (1/d.mixedmvSN(yi, pii, mu, Sigma, lambda))*( dmvSN(yi, mu[[j]], Sigma[[j]], lambda[[j]]) - dmvSN(yi, mu[[g]], Sigma[[g]], lambda[[g]]))    
        S <- c(S, Sipi)  
      }
      soma <- soma + S%*%t(S)
#      soma2 <- soma2 + S
    }
    NAME <- piiN <- c()
    for(i in 1:g){
        SigmaN <- muN <- c()
        for (k in 1:p) muN <- c(muN,paste("mu",i,"_",k,sep=""))
        l <- m <- 1
        for(k in 1:((p+1)*p/2)) {
           Vis <- FALSE
           SigmaN <- c(SigmaN,paste("Sigma",i,"_",l,m,sep=""))
           if(((l*m - p*floor((l*m)/p)) == 0) && (l != m)) {
              l <- l+1
              m <- l
              Vis <- TRUE
           }
           if(!Vis) m <- m+1
        }
       NAME <- c(NAME,muN,SigmaN)
    }
    if(g > 1){
       for(i in 1:(g-1)) piiN <- c(piiN, paste("pii",i,sep=""))
       NAME <- c(NAME,piiN)
    }
    dimnames(soma)[[1]] <- NAME
    dimnames(soma)[[2]] <- NAME   
#    dimnames(soma)[[1]] <- c("mu1_1","mu1_2","Sigma1_11","Sigma1_12","Sigma1_22","mu2_1","mu2_2","Sigma2_11","Sigma2_12","Sigma2_22","pii")
#    dimnames(soma)[[2]] <- c("mu1_1","mu1_2","Sigma1_11","Sigma1_12","Sigma1_22","mu2_1","mu2_2","Sigma2_11","Sigma2_12","Sigma2_22","pii")

  }

  return(list(IM=soma))  
#  return(list(soma, soma2))
}
