####################   spv function   #####################

spv <- function(design.matrix, design.matrix.2 = NULL, design.matrix.3 = NULL, 
                  des.names = c("Design 1","Design 2","Design 3"),
                  scale = TRUE, add.pts = TRUE, label = "ON"){  
  norm2 <- function(x){return(sqrt(sum(x^2)))}
  shuffle.fun <- function(row.vec){
    num.var <- length(row.vec)
    row.vec[shuffle(num.var)]
  }
  
  n.grid = 15
  n.var <- ncol(design.matrix)
  
  if(n.var > 7){stop("The maximum number of design factors allowed is 7 \n")}
  
  if(scale == TRUE){ 
    temp.radii <- apply(design.matrix, MARGIN = 1, norm2)
    scale.factor <- sqrt(n.var)/max(temp.radii)
    design.matrix <- scale.factor*design.matrix
  }
  
  design.matrix<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix))),byrow=FALSE,ncol=n.var))  
  N.generated <- 15^n.var   # CHANGE  !!!
  model.X <- model.matrix( ~quad(.) , design.matrix)
  nrow <- nrow(model.X)
  ncol <- ncol(model.X)
  M.inv<- nrow*solve(t(model.X)%*%model.X)
  
  # Design 2
  if(!is.null(design.matrix.2)){
    n.var.2 <- ncol(design.matrix.2)
    if(scale == TRUE){ 
      temp.radii.2 <- apply(design.matrix.2, MARGIN = 1, norm2)
      scale.factor.2 <- sqrt(n.var.2)/max(temp.radii.2)
      design.matrix.2 <- scale.factor.2*design.matrix.2
    }
    design.matrix.2<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix.2))),byrow=FALSE,ncol=n.var.2))  
    model.X.2 <- model.matrix( ~quad(.) , design.matrix.2)
    nrow.2 <- nrow(model.X.2)
    ncol.2 <- ncol(model.X.2)
    if(ncol != ncol.2){
      stop("Designs need to have the same number of design factors")
    }
    M.inv.2<- nrow.2*solve(t(model.X.2)%*%model.X.2)   
  } # Design 2
  # Design 3
  if(!is.null(design.matrix.3)){
    
    n.var.3 <- ncol(design.matrix.3)
    if(scale == TRUE){ 
      temp.radii.3 <- apply(design.matrix.3, MARGIN = 1, norm2)
      scale.factor.3 <- sqrt(n.var.3)/max(temp.radii.3)
      design.matrix.3 <- scale.factor.3*design.matrix.3
    }
    design.matrix.3<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix.3))),byrow=FALSE,ncol=n.var.3))  
    model.X.3 <- model.matrix( ~quad(.) , design.matrix.3)
    nrow.3 <- nrow(model.X.3)
    ncol.3 <- ncol(model.X.3)
    if(ncol != ncol.3){
      stop("Designs need to have the same number of design factors")
    }
    M.inv.3<- nrow.3*solve(t(model.X.3)%*%model.X.3)   
  } # Design 3  
  
  
  Matrix.V <- matrix(numeric(0), ncol = 4, nrow = n.grid)
  colnames(Matrix.V) <- c("Radius","  max.V","  min.V"," average.V")
  
  if(!is.null(design.matrix.2)){
    Matrix.V.2 <- matrix(numeric(0), ncol = 4, nrow = n.grid)
    colnames(Matrix.V.2) <- c("Radius","  max.V","  min.V"," average.V")  
  }# Design 2
  
  if(!is.null(design.matrix.3)){
    Matrix.V.3 <- matrix(numeric(0), ncol = 4, nrow = n.grid)
    colnames(Matrix.V.3) <- c("Radius","  max.V","  min.V"," average.V")  
  }# Design 3
  
  
  # Generate points on a sphere
  #set.seed(1972264)
  #n.var=2; N.generated=1000;
  Urand <- rnorm(n=(n.var*N.generated),mean = 0, sd = 1)
  #Urand <- ifelse(Urand < (-sqrt(n.var)+.001), -sqrt(n.var), Urand)
  #Urand <- ifelse(Urand > ( sqrt(n.var)-.001),  sqrt(n.var), Urand)
  # Urand <- ifelse(Urand > -.001 & Urand < .001,  0, Urand)
  
  rand.sphere <- matrix(Urand,byrow= TRUE, ncol=n.var)
  radius <- seq(from = 0, to = sqrt(n.var), length.out = n.grid)
  norm2.rand <- apply(rand.sphere, MARGIN = 1, norm2)
  norm2.rand.rep <- matrix(rep(norm2.rand, each = n.var),byrow= TRUE, ncol=n.var)
  gen.array <- array(rand.sphere/norm2.rand.rep,dim = c(N.generated,n.var,n.grid))
  Array.with.R <- rep(radius,each = (N.generated*n.var))*gen.array
  
  #plot(as.data.frame(Array.with.R[,,15]), pch = ".")
  
  for(ii in 1:n.grid){
    pred.model  <- model.matrix( ~quad(.), as.data.frame(Array.with.R[,,ii]))
    
    Var.pred <- numeric(N.generated)
    for(jj in 1:N.generated){
      each.obs <- as.vector(pred.model[jj, ])
      Var.pred[jj] <- as.numeric(t(each.obs)%*%M.inv%*%each.obs)
    }
    
    minimum <- min(Var.pred)
    maximum <- max(Var.pred)
    average <- mean(Var.pred)
    R.rho <- radius[ii]
    Matrix.V[ii,] <- c(R.rho,maximum,minimum,average)  
    
    # Design 2
    if(!is.null(design.matrix.2)){
      Var.pred.2 <- numeric(N.generated)
      for(jj in 1:N.generated){
        each.obs <- as.vector(pred.model[jj, ])
        Var.pred.2[jj] <- as.numeric(t(each.obs)%*%M.inv.2%*%each.obs)
      }
      minimum.2 <- min(Var.pred.2)
      maximum.2 <- max(Var.pred.2)
      average.2 <- mean(Var.pred.2)
      Matrix.V.2[ii,] <- c(R.rho,maximum.2,minimum.2,average.2)     
    } 
    # Design 3
    if(!is.null(design.matrix.3)){
      Var.pred.3 <- numeric(N.generated)
      for(jj in 1:N.generated){
        each.obs <- as.vector(pred.model[jj, ])
        Var.pred.3[jj] <- as.numeric(t(each.obs)%*%M.inv.3%*%each.obs)
      }
      minimum.3 <- min(Var.pred.3)
      maximum.3 <- max(Var.pred.3)
      average.3 <- mean(Var.pred.3)
      Matrix.V.3[ii,] <- c(R.rho,maximum.3,minimum.3,average.3)     
    } 
  }
  
  # Find lim.max
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    lim.max <- max(c(Matrix.V[,2], Matrix.V.2[,2], Matrix.V.3[,2]))
  }else{
    if(!is.null(design.matrix.2)){
      lim.max <- max(c(Matrix.V[,2], Matrix.V.2[,2]))
    }else{ 
      lim.max <- max(Matrix.V[,2])
    }
  }
  # Find lim.min
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    lim.min <- min(c(Matrix.V[,3], Matrix.V.2[,3], Matrix.V.3[,3]))
  }else{
    if(!is.null(design.matrix.2) ){
      lim.min <- min(c(Matrix.V[,3], Matrix.V.2[,3]))
    }else{
      lim.min <- min(Matrix.V[,3])
    }
  }
  
  #windows(width = 5, height = 5)
  #par(mfrow=c(1,2),
  #    mai = c(1, 1, 1, 0.5),
  #    omi = c(0.2, 0.1, 0.1, 0.1),
  #    mgp = c(2, 1, 0),
  #    xpd = FALSE)
  par(mfrow=c(1,1),
      mai = c(0.75, 0.75, 0.75, 0.375),
      omi = c(0.075, 0.0375, 0.0375, 0.0375),
      mgp = c(2, 1, 0),
      xpd = FALSE)
  
  if(label == "OFF"){
    plot(Matrix.V[,1],Matrix.V[,2], type = "l", lwd = 2, lty = 1, 
         ylim = c(lim.min-3, lim.max + 3), col = "#2E2E2E",
         xlab = "Radius", 
         ylab = "Scaled Variance",
         panel.first = grid())
    lines(Matrix.V[,1],Matrix.V[,3], lwd = 2, col = "#2E2E2E", lty = 6)
    lines(Matrix.V[,1],Matrix.V[,4], lwd = 2, col = "#2E2E2E", lty = 3)
  }
  if(label == "ON"){
    plot(Matrix.V[,1],Matrix.V[,2], type = "l", lwd = 2, lty = 1, 
         ylim = c(lim.min-4, lim.max + 6), col = "#2E2E2E",
         xlab = "Radius", 
         ylab = "Scaled Variance",
         panel.first = grid())
    lines(Matrix.V[,1],Matrix.V[,3], lwd = 2, col = "#2E2E2E", lty = 6)
    lines(Matrix.V[,1],Matrix.V[,4], lwd = 2, col = "#2E2E2E", lty = 3)
    legend(x = 0, y = lim.max + 6,  legend = c("Max","Avg","Min"),lty=c(1,3,6),
           lwd=c(2,2,2),col=c(1,1,1),
           inset = 0.05, bg="transparent", bty = "n")
  }
  #help(legend)
  # Design 2
  if(!is.null(design.matrix.2)){
    lines(Matrix.V.2[,1],Matrix.V.2[,2], lwd = 2, col = "#EE4000", lty = 1)
    lines(Matrix.V.2[,1],Matrix.V.2[,3], lwd = 2, col = "#EE4000", lty = 6)
    lines(Matrix.V.2[,1],Matrix.V.2[,4], lwd = 2, col = "#EE4000", lty = 3)
  }
  # Design 3
  if(!is.null(design.matrix.3)){
    lines(Matrix.V.3[,1],Matrix.V.3[,2], lwd = 2, col = "#3A5FCD", lty = 1)
    lines(Matrix.V.3[,1],Matrix.V.3[,3], lwd = 2, col = "#3A5FCD", lty = 6)
    lines(Matrix.V.3[,1],Matrix.V.3[,4], lwd = 2, col = "#3A5FCD", lty = 3)
  }
  
  if(label == "ON"){
    # Put legend Design names
    if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
      legend(x = 0.6, y = lim.max + 6,  
             legend = c(des.names[1],des.names[2],des.names[3]),lty=c(1,1,1),
             lwd=c(2,2,2),col=c("#2E2E2E","#EE4000","#3A5FCD"),
             inset = 0.05, bg="transparent", bty = "n")
    }else{
      if(!is.null(design.matrix.2)){
        legend(x = 0.6, y = lim.max + 6,  
               legend = c(des.names[1],des.names[2]),lty=c(1,1),
               lwd=c(2,2),col=c("#2E2E2E","#EE4000"),
               inset = 0.05, bg="transparent", bty = "n")
      }else{
        legend(x = 0.6, y = lim.max + 6,  
               legend = c(des.names[1]),lty=c(1),
               lwd=c(2),col=c("#2E2E2E"),
               inset = 0.05, bg="transparent", bty = "n")
        
      }}
    
    # abline
    abline(h = ncol, col = "#8B8682",  lwd = 1)
    text(sqrt(n.var)-0.15, lim.min -3 , paste("p =",ncol),col = "#8B8682")
  }
  
  
  if(add.pts == TRUE){
    # Filling points within max and min
    
    if(n.var <= 4){
      N.of.pts.gen <- N.generated
      set.seed(1972264)
      #help(rnorm)
      Urand2 <- rnorm(n=(n.var*N.generated), mean = 0, sd = 1 )        
      #Urand2 <- runif(n=(4*n.var*N.generated),min=-sqrt(n.var),max=sqrt(n.var))
      ##Urand2 <- ifelse(Urand2 < (-sqrt(n.var)+.001), -sqrt(n.var), Urand2)
      #Urand2 <- ifelse(Urand2 > ( sqrt(n.var)-.001),  sqrt(n.var), Urand2)
      #Urand2 <- ifelse(Urand2 > -.001 & Urand2 < .001,  0, Urand2)
      
      rand.sphere.2times <- matrix(Urand2,byrow= TRUE, ncol=n.var)
      norm2.rand.2times <- apply(rand.sphere.2times, MARGIN = 1, norm2)
      norm2.rand.rep.2times <- matrix(rep(norm2.rand.2times, each = n.var),byrow= TRUE, ncol=n.var)   
      
      Urand3 <-      runif(n = N.generated, min = 0, max = 1)
      radius.random    <- sqrt(n.var)*Urand3^(1/n.var)
      
      #Urand3 <- runif(n = 4*N.generated, min = 0, max = sqrt(n.var))
      #Urand3 <- ifelse(Urand3 < .001, 0, Urand3)
      #Urand3 <- ifelse(Urand3 > ( sqrt(n.var)-.001),  sqrt(n.var), Urand3)  
      #radius.random      <- Urand3
      Surface.pts.random <- radius.random*as.data.frame(rand.sphere.2times/norm2.rand.rep.2times)
      pred.model.random  <- model.matrix( ~quad(.), Surface.pts.random)
      
    }
    
    if(n.var > 4){
      N.of.pts.gen <- N.generated
      set.seed(1972264)
      
      Urand2 <- rnorm(n=(n.var*N.generated), mean = 0, sd = 1 )        
      #Urand2 <- runif(n=(2*n.var*N.generated), min=-sqrt(n.var),max=sqrt(n.var))
      #Urand2 <- ifelse(Urand2 < (-sqrt(n.var)+.001), -sqrt(n.var), Urand2)
      #Urand2 <- ifelse(Urand2 > ( sqrt(n.var)-.001),  sqrt(n.var), Urand2)
      #Urand2 <- ifelse(Urand2 > -.001 & Urand2 < .001,  0, Urand2)
      
      
      rand.sphere.2times <- matrix(Urand2,byrow= TRUE, ncol=n.var)
      norm2.rand.2times <- apply(rand.sphere.2times, MARGIN = 1, norm2)
      norm2.rand.rep.2times <- matrix(rep(norm2.rand.2times, each = n.var),byrow= TRUE, ncol=n.var)   
      
      #Urand3 <- runif(n = 2*N.generated, min = 0, max = sqrt(n.var))
      #Urand3 <- ifelse(Urand3 < .001, 0, Urand3)
      #Urand3 <- ifelse(Urand3 > ( sqrt(n.var)-.001),  sqrt(n.var), Urand3)
      #radius.random      <- Urand3 
      
      Urand3 <-  runif(n = N.generated, min = 0, max = 1)
      radius.random  <- sqrt(n.var)*Urand3^(1/n.var)
      
      Surface.pts.random <- radius.random*as.data.frame(rand.sphere.2times/norm2.rand.rep.2times)
      pred.model.random  <- model.matrix( ~quad(.), Surface.pts.random)
    }
    
    Var.pred.random<- numeric(N.of.pts.gen)
    for(kk in 1:N.of.pts.gen){
      each.obs <- as.vector(pred.model.random[kk, ])
      Var.pred.random[kk] <- as.numeric(t(each.obs)%*%M.inv%*%each.obs)
    }
    points(radius.random,Var.pred.random, pch = ".", col = "#2E2E2E")
    
    # Filling points within max and min for Design #2
    if(!is.null(design.matrix.2)){
      Var.pred.random.2<- numeric(N.of.pts.gen)
      for(kk in 1:N.of.pts.gen){
        each.obs <- as.vector(pred.model.random[kk, ])
        Var.pred.random.2[kk] <- as.numeric(t(each.obs)%*%M.inv.2%*%each.obs)
      }
      points(radius.random,Var.pred.random.2, pch = ".", col = "#EE4000")
    }
    # Filling points within max and min for Design #3
    if(!is.null(design.matrix.3)){
      Var.pred.random.3<- numeric(N.of.pts.gen)
      for(kk in 1:N.of.pts.gen){
        each.obs <- as.vector(pred.model.random[kk, ])
        Var.pred.random.3[kk] <- as.numeric(t(each.obs)%*%M.inv.3%*%each.obs)
      }
      points(radius.random,Var.pred.random.3, pch = ".", col = "#3A5FCD")
    }
  }
  
  par(mfrow=c(1,1))
  
  # Return Values
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    return(list(design.1 = Matrix.V, design.2 = Matrix.V.2, 
                design.3 = Matrix.V.3))
  }else{
    if(!is.null(design.matrix.2)){
      return(list(design.1 = Matrix.V, design.2 = Matrix.V.2))
    }else{
      return(Matrix.V)      
    }
    
  }
}

######################### fds.sphere function  ########################


fds.sphere <- function(design.matrix, design.matrix.2 = NULL, design.matrix.3 = NULL, 
                         des.names = c("Design 1","Design 2","Design 3"),
                         scale = TRUE, label = "ON"){
  
  norm2 <- function(x){return(sqrt(sum(x^2)))}
  shuffle.fun <- function(row.vec){
    num.var <- length(row.vec)
    row.vec[shuffle(num.var)]
  }
  
  n.var <- ncol(design.matrix)
  
  if(n.var > 7){stop("The maximum number of design factors allowed is 7 \n")}
  
  
  if(scale == TRUE){ 
    temp.radii <- apply(design.matrix, MARGIN = 1, norm2)
    scale.factor <- sqrt(n.var)/max(temp.radii)
    design.matrix <- scale.factor*design.matrix
  }
  
  design.matrix<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix))),byrow=FALSE,ncol=n.var))  
  N.generated <- 60^n.var
  model.X <- model.matrix( ~quad(.) , design.matrix)
  nrow <- nrow(model.X)
  ncol <- ncol(model.X)
  M.inv<- nrow*solve(t(model.X)%*%model.X)
  
  # Design 2
  if(!is.null(design.matrix.2)){
    n.var.2 <- ncol(design.matrix.2)
    if(scale == TRUE){ 
      temp.radii.2 <- apply(design.matrix.2, MARGIN = 1, norm2)
      scale.factor.2 <- sqrt(n.var.2)/max(temp.radii.2)
      design.matrix.2 <- scale.factor.2*design.matrix.2
    }
    design.matrix.2<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix.2))),byrow=FALSE,ncol=n.var.2))  
    model.X.2 <- model.matrix( ~quad(.) , design.matrix.2)
    nrow.2 <- nrow(model.X.2)
    ncol.2 <- ncol(model.X.2)
    if(ncol != ncol.2){
      stop("Designs need to have the same number of design factors")
    }
    M.inv.2<- nrow.2*solve(t(model.X.2)%*%model.X.2)   
  } # Design 2
  # Design 3
  if(!is.null(design.matrix.3)){
    
    n.var.3 <- ncol(design.matrix.3)
    if(scale == TRUE){ 
      temp.radii.3 <- apply(design.matrix.3, MARGIN = 1, norm2)
      scale.factor.3 <- sqrt(n.var.3)/max(temp.radii.3)
      design.matrix.3 <- scale.factor.3*design.matrix.3
    }
    design.matrix.3<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix.3))),byrow=FALSE,ncol=n.var.3))  
    model.X.3 <- model.matrix( ~quad(.) , design.matrix.3)
    nrow.3 <- nrow(model.X.3)
    ncol.3 <- ncol(model.X.3)
    if(ncol != ncol.3){
      stop("Designs need to have the same number of design factors")
    }
    M.inv.3<- nrow.3*solve(t(model.X.3)%*%model.X.3)   
  } # Design 3  
  
  # Filling points within max and min
  
  if(n.var <= 4){
    N.of.pts.gen <- N.generated
    # set.seed(1972264)
    
    Urand <- rnorm(n=(n.var*N.generated), mean = 0, sd = 1)
    #Urand <- runif(n=(4*n.var*N.generated),min=-sqrt(n.var),max=sqrt(n.var))
    #Urand <- ifelse(Urand < (-sqrt(n.var)+.001), -sqrt(n.var), Urand)
    #Urand <- ifelse(Urand > ( sqrt(n.var)-.001),  sqrt(n.var), Urand)
    #Urand <- ifelse(Urand > -.001 & Urand < .001,  0, Urand)
    
    rand.sphere.2times <- matrix(Urand,byrow= TRUE, ncol=n.var)
    norm2.rand.2times <- apply(rand.sphere.2times, MARGIN = 1, norm2)
    norm2.rand.rep.2times <- matrix(rep(norm2.rand.2times, each = n.var),byrow= TRUE, ncol=n.var) 
    
    #Urand2 <- runif(n = N.of.pts.gen, min = 0, max = sqrt(n.var))
    #Urand2 <- ifelse(Urand2 < .001, 0, Urand2)
    #Urand2 <- ifelse(Urand2 > ( sqrt(n.var)-.001),  sqrt(n.var), Urand2)
    Urand2 <-      runif(n = N.of.pts.gen, min = 0, max = 1)
    radius.random      <- sqrt(n.var)*Urand2^(1/n.var)
    
    #radius.random      <- Urand2
    Surface.pts.random <- radius.random*as.data.frame(rand.sphere.2times/norm2.rand.rep.2times)
    pred.model.random  <- model.matrix( ~quad(.), Surface.pts.random)
  }
  
  if(n.var > 4){
    N.of.pts.gen <- N.generated
    set.seed(1972264)
    
    Urand <- rnorm(n=(n.var*N.generated), mean = 0, sd = 1)
    #Urand <- runif(n=(3*n.var*N.generated),min=-sqrt(n.var),max=sqrt(n.var))
    #Urand <- ifelse(Urand < (-sqrt(n.var)+.001), -sqrt(n.var), Urand)
    #Urand <- ifelse(Urand > ( sqrt(n.var)-.001),  sqrt(n.var), Urand)
    #Urand <- ifelse(Urand > -.001 & Urand < .001,  0, Urand)
    
    rand.sphere.2times <- matrix(Urand,byrow= TRUE, ncol=n.var)
    norm2.rand.2times <- apply(rand.sphere.2times, MARGIN = 1, norm2)
    norm2.rand.rep.2times <- matrix(rep(norm2.rand.2times, each = n.var),byrow= TRUE, ncol=n.var) 
    
    #Urand2 <- runif(n = N.of.pts.gen, min = 0, max = sqrt(n.var))
    #Urand2 <- ifelse(Urand2 < .001, 0, Urand2)
    #Urand2 <- ifelse(Urand2 > ( sqrt(n.var)-.001),  sqrt(n.var), Urand2)
    Urand2 <-      runif(n = N.of.pts.gen, min = 0, max = 1)
    radius.random      <- sqrt(n.var)*Urand2^(1/n.var)
    
    #radius.random      <-  Urand2
    Surface.pts.random <- radius.random*as.data.frame(rand.sphere.2times/norm2.rand.rep.2times)
    pred.model.random  <- model.matrix( ~ quad(.), Surface.pts.random)
  }
  
  Var.pred.random<- numeric(N.of.pts.gen)
  for(kk in 1:N.of.pts.gen){
    each.obs <- as.vector(pred.model.random[kk, ])
    Var.pred.random[kk] <- as.numeric(t(each.obs)%*%M.inv%*%each.obs)
  }
  
  # Design #2
  if(!is.null(design.matrix.2)){
    Var.pred.random.2<- numeric(N.of.pts.gen)
    for(kk in 1:N.of.pts.gen){
      each.obs <- as.vector(pred.model.random[kk, ])
      Var.pred.random.2[kk] <- as.numeric(t(each.obs)%*%M.inv.2%*%each.obs)
    }
  }
  #Design #3
  if(!is.null(design.matrix.3)){
    Var.pred.random.3<- numeric(N.of.pts.gen)
    for(kk in 1:N.of.pts.gen){
      each.obs <- as.vector(pred.model.random[kk, ])
      Var.pred.random.3[kk] <- as.numeric(t(each.obs)%*%M.inv.3%*%each.obs)
    }    
  }
  
  # Find lim.max
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    lim.max <- max(c(Var.pred.random, Var.pred.random.2, Var.pred.random.3))
  }else{
    if(!is.null(design.matrix.2)){
      lim.max <- max(c(Var.pred.random, Var.pred.random.2))
    }else{ 
      lim.max <- max(Var.pred.random)
    }
  }
  # Find lim.min
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    lim.min <- min(c(Var.pred.random, Var.pred.random.2, Var.pred.random.3))
  }else{
    if(!is.null(design.matrix.2)){
      lim.min <- min(c(Var.pred.random, Var.pred.random.2))
    }else{ 
      lim.min <- min(Var.pred.random)
    }
  }
  
  #windows(width = 5, height = 5)
  #par(mfrow=c(1,2),
  #    mai = c(1, 1, 1, 0.5),
  #    omi = c(0.2, 0.1, 0.1, 0.1),
  #    mgp = c(2, 1, 0),
  #    xpd = FALSE)
  par(mfrow=c(1,1),
      mai = c(0.85, 0.75, 0.75, 0.375),
      omi = c(0.075, 0.0375, 0.0375, 0.0375),
      mgp = c(2, 1, 0),
      xpd = FALSE)
  #   par(mfrow=c(1,1),
  #       mai = c(0.75, 0.75, 0.75, 0.375),
  #       omi = c(0.075, 0.0375, 0.0375, 0.0375),
  #       mgp = c(2, 1, 0),
  #       xpd = FALSE)
  # Make FDS plots
  
  if(label == "ON"){
    plot(seq(0, 1, 0.01),quantile(Var.pred.random, probs = seq(0, 1, 0.01), type = 4),
         type = "l",lwd = 2, col = "#2E2E2E",
         xlab = "Fraction of design space",ylab = "Scaled Variance",
         ylim = c(lim.min-4, lim.max + 6), panel.first = grid())
    
    # Make FDS plots for Design 2
    if(!is.null(design.matrix.2)){
      lines(seq(0, 1, 0.01),quantile(Var.pred.random.2, probs = seq(0, 1, 0.01), type = 4), lwd = 2, col = "#EE4000")
    }
    # Make FDS plots for Design 3
    if(!is.null(design.matrix.3)){
      lines(seq(0, 1, 0.01),quantile(Var.pred.random.3, probs = seq(0, 1, 0.01), type = 4), lwd = 2, col = "#3A5FCD")
    }
    
    # Put legend in VDG plots
    if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
      legend(x = 0, y = lim.max + 6,  
             legend = c(des.names[1],des.names[2],des.names[3]),lty=c(1,1,1),
             lwd=c(2,2,2),col=c("#2E2E2E","#EE4000","#3A5FCD"),
             inset = 0.05, bg="transparent", bty = "n")
    }else{
      if(!is.null(design.matrix.2)){
        legend(x = 0, y = lim.max + 6,  
               legend = c(des.names[1],des.names[2]),lty=c(1,1),
               lwd=c(2,2),col=c("#2E2E2E","#EE4000"),
               inset = 0.05, bg="transparent", bty = "n")
      }else{
        legend(x = 0, y = lim.max + 6,  
               legend = c(des.names[1]),lty=c(1),
               lwd=c(2),col=c("#2E2E2E"),
               inset = 0.05, bg="transparent", bty = "n")
      }} 
  }
  
  if(label == "OFF"){
    plot(seq(0, 1, 0.01),quantile(Var.pred.random, probs = seq(0, 1, 0.01), type = 4),
         type = "l",lwd = 2, col = "#2E2E2E",
         xlab = "Fraction of design space",ylab = "Scaled Variance",
         ylim = c(lim.min-3, lim.max + 3), panel.first = grid())
    
    # Make FDS plots for Design 2
    if(!is.null(design.matrix.2)){
      lines(seq(0, 1, 0.01),quantile(Var.pred.random.2, probs = seq(0, 1, 0.01), type = 4), lwd = 2, col = "#EE4000")
    }
    # Make FDS plots for Design 3
    if(!is.null(design.matrix.3)){
      lines(seq(0, 1, 0.01),quantile(Var.pred.random.3, probs = seq(0, 1, 0.01), type = 4), lwd = 2, col = "#3A5FCD")
    }
  }
  
  
  
  # abline
  abline(h = ncol, col = "#8B8682",  lwd = 1)
  text(sqrt(n.var)-0.15, lim.min -3 , paste("p =",ncol),col = "#8B8682")
  
}

#########################  cpv function ####################

cpv<- function(design.matrix, design.matrix.2 = NULL, des.names = c("Design 1","Design 2"), add.pts = TRUE){
  
  norm2 <- function(x){return(sqrt(sum(x^2)))}
  shuffle.fun <- function(row.vec){
    num.var <- length(row.vec)
    row.vec[shuffle(num.var)]
  }
  
  n.var <- ncol(design.matrix)
  if(n.var > 6){stop("The maximum number of design factors allowed is 6 \n")}
  critical.R <- numeric(0)
  
  for(ii in 0:n.var){
    critical.R <- append(critical.R,sqrt(ii))
  }
  LR <- length(critical.R)
  
  if(n.var > 4){
    for(jj in 1:(LR-1)){
      critical.R <- 
        append(critical.R,
               seq(from = critical.R[jj], to = critical.R[jj+1], length.out = 3))
    }
  }else{
    for(jj in 1:(LR-1)){
      critical.R <- 
        append(critical.R,
               seq(from = critical.R[jj], to = critical.R[jj+1], length.out = 8))
    }
  }
  
  critical.R<- unique(sort(critical.R))
  
  if(n.var <= 4){
    n.grid <- length(critical.R)
    N.generated <- 10000
    model.X <- model.matrix( ~quad(.) , design.matrix)
    nrow <- nrow(model.X)
    ncol <- ncol(model.X)
    M.inv<- nrow*solve(t(model.X)%*%model.X)
    
    Matrix.V <- matrix(numeric(0), ncol = 4, nrow = n.grid)
    colnames(Matrix.V) <- c("Radius","  max.V","  min.V"," average.V")
    
    if(!is.null(design.matrix.2)){
      Matrix.V.2 <- matrix(numeric(0), ncol = 4, nrow = n.grid)
      colnames(Matrix.V.2) <- c("Radius","  max.V","  min.V"," average.V")  
    }# Design 2
    
    # Design 2
    if(!is.null(design.matrix.2)){
      n.grid <- length(critical.R)
      N.generated <- 10000 
      model.X.2 <- model.matrix( ~quad(.) , design.matrix.2)
      nrow.2 <- nrow(model.X.2)
      ncol.2 <- ncol(model.X.2)
      if(ncol != ncol.2){
        stop("Designs need to have the same number of design factors")
      }
      M.inv.2<- nrow.2*solve(t(model.X.2)%*%model.X.2)   
    } # Design 2
    
    
    # Generate points on a sphere
    set.seed(1234567)
    
    Urand <- rnorm(n=(n.var*N.generated),mean = 0, sd = 1)
    #Urand <- runif(n=(n.var*N.generated),min=-sqrt(n.var),max=sqrt(n.var))
    #Urand <- ifelse(Urand < (-sqrt(n.var)+.001), -sqrt(n.var), Urand)
    #Urand <- ifelse(Urand > ( sqrt(n.var)-.001),  sqrt(n.var), Urand)
    #Urand <- ifelse(Urand > -.001 & Urand < .001,  0, Urand)
    
    rand.sphere <- matrix(Urand,byrow= TRUE, ncol=n.var)
    norm2.rand <- apply(rand.sphere, MARGIN = 1, norm2)
    norm2.rand.rep <- matrix(rep(norm2.rand, each = n.var),byrow= TRUE, ncol=n.var)
    generated.pts<- as.data.frame(rand.sphere/norm2.rand.rep)
    
    count <- 0     
    for(R in critical.R){
      count <- count + 1
      Surface.pts <- R*generated.pts
      Conbind.col<- cbind(Surface.pts,rowSums(abs(Surface.pts[,])>1))
      Surface.pts <- Conbind.col[which(Conbind.col[,(n.var+1)] == 0),]
      Surface.pts <- Surface.pts[,-(n.var+1)]    
      
      nrow.surf <- dim(Surface.pts)[1]
      Var.pred <- numeric(nrow.surf)
      pred.model  <- model.matrix( ~quad(.), Surface.pts)
      if( dim(Surface.pts)[1] == 0){ break }
      
      for(kk in 1:nrow.surf){
        each.obs <- as.vector(pred.model[kk, ])
        Var.pred[kk] <- as.numeric(t(each.obs)%*%M.inv%*%each.obs)
      }
      
      minimum <- min(Var.pred)
      maximum <- max(Var.pred)
      average <- mean(Var.pred)
      Matrix.V[count,] <- c(R,maximum,minimum,average) 
      
      # Design 2
      if(!is.null(design.matrix.2)){
        
        Var.pred.2 <- numeric(nrow.surf)
        for(kk in 1:nrow.surf){
          each.obs <- as.vector(pred.model[kk, ])
          Var.pred.2[kk] <- as.numeric(t(each.obs)%*%M.inv.2%*%each.obs)
        }
        
        minimum.2 <- min(Var.pred.2)
        maximum.2 <- max(Var.pred.2)
        average.2 <- mean(Var.pred.2)
        Matrix.V.2[count,] <- c(R,maximum.2 ,minimum.2 ,average.2) 
        
      }       
      
    } 
    
    Matrix.V<- matrix(Matrix.V[is.na(Matrix.V) == FALSE],byrow = FALSE, ncol = 4)
    colnames(Matrix.V) <- c("Radius","  max.V","  min.V"," average.V")
    
    # Design 2
    if(!is.null(design.matrix.2)){
      Matrix.V.2<- matrix(Matrix.V.2[is.na(Matrix.V.2) == FALSE],byrow = FALSE, ncol = 4)
      colnames(Matrix.V.2) <- c("Radius","  max.V","  min.V"," average.V")      
    }
    ### End Design 2
    R <- critical.R[n.grid]
    Surface.pts<- as.data.frame(matrix(sample(c(-1,1),size = n.var*50,replace = TRUE),byrow = TRUE, ncol = n.var))
    
    Npt.extream <- dim(Surface.pts)[1]
    pred.model  <- model.matrix( ~quad(.), Surface.pts)
    
    
    Var.pred <- numeric(Npt.extream)
    for(kk in 1:Npt.extream){
      each.obs <- as.vector(pred.model[kk, ])
      Var.pred[kk] <- as.numeric(t(each.obs)%*%M.inv%*%each.obs)
    }
    
    minimum <- min(Var.pred)
    maximum <- max(Var.pred)
    average <- mean(Var.pred)
    Matrix.V<- rbind(Matrix.V,c(R,maximum,minimum,average))
    
    # Design 2
    if(!is.null(design.matrix.2)){
      Var.pred.2 <- numeric(Npt.extream)
      for(kk in 1:Npt.extream){
        each.obs <- as.vector(pred.model[kk, ])
        Var.pred.2[kk] <- as.numeric(t(each.obs)%*%M.inv.2%*%each.obs)
      }
      
      minimum.2 <- min(Var.pred.2)
      maximum.2 <- max(Var.pred.2)
      average.2 <- mean(Var.pred.2)
      Matrix.V.2<- rbind(Matrix.V.2, c(R,maximum.2,minimum.2,average.2))     
    }
    ### End Design 2
    
    
    #lim.max <- max(Matrix.V[,2],na.rm = TRUE)
    #lim.min <- min(Matrix.V[,3],na.rm = TRUE)
    
    # Find lim.max
    lim.max <- max(Matrix.V[,2])
    
    # Find lim.max
    if(!is.null(design.matrix.2)){
      lim.max <- max(c(Matrix.V[,2], Matrix.V.2[,2]))
    }else{ 
      lim.max <- max(Matrix.V[,2])
    }
    # Find lim.min
    if(!is.null(design.matrix.2)){
      lim.min <- min(c(Matrix.V[,3], Matrix.V.2[,3]))
    }else{ 
      lim.min <- min(Matrix.V[,3])
    }
    
    
    
    #windows(width = 5, height = 5)
    par(mfrow=c(1,1),
        mai = c(0.75, 0.75, 0.75, 0.375),
        omi = c(0.075, 0.0375, 0.0375, 0.0375),
        mgp = c(2, 1, 0),
        xpd = FALSE)
    
    plot(Matrix.V[,1],Matrix.V[,2], type = "l", lwd = 2, lty = 4, 
         ylim = c(lim.min-4, lim.max + 12), col = "#2E2E2E",
         xlab = "Radius", 
         ylab = "Scaled  Variance",
         panel.first = grid())
    lines(Matrix.V[,1],Matrix.V[,3], lwd = 2, col = "#2E2E2E", lty = 2)
    legend(x = 0, y = lim.max + 13,  legend = c("Max","Min"),lty=c(4,2),
           lwd=c(2,2),col=c(1,1),
           inset = 0.05, bg="transparent", bty = "n")
    abline(h = ncol, col = "#8B8682",  lwd = 1)
    text(sqrt(n.var)-0.15, lim.min -3 , paste("p =",ncol),col = "#8B8682")
    
    # Design 2
    if(!is.null(design.matrix.2)){
      lines(Matrix.V.2[,1],Matrix.V.2[,2], lwd = 2, col = "#EE4000", lty = 4)
      lines(Matrix.V.2[,1],Matrix.V.2[,3], lwd = 2, col = "#EE4000", lty = 2)
    }
    # Filling points
    if(add.pts == TRUE){
      set.seed(1972264)
      
      Urand <- runif(n = N.generated*n.var, min = -1, max = 1)
      Urand <- ifelse(Urand < (-1+.001), -1, Urand)
      Urand <- ifelse(Urand > ( 1-.001),  1, Urand)
      Urand <- ifelse(Urand > -.001 & Urand < .001,  0, Urand)
      
      Cuboidal.pts.random  <- matrix(Urand, byrow = TRUE, ncol = n.var)      
      Cuboidal.pts.random  <- as.data.frame(Cuboidal.pts.random)
      pred.model.random  <- model.matrix( ~quad(.), Cuboidal.pts.random)
      radius.random<- apply(Cuboidal.pts.random, MARGIN = 1, norm2)
      
      size.pred.rand <- dim(pred.model.random)[1]
      Var.pred.random <- numeric(size.pred.rand)
      for(kk in 1:size.pred.rand){
        each.obs <- as.vector(pred.model.random[kk, ])
        Var.pred.random[kk] <- as.numeric(t(each.obs)%*%M.inv%*%each.obs)
      }
      points(radius.random,Var.pred.random, pch = ".", col = "#2E2E2E")
      
      if(!is.null(design.matrix.2)){
        Var.pred.random.2 <- numeric(size.pred.rand)
        for(kk in 1:size.pred.rand){
          each.obs <- as.vector(pred.model.random[kk, ])
          Var.pred.random.2[kk] <- as.numeric(t(each.obs)%*%M.inv.2%*%each.obs)
        }
        points(radius.random,Var.pred.random.2, pch = ".", col = "#EE4000")
      }   
      
    }
    
    # Put designs' name
    if(!is.null(design.matrix.2)){
      legend(x = 0.4, y = lim.max + 13,  
             legend = c(des.names[1],des.names[2]),lty=c(1,1),
             lwd=c(2,2),col=c("#2E2E2E","#EE4000"),
             inset = 0.05, bg="transparent", bty = "n")
    }else{
      legend(x = 0.4, y = lim.max + 13,  
             legend = c(des.names[1]),lty=c(1),
             lwd=c(2),col=c("#2E2E2E"),
             inset = 0.05, bg="transparent", bty = "n")
      
    }
    # Return Values          
    if(!is.null(design.matrix.2)){
      return(list(design.1 = Matrix.V, design.2 = Matrix.V.2))
    }else{
      return(Matrix.V)      
    }       
    
    
  }
  if(n.var > 4 && n.var <= 6){
    n.grid <- length(critical.R)
    N.generated <- 15000
    model.X <- model.matrix( ~quad(.) , design.matrix)
    nrow <- nrow(model.X)
    ncol <- ncol(model.X)
    M.inv<- nrow*solve(t(model.X)%*%model.X)
    
    Matrix.V <- matrix(NA, ncol = n.grid, nrow = N.generated)
    
    # Generate points on a sphere
    set.seed(1234567)
    
    Urand <- rnorm(n=(n.var*N.generated), mean = 0, sd = 1)
    #Urand <- runif(n=(n.var*N.generated),min=-sqrt(n.var),max=sqrt(n.var))
    #Urand <- ifelse(Urand < (-sqrt(n.var)+.001), -sqrt(n.var), Urand)
    #Urand <- ifelse(Urand > ( sqrt(n.var)-.001),  sqrt(n.var), Urand)
    #Urand <- ifelse(Urand > -.001 & Urand < .001,  0, Urand)
    
    rand.sphere <- matrix(Urand, byrow= TRUE, ncol=n.var)
    norm2.rand <- apply(rand.sphere, MARGIN = 1, norm2)
    norm2.rand.rep <- matrix(rep(norm2.rand, each = n.var),byrow= TRUE, ncol=n.var)
    generated.pts<- as.data.frame(rand.sphere/norm2.rand.rep)
    
    count <- 0 
    for(R in critical.R[1:(n.grid - 1)]){
      count <- count + 1
      Surface.pts <- R*generated.pts
      Conbind.col<- cbind(Surface.pts,rowSums(abs(Surface.pts[,])>1))
      Surface.pts <- Conbind.col[which(Conbind.col[,(n.var+1)] == 0),]
      Surface.pts <- Surface.pts[,-(n.var+1)]    
      
      #if( dim(Surface.pts)[1] == 0){ break }
      pred.model  <- model.matrix( ~quad(.), Surface.pts)
      
      size.pred.model <- dim(pred.model)[1]
      Answer <- numeric(size.pred.model)
      for(kk in 1:size.pred.model){
        each.obs <- as.vector(pred.model[kk, ])
        Answer[kk] <- as.numeric(t(each.obs)%*%M.inv%*%each.obs)
      }
      
      Matrix.V[1:length(Answer),count]    <- Answer
    }
    
    R <- critical.R[n.grid]
    Surface.pts<- as.data.frame(matrix(sample(c(-1,1),size = n.var*50,replace = TRUE),byrow = TRUE, ncol = n.var))
    
    pred.model  <- model.matrix( ~quad(.), Surface.pts)
    
    size.pred.m <- dim(pred.model)[1]
    Answer <- numeric(size.pred.m)
    for(kk in 1:size.pred.m){
      each.obs <- as.vector(pred.model[kk, ])
      Answer[kk] <- as.numeric(t(each.obs)%*%M.inv%*%each.obs)
    }
    
    Matrix.V[1:length(Answer),n.grid]    <- Answer
    #windows(width = 5, height = 5)
    par(mfrow=c(1,1),
        mai = c(0.75, 0.75, 0.75, 0.375),
        omi = c(0.075, 0.0375, 0.0375, 0.0375),
        mgp = c(2, 1, 0),
        xpd = FALSE)
    boxplot(Matrix.V, use.cols = TRUE, outline = TRUE, range = 9.5, col = "lightgray",
            whisklty = 1, staplelty = 0, ylab = "Scaled  Variace",
            names = round(critical.R,2), xlab = "Radius", lwd = 2)
  }
  
}

########################### fds.cube function  ########################

fds.cube <- function(design.matrix, design.matrix.2 = NULL, design.matrix.3 = NULL, 
                       des.names = c("Design 1","Design 2","Design 3")){
  
  norm2 <- function(x){return(sqrt(sum(x^2)))}
  shuffle.fun <- function(row.vec){
    num.var <- length(row.vec)
    row.vec[shuffle(num.var)]
  }
  
  n.var <- ncol(design.matrix)
  
  if(n.var > 6){stop("The maximum number of design factors allowed is 6 \n")}
  
  
  design.matrix<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix))),byrow=FALSE,ncol=n.var))  
  model.X <- model.matrix( ~quad(.) , design.matrix)
  nrow <- nrow(model.X)
  ncol <- ncol(model.X)
  M.inv<- nrow*solve(t(model.X)%*%model.X)
  
  # Design 2
  if(!is.null(design.matrix.2)){
    n.var.2 <- ncol(design.matrix.2)
    
    design.matrix.2<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix.2))),byrow=FALSE,ncol=n.var.2))  
    model.X.2 <- model.matrix( ~quad(.) , design.matrix.2)
    nrow.2 <- nrow(model.X.2)
    ncol.2 <- ncol(model.X.2)
    if(ncol != ncol.2){
      stop("Designs need to have the same number of design factors")
    }
    M.inv.2<- nrow.2*solve(t(model.X.2)%*%model.X.2)   
  } # Design 2
  # Design 3
  if(!is.null(design.matrix.3)){
    
    n.var.3 <- ncol(design.matrix.3)
    
    design.matrix.3<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix.3))),byrow=FALSE,ncol=n.var.3))  
    model.X.3 <- model.matrix( ~quad(.) , design.matrix.3)
    nrow.3 <- nrow(model.X.3)
    ncol.3 <- ncol(model.X.3)
    if(ncol != ncol.3){
      stop("Designs need to have the same number of design factors")
    }
    M.inv.3<- nrow.3*solve(t(model.X.3)%*%model.X.3)   
  } # Design 3 
  
  
  if(n.var <= 4){
    set.seed(1972264)
    
    Urand <- runif(n = 7000*n.var, min = -1, max = 1)
    Urand <- ifelse(Urand < (-1+.001), -1, Urand)
    Urand <- ifelse(Urand > ( 1-.001),  1, Urand)
    Urand <- ifelse(Urand > -.001 & Urand < .001,  0, Urand)
    
    Cuboidal.pts.random  <- matrix(Urand, byrow = TRUE, ncol = n.var)      
    Cuboidal.pts.random  <- as.data.frame(Cuboidal.pts.random)
    pred.model.random  <- model.matrix( ~quad(.), Cuboidal.pts.random)
    radius.random<- apply(Cuboidal.pts.random, MARGIN = 1, norm2)
  }
  
  if(n.var > 4){
    set.seed(1972264)
    
    Urand <- runif(n = 12000*n.var, min = -1, max = 1)
    Urand <- ifelse(Urand < (-1+.001), -1, Urand)
    Urand <- ifelse(Urand > ( 1-.001),  1, Urand)
    Urand <- ifelse(Urand > -.001 & Urand < .001,  0, Urand)
    
    Cuboidal.pts.random  <- matrix(Urand, byrow = TRUE, ncol = n.var)      
    Cuboidal.pts.random  <- as.data.frame(Cuboidal.pts.random)
    pred.model.random  <- model.matrix( ~quad(.), Cuboidal.pts.random)
    radius.random<- apply(Cuboidal.pts.random, MARGIN = 1, norm2)
  }
  
  nrow.pred.mo <- dim(pred.model.random)[1]
  Var.pred.random <- numeric(nrow.pred.mo)
  for(kk in 1:nrow.pred.mo){
    each.obs <- as.vector(pred.model.random[kk, ])
    Var.pred.random[kk] <- as.numeric(t(each.obs)%*%M.inv%*%each.obs)
  }
  
  # Design #2
  if(!is.null(design.matrix.2)){
    Var.pred.random.2 <- numeric(nrow.pred.mo)
    for(kk in 1:nrow.pred.mo){
      each.obs <- as.vector(pred.model.random[kk, ])
      Var.pred.random.2[kk] <- as.numeric(t(each.obs)%*%M.inv.2%*%each.obs)
    }
  }
  #Design #3
  if(!is.null(design.matrix.3)){
    Var.pred.random.3 <- numeric(nrow.pred.mo)
    for(kk in 1:nrow.pred.mo){
      each.obs <- as.vector(pred.model.random[kk, ])
      Var.pred.random.3[kk] <- as.numeric(t(each.obs)%*%M.inv.3%*%each.obs)
    }
  }
  
  # Find lim.max
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    lim.max <- max(c(Var.pred.random, Var.pred.random.2, Var.pred.random.3))
  }else{
    if(!is.null(design.matrix.2)){
      lim.max <- max(c(Var.pred.random, Var.pred.random.2))
    }else{ 
      lim.max <- max(Var.pred.random)
    }
  }
  # Find lim.min
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    lim.min <- min(c(Var.pred.random, Var.pred.random.2, Var.pred.random.3))
  }else{
    if(!is.null(design.matrix.2)){
      lim.min <- min(c(Var.pred.random, Var.pred.random.2))
    }else{ 
      lim.min <- min(Var.pred.random)
    }
  }
  
  #windows(width = 5, height = 5)
  par(mfrow=c(1,1),
      mai = c(0.75, 0.75, 0.75, 0.375),
      omi = c(0.075, 0.0375, 0.0375, 0.0375),
      mgp = c(2, 1, 0),
      xpd = FALSE)
  
  # Make FDS plots
  plot(seq(0, 1, 0.01),quantile(Var.pred.random, probs = seq(0, 1, 0.01), type = 4),
       type = "l",lwd = 2, col = "#2E2E2E",
       xlab = "Fraction of design space",ylab = "Scaled Variance",
       ylim = c(lim.min-4, lim.max + 12), panel.first = grid())
  
  # Make FDS plots for Design 2
  if(!is.null(design.matrix.2)){
    lines(seq(0, 1, 0.01),quantile(Var.pred.random.2, probs = seq(0, 1, 0.01), type = 4), lwd = 2, col = "#EE4000")
  }
  # Make FDS plots for Design 3
  if(!is.null(design.matrix.3)){
    lines(seq(0, 1, 0.01),quantile(Var.pred.random.3, probs = seq(0, 1, 0.01), type = 4), lwd = 2, col = "#3A5FCD")
  }
  
  # Put legend in VDG plots
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    legend(x = 0, y = lim.max + 13,  
           legend = c(des.names[1],des.names[2],des.names[3]),lty=c(1,1,1),
           lwd=c(2,2,2),col=c("#2E2E2E","#EE4000","#3A5FCD"),
           inset = 0.05, bg="transparent", bty = "n")
  }else{
    if(!is.null(design.matrix.2)){
      legend(x = 0, y = lim.max + 13,  
             legend = c(des.names[1],des.names[2]),lty=c(1,1),
             lwd=c(2,2),col=c("#2E2E2E","#EE4000"),
             inset = 0.05, bg="transparent", bty = "n")
    }else{
      legend(x = 0, y = lim.max + 13,  
             legend = c(des.names[1]),lty=c(1),
             lwd=c(2),col=c("#2E2E2E"),
             inset = 0.05, bg="transparent", bty = "n")
    }}
  # abline
  abline(h = ncol, col = "#8B8682",  lwd = 1)
  text(0.90, lim.min -3 , paste("p =",ncol),col = "#8B8682")
}

##########################  hyper.vdg #######################

hyperarcs.vdg<- function(design.matrix, design.matrix.2 = NULL, design.matrix.3 = NULL, 
                           des.names = c("Design 1","Design 2","Design 3")){
  
  norm2 <- function(x){return(sqrt(sum(x^2)))}
  shuffle.fun <- function(row.vec){
    num.var <- length(row.vec)
    row.vec[shuffle(num.var)]
  }
  
  n.var <- ncol(design.matrix)
  if(n.var > 6){stop("The maximum number of design factors allowed is 6 \n")}
  n.grid <- 15
  design.matrix<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix))),byrow=FALSE,ncol=n.var))  
  
  model.X <- model.matrix( ~quad(.) , design.matrix)
  nrow <- nrow(model.X)
  ncol <- ncol(model.X)
  M.inv<- nrow*solve(t(model.X)%*%model.X)
  
  # Design 2
  if(!is.null(design.matrix.2)){
    n.var.2 <- ncol(design.matrix.2)
    
    design.matrix.2<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix.2))),byrow=FALSE,ncol=n.var.2))  
    model.X.2 <- model.matrix( ~quad(.) , design.matrix.2)
    nrow.2 <- nrow(model.X.2)
    ncol.2 <- ncol(model.X.2)
    if(ncol != ncol.2){
      stop("Designs need to have the same number of design factors")
    }
    M.inv.2<- nrow.2*solve(t(model.X.2)%*%model.X.2)   
  } # Design 2
  # Design 3
  if(!is.null(design.matrix.3)){
    
    n.var.3 <- ncol(design.matrix.3)
    design.matrix.3<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix.3))),byrow=FALSE,ncol=n.var.3))  
    model.X.3 <- model.matrix( ~quad(.) , design.matrix.3)
    nrow.3 <- nrow(model.X.3)
    ncol.3 <- ncol(model.X.3)
    if(ncol != ncol.3){
      stop("Designs need to have the same number of design factors")
    }
    M.inv.3<- nrow.3*solve(t(model.X.3)%*%model.X.3)   
  } # Design 3  
  
  Matrix.V.A <- matrix(numeric(0), ncol = 4, nrow = n.grid )
  colnames(Matrix.V.A) <- c("Hypercube Radius","  max.V","  min.V"," average.V")
  
  if(!is.null(design.matrix.2)){
    Matrix.V.A.2 <- matrix(numeric(0), ncol = 4, nrow = n.grid)
    colnames(Matrix.V.A.2) <- c("Hypercube Radius","  max.V","  min.V"," average.V")
  }# Design 2
  
  if(!is.null(design.matrix.3)){
    Matrix.V.A.3 <- matrix(numeric(0), ncol = 4, nrow = n.grid)
    colnames(Matrix.V.A.3) <- c("Hypercube Radius","  max.V","  min.V"," average.V")
  }# Design 3
  
  radius <- seq(from = 0, to = 1, length.out = n.grid)
  N.generated <- 25*n.var*2^(n.var)
  
  count<- 0   
  for(R in radius){
    count <- count + 1
    set.seed(1972264)
    mat.shuffle <- as.matrix(shuffleSet(n = n.var, nset = N.generated,check = FALSE))
    sample(c(-R,R), size = 5,replace = TRUE)
    full.basket <- matrix(c(sample(c(-R,R), size = N.generated, replace = TRUE)),
                          ncol = 1, byrow = TRUE)
    if(n.var == 2){
      Urand<- runif(n = N.generated, min = -R, max = R)
      Urand <- ifelse(Urand < (-R+.001), -R, Urand)
      Urand <- ifelse(Urand > ( R-.001),  R, Urand)
      Urand <- ifelse(-.001 < Urand & Urand < .001 ,  0, Urand)  
      mat.before <- cbind(full.basket, Urand)
    }
    if(n.var == 3){
      Urand.1<- runif(n = N.generated, min = -R, max = R)
      Urand.1 <- ifelse(Urand.1 < (-R+.001), -R, Urand.1)
      Urand.1 <- ifelse(Urand.1 > ( R-.001),  R, Urand.1)
      Urand.1 <- ifelse(-.001 < Urand.1 & Urand.1 < .001 ,  0, Urand.1)  
      
      Urand.2<- runif(n = N.generated, min = -R, max = R)
      Urand.2<- ifelse(Urand.2 < (-R+.001), -R, Urand.2)
      Urand.2 <- ifelse(Urand.2 > ( R-.001),  R, Urand.2)
      Urand.2 <- ifelse(-.001 < Urand.2 & Urand.2 < .001 ,  0, Urand.2)  
      
      mat.before <- cbind(full.basket,  Urand.1,  Urand.2)
    }
    if(n.var == 4){
      Urand.1<- runif(n = N.generated, min = -R, max = R)
      Urand.1 <- ifelse(Urand.1 < (-R+.001), -R, Urand.1)
      Urand.1 <- ifelse(Urand.1 > ( R-.001),  R, Urand.1)
      Urand.1 <- ifelse(-.001 < Urand.1 & Urand.1 < .001 ,  0, Urand.1)  
      
      Urand.2<- runif(n = N.generated, min = -R, max = R)
      Urand.2<- ifelse(Urand.2 < (-R+.001), -R, Urand.2)
      Urand.2 <- ifelse(Urand.2 > ( R-.001),  R, Urand.2)
      Urand.2 <- ifelse(-.001 < Urand.2 & Urand.2 < .001 ,  0, Urand.2)  
      
      Urand.3<- runif(n = N.generated, min = -R, max = R)
      Urand.3<- ifelse(Urand.3 < (-R+.001), -R, Urand.3)
      Urand.3 <- ifelse(Urand.3 > ( R-.001),  R, Urand.3)
      Urand.3 <- ifelse(-.001 < Urand.3 & Urand.3 < .001 ,  0, Urand.3)  
      
      mat.before <- cbind(full.basket,Urand.1, Urand.2, Urand.3)
    }
    if(n.var == 5){
      Urand.1<- runif(n = N.generated, min = -R, max = R)
      Urand.1 <- ifelse(Urand.1 < (-R+.001), -R, Urand.1)
      Urand.1 <- ifelse(Urand.1 > ( R-.001),  R, Urand.1)
      Urand.1 <- ifelse(-.001 < Urand.1 & Urand.1 < .001 ,  0, Urand.1)  
      
      Urand.2<- runif(n = N.generated, min = -R, max = R)
      Urand.2<- ifelse(Urand.2 < (-R+.001), -R, Urand.2)
      Urand.2 <- ifelse(Urand.2 > ( R-.001),  R, Urand.2)
      Urand.2 <- ifelse(-.001 < Urand.2 & Urand.2 < .001 ,  0, Urand.2)  
      
      Urand.3<- runif(n = N.generated, min = -R, max = R)
      Urand.3<- ifelse(Urand.3 < (-R+.001), -R, Urand.3)
      Urand.3 <- ifelse(Urand.3 > ( R-.001),  R, Urand.3)
      Urand.3 <- ifelse(-.001 < Urand.3 & Urand.3 < .001 ,  0, Urand.3)  
      
      Urand.4<- runif(n = N.generated, min = -R, max = R)
      Urand.4<- ifelse(Urand.4 < (-R+.001), -R, Urand.4)
      Urand.4 <- ifelse(Urand.4 > ( R-.001),  R, Urand.4)
      Urand.4 <- ifelse(-.001 < Urand.4 & Urand.4 < .001 ,  0, Urand.4)  
      
      mat.before <- cbind(full.basket, Urand.1 , Urand.2 , Urand.3 , Urand.4 )
    }
    if(n.var == 6){
      Urand.1<- runif(n = N.generated, min = -R, max = R)
      Urand.1 <- ifelse(Urand.1 < (-R+.001), -R, Urand.1)
      Urand.1 <- ifelse(Urand.1 > ( R-.001),  R, Urand.1)
      Urand.1 <- ifelse(-.001 < Urand.1 & Urand.1 < .001 ,  0, Urand.1)  
      
      Urand.2<- runif(n = N.generated, min = -R, max = R)
      Urand.2<- ifelse(Urand.2 < (-R+.001), -R, Urand.2)
      Urand.2 <- ifelse(Urand.2 > ( R-.001),  R, Urand.2)
      Urand.2 <- ifelse(-.001 < Urand.2 & Urand.2 < .001 ,  0, Urand.2)  
      
      Urand.3<- runif(n = N.generated, min = -R, max = R)
      Urand.3<- ifelse(Urand.3 < (-R+.001), -R, Urand.3)
      Urand.3 <- ifelse(Urand.3 > ( R-.001),  R, Urand.3)
      Urand.3 <- ifelse(-.001 < Urand.3 & Urand.3 < .001 ,  0, Urand.3)  
      
      Urand.4<- runif(n = N.generated, min = -R, max = R)
      Urand.4<- ifelse(Urand.4 < (-R+.001), -R, Urand.4)
      Urand.4 <- ifelse(Urand.4 > ( R-.001),  R, Urand.4)
      Urand.4 <- ifelse(-.001 < Urand.4 & Urand.4 < .001 ,  0, Urand.4)  
      
      Urand.5<- runif(n = N.generated, min = -R, max = R)
      Urand.5<- ifelse(Urand.5 < (-R+.001), -R, Urand.5)
      Urand.5 <- ifelse(Urand.5 > ( R-.001),  R, Urand.5)
      Urand.5 <- ifelse(-.001 < Urand.5 & Urand.5 < .001 ,  0, Urand.5) 
      
      mat.before <- cbind(full.basket, Urand.1 , Urand.2 , Urand.3 , Urand.4 , Urand.5)
    }
    Cuboidal.pts <- as.data.frame(t(apply(mat.before, MARGIN = 1, shuffle.fun)))
    pred.model  <- model.matrix( ~quad(.), Cuboidal.pts)
    
    N.pt.pred.model <- dim(pred.model)[1]
    Var.pred <- numeric(N.pt.pred.model)
    for(kk in 1:N.pt.pred.model){
      each.obs <- as.vector(pred.model[kk, ])
      Var.pred[kk] <- as.numeric(t(each.obs)%*%M.inv%*%each.obs)
    }
    
    minimum <- min(Var.pred)
    maximum <- max(Var.pred)
    average <- mean(Var.pred)
    Matrix.V.A[count,] <- c(R,maximum,minimum,average) 
    
    ###################################
    # Design 2
    if(!is.null(design.matrix.2)){
      
      Var.pred.2 <- numeric(N.pt.pred.model)
      for(kk in 1:N.pt.pred.model){
        each.obs <- as.vector(pred.model[kk, ])
        Var.pred.2[kk] <- as.numeric(t(each.obs)%*%M.inv.2%*%each.obs)
      }
      
      minimum.2 <- min(Var.pred.2)
      maximum.2 <- max(Var.pred.2)
      average.2 <- mean(Var.pred.2)
      Matrix.V.A.2[count,] <- c(R,maximum.2,minimum.2,average.2)     
    } 
    # Design 3
    if(!is.null(design.matrix.3)){
      
      Var.pred.3 <- numeric(N.pt.pred.model)
      for(kk in 1:N.pt.pred.model){
        each.obs <- as.vector(pred.model[kk, ])
        Var.pred.3[kk] <- as.numeric(t(each.obs)%*%M.inv.3%*%each.obs)
      }
      
      minimum.3 <- min(Var.pred.3)
      maximum.3 <- max(Var.pred.3)
      average.3 <- mean(Var.pred.3)
      Matrix.V.A.3[count,] <- c(R,maximum.3,minimum.3,average.3)      
    } 
  }
  # Find lim.max
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    lim.max <- max(c(Matrix.V.A[,2], Matrix.V.A.2[,2], Matrix.V.A.3[,2]))
  }else{
    if(!is.null(design.matrix.2)){
      lim.max <- max(c(Matrix.V.A[,2], Matrix.V.A.2[,2]))
    }else{ 
      lim.max <- max(Matrix.V.A[,2])
    }
  }
  # Find lim.min
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    lim.min <- min(c(Matrix.V.A[,3], Matrix.V.A.2[,3], Matrix.V.A.3[,3]))
  }else{
    if(!is.null(design.matrix.2) ){
      lim.min <- min(c(Matrix.V.A[,3], Matrix.V.A.2[,3]))
    }else{
      lim.min <- min(Matrix.V.A[,3])
    }
  }
  
  #windows(width = 5, height = 5)
  par(mfrow=c(1,1),
      mai = c(0.75, 0.75, 0.75, 0.375),
      omi = c(0.075, 0.0375, 0.0375, 0.0375),
      mgp = c(2, 1, 0),
      xpd = FALSE)
  plot(Matrix.V.A[,1],Matrix.V.A[,2], type = "l", lwd = 2, lty = 1, 
       ylim = c(lim.min-4, lim.max + 12), col = "#2E2E2E",
       xlab = "Hypercube Radius", 
       ylab = "Scaled Variance",
       panel.first = grid())
  lines(Matrix.V.A[,1],Matrix.V.A[,3], lwd = 2, col = "#2E2E2E", lty = 6)
  lines(Matrix.V.A[,1],Matrix.V.A[,4], lwd = 2, col = "#2E2E2E", lty = 3)
  
  
  legend(x = 0, y = lim.max + 13,  legend = c("Max","Avg","Min"),lty=c(1,3,6),
         lwd=c(2,2,2),col=c(1,1,1),
         inset = 0.05, bg="transparent", bty = "n")
  
  # Design 2
  if(!is.null(design.matrix.2)){
    lines(Matrix.V.A.2[,1],Matrix.V.A.2[,2], lwd = 2, col = "#EE4000", lty = 1)
    lines(Matrix.V.A.2[,1],Matrix.V.A.2[,3], lwd = 2, col = "#EE4000", lty = 6)
    lines(Matrix.V.A.2[,1],Matrix.V.A.2[,4], lwd = 2, col = "#EE4000", lty = 3)
  }
  # Design 3
  if(!is.null(design.matrix.3)){
    lines(Matrix.V.A.3[,1],Matrix.V.A.3[,2], lwd = 2, col = "#3A5FCD", lty = 1)
    lines(Matrix.V.A.3[,1],Matrix.V.A.3[,3], lwd = 2, col = "#3A5FCD", lty = 6)
    lines(Matrix.V.A.3[,1],Matrix.V.A.3[,4], lwd = 2, col = "#3A5FCD", lty = 3)
  }
  # Put legend Design names
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    legend(x = 0.4, y = lim.max + 13,  
           legend = c(des.names[1],des.names[2],des.names[3]),lty=c(1,1,1),
           lwd=c(2,2,2),col=c("#2E2E2E","#EE4000","#3A5FCD"),
           inset = 0.05, bg="transparent", bty = "n")
  }else{
    if(!is.null(design.matrix.2)){
      legend(x = 0.4, y = lim.max + 13,  
             legend = c(des.names[1],des.names[2]),lty=c(1,1),
             lwd=c(2,2),col=c("#2E2E2E","#EE4000"),
             inset = 0.05, bg="transparent", bty = "n")
    }else{
      legend(x = 0.4, y = lim.max + 13,  
             legend = c(des.names[1]),lty=c(1),
             lwd=c(2),col=c("#2E2E2E"),
             inset = 0.05, bg="transparent", bty = "n")
      
    }}
  
  # Return Values
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    return(list(design.1 = Matrix.V.A, design.2 = Matrix.V.A.2, 
                design.3 = Matrix.V.A.3))
  }else{
    if(!is.null(design.matrix.2)){
      return(list(design.1 = Matrix.V.A, design.2 = Matrix.V.A.2))
    }else{
      return(Matrix.V.A)      
    }
    
  }
}


#########################  generate Factorial design ####################
gen.Factr <-function(n.vars, n.levels, varNames = NULL, scale = TRUE){
  if (n.vars <= 0 | n.vars%%1 != 0){stop("n.vars has to be a positive integer")}
  N.factorial <- n.levels^(n.vars)
  factorial.matrix  <- matrix(0,N.factorial,n.vars)
  point.co <- matrix(NA,n.vars,N.factorial)
  for(i in 1:n.vars){
    point.co[i,1:(N.factorial/(N.factorial/n.levels^i))] <- N.factorial/(n.levels^i)
  }  
  
  for(j in 1:n.vars){
    m <- 1
    for(k in 1:length(subset(point.co[j,], is.na(point.co[j,]) == FALSE))){
      if(n.levels%%2 == 1){
        seq.val <- seq(-(n.levels%/%2),(n.levels%/%2),1)
        factorial.matrix[m:(k*N.factorial/(n.levels^j)),j] <- rep(rep(seq.val,N.factorial)[k],subset(point.co[j,], is.na(point.co[j,]) == FALSE)[k])
        m <- 1+(k*N.factorial/(n.levels^j))
      }
      if(n.levels%%2 == 0){
        seq.val <- seq(-(n.levels-1),(n.levels-1),2)
        factorial.matrix[m:(k*N.factorial/(n.levels^j)),j] <- rep(rep(seq.val,N.factorial)[k],subset(point.co[j,], is.na(point.co[j,]) == FALSE)[k])
        m <- 1+(k*N.factorial/(n.levels^j))
      }
    }
  }         
  
  Factr.Dmatrix<- as.data.frame(factorial.matrix)
  
  if(!missing(varNames) && length(varNames == n.vars)) colnames(Factr.Dmatrix) <- varNames
  else colnames(Factr.Dmatrix) <- paste("X", 1:n.vars, sep = "")
  
  if(scale == TRUE){  
    mid.point <- (max(Factr.Dmatrix[,1])+min(Factr.Dmatrix[,1]))/2
    half.range <- (max(Factr.Dmatrix[,1])-min(Factr.Dmatrix[,1]))/2
    return((Factr.Dmatrix - mid.point)/half.range)
  } else return(Factr.Dmatrix)
}

########################################
############### Generate CCD ###########
########################################
gen.CCD <- function(n.vars, n.center, alpha, varNames){
  if (n.vars <= 0 | n.vars%%1 != 0){stop("n.vars has to be a positive integer")}
  if (n.center < 0){stop("n.center cannot be a negative value")}
  N.factorial <- 2^(n.vars)
  if(n.center > 0){center.matrix <- matrix(0,n.center,n.vars)}else{center.matrix <- NULL}  
  factorial.matrix  <- matrix(0,N.factorial,n.vars)
  point.co <- matrix(NA,n.vars,N.factorial)
  for(i in 1:n.vars){
    point.co[i,1:(N.factorial/(N.factorial/2^i))] <- N.factorial/(2^i)
  }  
  
  for(j in 1:n.vars){
    m <- 1
    for(k in 1:length(subset(point.co[j,], is.na(point.co[j,]) == FALSE))){
      factorial.matrix[m:(k*N.factorial/(2^j)),j] <- rep((-1)^k,subset(point.co[j,], is.na(point.co[j,]) == FALSE)[k])
      m <- 1+(k*N.factorial/(2^j))
    }
  }     
  
  axial.matrix <- matrix(0,2*n.vars,n.vars)     
  j <- 1
  for(i in 1:n.vars){
    axial.matrix[j,i] <- (-1)*alpha
    axial.matrix[j+1,i] <- alpha
    j <- j+2
  }
  
  CCD.Dmatrix<- rbind(factorial.matrix,axial.matrix,center.matrix )
  CCD.Dmatrix<- as.data.frame(CCD.Dmatrix)
  if(!missing(varNames) && length(varNames == n.vars)) colnames(CCD.Dmatrix) <- varNames
  else colnames(CCD.Dmatrix) <- paste("X", 1:n.vars, sep = "")
  return(CCD.Dmatrix)
}

####################################################################
####        generate Hartley's Small Composite Designs          ####
####################################################################
gen.HSCD<- function(k, alpha ="rotatable", n.center = 0){
  if(!(k >=2 && k <=7) | !(k%%1==0)){
    print("This function now provides Hartley's Small Composite Designs for k = 2 to 7 ")
    stop
  }  
  if(k == 2){
    if(alpha == "rotatable") a <- (2)^(1/4)
    else{if(alpha == "face-center") a <- 1
         else a<- alpha}
    HSCD2 <- matrix(c(
      -1,   -1,
      1,    1,
      -a,    0,
      a,    0,
      0,   -a,
      0,    a),
      byrow = TRUE, ncol = k)
    HSCD2 <- as.data.frame(rbind(HSCD2,matrix(rep(0,k*n.center),ncol=k)))
    names(HSCD2) <- paste("X",1:k,sep="")
    return(HSCD2)
  }
  
  if(k == 3){
    if(alpha == "rotatable") a <- (4)^(1/4)
    else{if(alpha == "face-center") a <- 1
         else a<- alpha}
    HSCD3 <- matrix(c(
      -1,   -1,   1,
      1,   -1,  -1,
      -1,    1,  -1,
      1,    1,   1,
      -a,    0,   0,
      a,    0,   0,
      0,   -a,   0,
      0,    a,   0,
      0,    0,  -a,
      0,    0,   a),
      byrow = TRUE, ncol = k)
    HSCD3 <- as.data.frame(rbind(HSCD3,matrix(rep(0,k*n.center),ncol=k)))
    names(HSCD3) <- paste("X",1:k,sep="")
    return(HSCD3)
  }
  
  if(k == 4){
    if(alpha == "rotatable") a <- (8)^(1/4)
    else{if(alpha == "face-center") a <- 1
         else a<- alpha}
    HSCD4 <- matrix(c(
      -1,   -1,  -1,   1,
      1,   -1,  -1,   1,
      -1,    1,  -1,  -1,
      1,    1,  -1,  -1,     
      -1,   -1,   1,  -1,
      1,   -1,   1,  -1,
      -1,    1,   1,   1,
      1,    1,   1,   1,          
      -a,    0,   0,   0,
      a,    0,   0,   0,
      0,   -a,   0,   0,
      0,    a,   0,   0,
      0,    0,  -a,   0,
      0,    0,   a,   0,
      0,    0,   0,  -a,
      0,    0,   0,   a),
      byrow = TRUE, ncol = k)
    HSCD4 <- as.data.frame(rbind(HSCD4,matrix(rep(0,k*n.center),ncol=k)))
    names(HSCD4) <- paste("X",1:k,sep="")
    return(HSCD4)
  }
  
  if(k == 5){
    if(alpha == "rotatable") a <- (16)^(1/4)
    else{if(alpha == "face-center") a <- 1
         else a<- alpha}
    HSCD5 <- matrix(c(
      -1,   -1,  -1,  -1,   1,
      1,   -1,  -1,  -1,  -1,
      -1,    1,  -1,  -1,  -1,
      1,    1,  -1,  -1,   1,       
      -1,   -1,   1,  -1,  -1,
      1,   -1,   1,  -1,   1,
      -1,    1,   1,  -1,   1,
      1,    1,   1,  -1,  -1,    
      -1,   -1,  -1,   1,  -1,
      1,   -1,  -1,   1,   1,
      -1,    1,  -1,   1,   1,
      1,    1,  -1,   1,  -1,    
      -1,   -1,   1,   1,   1,
      1,   -1,   1,   1,  -1,
      -1,    1,   1,   1,  -1,
      1,    1,   1,   1,   1,      
      -a,    0,   0,   0,   0,
      a,    0,   0,   0,   0,
      0,   -a,   0,   0,   0,
      0,    a,   0,   0,   0,
      0,    0,  -a,   0,   0,
      0,    0,   a,   0,   0,
      0,    0,   0,  -a,   0,
      0,    0,   0,   a,   0,
      0,    0,   0,   0,  -a,
      0,    0,   0,   0,   a),
      byrow = TRUE, ncol = k)
    HSCD5 <- as.data.frame(rbind(HSCD5,matrix(rep(0,k*n.center),ncol=k)))
    names(HSCD5) <- paste("X",1:k,sep="")
    return(HSCD5)
  }
  
  if(k == 6){
    if(alpha == "rotatable") a <- (16)^(1/4)
    else{if(alpha == "face-center") a <- 1
         else a<- alpha}
    HSCD6 <- matrix(c(
      -1,   -1,  -1,  -1,   1,   1,
      1,   -1,  -1,  -1,  -1,   1,
      -1,    1,  -1,  -1,  -1,   1,
      1,    1,  -1,  -1,   1,   1,          
      -1,   -1,   1,  -1,   1,  -1,
      1,   -1,   1,  -1,  -1,  -1,
      -1,    1,   1,  -1,  -1,  -1,
      1,    1,   1,  -1,   1,  -1,   
      -1,   -1,  -1,   1,   1,  -1,
      1,   -1,  -1,   1,  -1,  -1,
      -1,    1,  -1,   1,  -1,  -1,
      1,    1,  -1,   1,   1,  -1,   
      -1,   -1,   1,   1,   1,   1,
      1,   -1,   1,   1,  -1,   1,
      -1,    1,   1,   1,  -1,   1,
      1,    1,   1,   1,   1,   1,      
      -a,    0,   0,   0,   0,   0,
      a,    0,   0,   0,   0,   0,
      0,   -a,   0,   0,   0,   0,
      0,    a,   0,   0,   0,   0,
      0,    0,  -a,   0,   0,   0,
      0,    0,   a,   0,   0,   0,
      0,    0,   0,  -a,   0,   0,
      0,    0,   0,   a,   0,   0,
      0,    0,   0,   0,  -a,   0,
      0,    0,   0,   0,   a,   0,
      0,    0,   0,   0,   0,  -a,
      0,    0,   0,   0,   0,   a    
    ),byrow = TRUE, ncol = k)
    HSCD6 <- as.data.frame(rbind(HSCD6,matrix(rep(0,k*n.center),ncol=k)))
    names(HSCD6) <- paste("X",1:k,sep="")
    return(HSCD6)
  }
  
  if(k == 7){
    if(alpha == "rotatable") a <- (32)^(1/4)
    else{if(alpha == "face-center") a <- 1
         else a<- alpha}
    HSCD7 <- matrix(c(
      -1,   -1,  -1,  -1,  -1,   1,   1,
      1,   -1,  -1,  -1,  -1,  -1,   1,
      -1,    1,  -1,  -1,  -1,  -1,   1,
      1,    1,  -1,  -1,  -1,   1,   1,    
      -1,   -1,   1,  -1,  -1,   1,  -1,
      1,   -1,   1,  -1,  -1,  -1,  -1,
      -1,    1,   1,  -1,  -1,  -1,  -1,
      1,    1,   1,  -1,  -1,   1,  -1,    
      -1,   -1,  -1,   1,  -1,   1,  -1,
      1,   -1,  -1,   1,  -1,  -1,  -1,
      -1,    1,  -1,   1,  -1,  -1,  -1,
      1,    1,  -1,   1,  -1,   1,  -1,
      -1,   -1,   1,   1,  -1,   1,   1,
      1,   -1,   1,   1,  -1,  -1,   1,
      -1,    1,   1,   1,  -1,  -1,   1,
      1,    1,   1,   1,  -1,   1,   1,
      -1,   -1,  -1,  -1,   1,   1,   1,
      1,   -1,  -1,  -1,   1,  -1,   1,
      -1,    1,  -1,  -1,   1,  -1,   1,
      1,    1,  -1,  -1,   1,   1,   1,    
      -1,   -1,   1,  -1,   1,   1,  -1,
      1,   -1,   1,  -1,   1,  -1,  -1,
      -1,    1,   1,  -1,   1,  -1,  -1,
      1,    1,   1,  -1,   1,   1,  -1,  
      -1,   -1,  -1,   1,   1,   1,  -1,
      1,   -1,  -1,   1,   1,  -1,  -1,
      -1,    1,  -1,   1,   1,  -1,  -1,
      1,    1,  -1,   1,   1,   1,  -1, 
      -1,   -1,   1,   1,   1,   1,   1,
      1,   -1,   1,   1,   1,  -1,   1,
      -1,    1,   1,   1,   1,  -1,   1,
      1,    1,   1,   1,   1,   1,   1,          
      -a,    0,   0,   0,   0,   0,   0,
      a,    0,   0,   0,   0,   0,   0,
      0,   -a,   0,   0,   0,   0,   0,
      0,    a,   0,   0,   0,   0,   0,
      0,    0,  -a,   0,   0,   0,   0,
      0,    0,   a,   0,   0,   0,   0,
      0,    0,   0,  -a,   0,   0,   0,
      0,    0,   0,   a,   0,   0,   0,
      0,    0,   0,   0,  -a,   0,   0,
      0,    0,   0,   0,   a,   0,   0,
      0,    0,   0,   0,   0,  -a,   0,
      0,    0,   0,   0,   0,   a,   0,    
      0,    0,   0,   0,   0,   0,  -a,
      0,    0,   0,   0,   0,   0,   a   
    ),byrow = TRUE, ncol = k)
    HSCD7 <- as.data.frame(rbind(HSCD7,matrix(rep(0,k*n.center),ncol=k)))
    names(HSCD7) <- paste("X",1:k,sep="")
    return(HSCD7)
  }  
}

####################################################################
####              generate Hybrid Design Roquemore              ####
####################################################################
gen.Roquemore<- function(k,  n.center = 0){
  if(!(k == 3 || k == 4 || k == 6) | !(k%%1==0)){
    print("Roquemore's designs are allowed for k = 3, 4, or 6 only")
    stop
  }  
  if(k == 3){
    R310 <- matrix(c(
      0,     0,   1.2906,
      0,     0,   -.1360, 
      -1,    -1,    .6386, 
      1,    -1,    .6386, 
      -1,     1,    .6386, 
      1,     1,    .6386, 
      1.736,     0,   -.9273,
      -1.736,     0,   -.9273,
      0, 1.736,   -.9273,
      0,-1.736,   -.9273),
      byrow = TRUE, ncol = k)
    R310 <- as.data.frame(rbind(R310,matrix(rep(0,k*n.center),ncol=k)))
    names(R310) <- paste("X",1:k,sep="")
    
    R311A <- matrix(c(
      0,       0,     sqrt(2),
      0,       0,    -sqrt(2), 
      -1,      -1,   1/sqrt(2), 
      1,      -1,   1/sqrt(2),
      -1,       1,   1/sqrt(2),
      1,       1,   1/sqrt(2),
      sqrt(2),       0,  -1/sqrt(2),
      -sqrt(2),       0,  -1/sqrt(2),
      0, sqrt(2),  -1/sqrt(2),
      0,-sqrt(2),  -1/sqrt(2),
      0,       0,          0),
      byrow = TRUE, ncol = k)
    R311A <- as.data.frame(rbind(R311A,matrix(rep(0,k*n.center),ncol=k)))
    names(R311A) <- paste("X",1:k,sep="")
    
    a1.hb <- .7507;   a2.hb <- 2.1063;
    R311B <- matrix(c(
      0,       0,     sqrt(6),
      0,       0,    -sqrt(6), 
      -a1.hb,   a2.hb,           1, 
      a2.hb,   a1.hb,           1, 
      a1.hb,  -a2.hb,           1, 
      -a2.hb,  -a1.hb,           1, 
      a1.hb,   a2.hb,          -1, 
      a2.hb,  -a1.hb,          -1, 
      -a1.hb,  -a2.hb,          -1, 
      -a2.hb,   a1.hb,          -1,
      0,       0,          0), 
      byrow = TRUE, ncol = k)
    R311B <- as.data.frame(rbind(R311B,matrix(rep(0,k*n.center),ncol=k)))
    names(R311B) <- paste("X",1:k,sep="")
    
    Roquemore.3<- list(R310=R310,R311A=R311A,R311B=R311B)
    return(Roquemore.3)
  }  
  if(k == 4){
    R416A <- matrix(c(
      0,       0,       0,  1.7844,
      0,       0,       0, -1.4945, 
      -1,      -1,      -1,    .6444,
      1,      -1,      -1,    .6444,     
      -1,       1,      -1,    .6444,     
      1,       1,      -1,    .6444,       
      -1,      -1,       1,    .6444,
      1,      -1,       1,    .6444,     
      -1,       1,       1,    .6444,     
      1,       1,       1,    .6444,   
      1.6853,       0,       0,   -.9075,   
      -1.6853,       0,       0,   -.9075,    
      0,  1.6853,       0,   -.9075,   
      0, -1.6853,       0,   -.9075,    
      0,       0,  1.6853,   -.9075,   
      0,       0, -1.6853,   -.9075),  
      byrow = TRUE, ncol = k)
    R416A <- as.data.frame(rbind(R416A,matrix(rep(0,k*n.center),ncol=k)))
    names(R416A) <- paste("X",1:k,sep="")
    
    R416B <- matrix(c(
      0,       0,       0,  1.7317,
      0,       0,       0, -0.2692, 
      -1,      -1,      -1,    .6045,
      1,      -1,      -1,    .6045,     
      -1,       1,      -1,    .6045,     
      1,       1,      -1,    .6045,       
      -1,      -1,       1,    .6045,
      1,      -1,       1,    .6045,     
      -1,       1,       1,    .6045,     
      1,       1,       1,    .6045,   
      1.5177,       0,       0,  -1.0498,   
      -1.5177,       0,       0,  -1.0498,  
      0,  1.5177,       0,  -1.0498,   
      0, -1.5177,       0,  -1.0498, 
      0,       0,  1.5177,  -1.0498,   
      0,       0, -1.5177,  -1.0498),  
      byrow = TRUE, ncol = k)
    R416B <- as.data.frame(rbind(R416B,matrix(rep(0,k*n.center),ncol=k)))
    names(R416B) <- paste("X",1:k,sep="")
    
    R416C <- matrix(c(
      0,       0,       0,   1.7654,
      -1,      -1,      -1,    .5675,
      1,      -1,      -1,    .5675,     
      -1,       1,      -1,    .5675,     
      1,       1,      -1,    .5675,       
      -1,      -1,       1,    .5675,
      1,      -1,       1,    .5675,     
      -1,       1,       1,    .5675,     
      1,       1,       1,    .5675,   
      1.4697,       0,       0,  -1.0509,   
      -1.4697,       0,       0,  -1.0509,  
      0,  1.4697,       0,  -1.0509,   
      0, -1.4697,       0,  -1.0509, 
      0,       0,  1.4697,  -1.0509,   
      0,       0, -1.4697,  -1.0509,
      0,       0,       0,        0),  
      byrow = TRUE, ncol = k)
    R416C <- as.data.frame(rbind(R416C,matrix(rep(0,k*n.center),ncol=k)))
    names(R416C) <- paste("X",1:k,sep="")
    
    Roquemore.4<- list(R416A=R416A,R416B=R416B,R416C=R416C)
    return(Roquemore.4)
  }
  if(k == 6){
    
    R628A <- matrix(c(
      0,   0,   0,  0,  0,  4/sqrt(3),
      -1,  -1,  -1, -1, -1,  1/sqrt(3),
      1,   1,  -1, -1, -1,  1/sqrt(3),
      1,  -1,   1, -1, -1,  1/sqrt(3),
      -1,   1,   1, -1, -1,  1/sqrt(3),
      1,  -1,  -1,  1, -1,  1/sqrt(3),
      -1,   1,  -1,  1, -1,  1/sqrt(3),
      -1,  -1,   1,  1, -1,  1/sqrt(3),
      1,   1,   1,  1, -1,  1/sqrt(3),          
      1,  -1,  -1, -1,  1,  1/sqrt(3),
      -1,   1,  -1, -1,  1,  1/sqrt(3),
      -1,  -1,   1, -1,  1,  1/sqrt(3),
      1,   1,   1, -1,  1,  1/sqrt(3),
      -1,  -1,  -1,  1,  1,  1/sqrt(3),
      1,   1,  -1,  1,  1,  1/sqrt(3),
      1,  -1,   1,  1,  1,  1/sqrt(3),
      -1,   1,   1,  1,  1,  1/sqrt(3),  
      2,   0,   0,  0,  0, -2/sqrt(3),
      -2,   0,   0,  0,  0, -2/sqrt(3), 
      0,   2,   0,  0,  0, -2/sqrt(3),
      0,  -2,   0,  0,  0, -2/sqrt(3),  
      0,   0,   2,  0,  0, -2/sqrt(3),
      0,   0,  -2,  0,  0, -2/sqrt(3), 
      0,   0,   0,  2,  0, -2/sqrt(3),
      0,   0,   0, -2,  0, -2/sqrt(3),  
      0,   0,   0,  0,  2, -2/sqrt(3),
      0,   0,   0,  0, -2, -2/sqrt(3), 
      0,   0,   0,  0,  0,          0),
      byrow = TRUE, ncol = k)
    
    R628A <- as.data.frame(rbind(R628A,matrix(rep(0,k*n.center),ncol=k)))
    names(R628A) <- paste("X",1:k,sep="")
    
    a3 <- 2.1749
    R628B <- matrix(c(
      0,   0,   0,  0,  0,   2.3677,
      -1,  -1,  -1, -1, -1,    .6096,
      1,   1,  -1, -1, -1,    .6096,
      1,  -1,   1, -1, -1,    .6096,
      -1,   1,   1, -1, -1,    .6096,
      1,  -1,  -1,  1, -1,    .6096,
      -1,   1,  -1,  1, -1,    .6096,
      -1,  -1,   1,  1, -1,    .6096,
      1,   1,   1,  1, -1,    .6096,     
      1,  -1,  -1, -1,  1,    .6096,
      -1,   1,  -1, -1,  1,    .6096,
      -1,  -1,   1, -1,  1,    .6096,
      1,   1,   1, -1,  1,    .6096,
      -1,  -1,  -1,  1,  1,    .6096,
      1,   1,  -1,  1,  1,    .6096,
      1,  -1,   1,  1,  1,    .6096,
      -1,   1,   1,  1,  1,    .6096,  
      a3,   0,   0,  0,  0,   -1.031,
      -a3,   0,   0,  0,  0,   -1.031,
      0,  a3,   0,  0,  0,   -1.031,
      0, -a3,   0,  0,  0,   -1.031,
      0,   0,  a3,  0,  0,   -1.031,
      0,   0, -a3,  0,  0,   -1.031,
      0,   0,   0, a3,  0,   -1.031,
      0,   0,   0,-a3,  0,   -1.031,
      0,   0,   0,  0, a3,   -1.031,
      0,   0,   0,  0,-a3,   -1.031,
      0,   0,   0,  0,  0,   -1.811),
      byrow = TRUE, ncol = k)
    
    R628B <- as.data.frame(rbind(R628B,matrix(rep(0,k*n.center),ncol=k)))
    names(R628B) <- paste("X",1:k,sep="")
    
    Roquemore.6<- list(R628A=R628A,R628B=R628B)
    return(Roquemore.6)
  }
}


####################################################################
####    generate Plackett-Burman Composite Designs (PBCD)        ###
####################################################################
gen.PBCD<- function(k, alpha ="rotatable", n.center = 0){
  if(!(k >=3 && k <=7) | !(k%%1==0)){
    print("This function now provides Plackett-Burman Composite Designs for k = 3 to 7 ")
    stop
  }  
  if(k == 3){
    if(alpha == "rotatable") a <- (4)^(1/4)
    else{if(alpha == "face-center") a <- 1
         else a<- alpha}
    PBCD3 <- matrix(c(
      -1,   -1,   1,
      -1,    1,  -1,
      1,   -1,  -1,
      1,    1,   1,
      -a,    0,   0,
      a,    0,   0,
      0,   -a,   0,
      0,    a,   0,
      0,    0,  -a,
      0,    0,   a),
      byrow = TRUE, ncol = k)
    PBCD3 <- as.data.frame(rbind(PBCD3,matrix(rep(0,k*n.center),ncol=k)))
    names(PBCD3) <- paste("X",1:k,sep="")
    return(PBCD3)
  }
  
  if(k == 4){
    if(alpha == "rotatable") a <- (8)^(1/4)
    else{if(alpha == "face-center") a <- 1
         else a<- alpha}
    PBCD4 <- matrix(c(
      -1,   -1,  -1,   1,
      1,   -1,  -1,   1,
      -1,    1,  -1,  -1,
      1,    1,  -1,  -1,     
      -1,   -1,   1,  -1,
      1,   -1,   1,  -1,
      -1,    1,   1,   1,
      1,    1,   1,   1,          
      -a,    0,   0,   0,
      a,    0,   0,   0,
      0,   -a,   0,   0,
      0,    a,   0,   0,
      0,    0,  -a,   0,
      0,    0,   a,   0,
      0,    0,   0,  -a,
      0,    0,   0,   a),
      byrow = TRUE, ncol = k)
    PBCD4 <- as.data.frame(rbind(PBCD4,matrix(rep(0,k*n.center),ncol=k)))
    names(PBCD4) <- paste("X",1:k,sep="")
    return(PBCD4)
  }
  
  if(k == 5){
    if(alpha == "rotatable") a <- (12)^(1/4)
    else{if(alpha == "face-center") a <- 1
         else a<- alpha}
    PBCD5 <- matrix(c(
      1,   -1,   1,   1,   1,
      1,    1,  -1,  -1,  -1,
      -1,    1,   1,  -1,   1,
      1,   -1,   1,  -1,   1,         
      1,    1,  -1,   1,   1,
      1,    1,   1,  -1,  -1,
      -1,    1,   1,   1,  -1,
      -1,   -1,   1,   1,  -1,       
      -1,   -1,  -1,  -1,   1,
      1,   -1,  -1,   1,  -1,
      -1,    1,  -1,   1,   1,
      -1,   -1,  -1,  -1,  -1,     
      -a,    0,   0,   0,   0,
      a,    0,   0,   0,   0,
      0,   -a,   0,   0,   0,
      0,    a,   0,   0,   0,
      0,    0,  -a,   0,   0,
      0,    0,   a,   0,   0,
      0,    0,   0,  -a,   0,
      0,    0,   0,   a,   0,
      0,    0,   0,   0,  -a,
      0,    0,   0,   0,   a),
      byrow = TRUE, ncol = k)
    PBCD5 <- as.data.frame(rbind(PBCD5,matrix(rep(0,k*n.center),ncol=k)))
    names(PBCD5) <- paste("X",1:k,sep="")
    return(PBCD5)
  }
  
  if(k == 6){
    if(alpha == "rotatable") a <- (16)^(1/4)
    else{if(alpha == "face-center") a <- 1
         else a<- alpha}
    PBCD6 <- matrix(c(
      -1,   -1,  -1,  -1,   1,   1,
      1,   -1,  -1,  -1,  -1,   1,
      -1,    1,  -1,  -1,  -1,   1,
      1,    1,  -1,  -1,   1,   1,          
      -1,   -1,   1,  -1,   1,  -1,
      1,   -1,   1,  -1,  -1,  -1,
      -1,    1,   1,  -1,  -1,  -1,
      1,    1,   1,  -1,   1,  -1,   
      -1,   -1,  -1,   1,   1,  -1,
      1,   -1,  -1,   1,  -1,  -1,
      -1,    1,  -1,   1,  -1,  -1,
      1,    1,  -1,   1,   1,  -1,   
      -1,   -1,   1,   1,   1,   1,
      1,   -1,   1,   1,  -1,   1,
      -1,    1,   1,   1,  -1,   1,
      1,    1,   1,   1,   1,   1,      
      -a,    0,   0,   0,   0,   0,
      a,    0,   0,   0,   0,   0,
      0,   -a,   0,   0,   0,   0,
      0,    a,   0,   0,   0,   0,
      0,    0,  -a,   0,   0,   0,
      0,    0,   a,   0,   0,   0,
      0,    0,   0,  -a,   0,   0,
      0,    0,   0,   a,   0,   0,
      0,    0,   0,   0,  -a,   0,
      0,    0,   0,   0,   a,   0,
      0,    0,   0,   0,   0,  -a,
      0,    0,   0,   0,   0,   a    
    ),byrow = TRUE, ncol = k)
    PBCD6 <- as.data.frame(rbind(PBCD6,matrix(rep(0,k*n.center),ncol=k)))
    names(PBCD6) <- paste("X",1:k,sep="")
    return(PBCD6)  
  }
  
  if(k == 7){
    if(alpha == "rotatable") a <- (24)^(1/4)
    else{if(alpha == "face-center") a <- 1
         else a<- alpha}
    PBCD7 <- matrix(c(
      1,   -1,  -1,  -1,   1,  -1,  -1,
      1,    1,  -1,  -1,  -1,   1,   1,
      1,    1,   1,  -1,  -1,  -1,  -1,
      1,    1,   1,  -1,  -1,  -1,   1,             
      1,    1,   1,   1,  -1,  -1,  -1,
      -1,    1,   1,   1,   1,  -1,  -1,
      1,   -1,   1,   1,   1,   1,  -1,
      -1,    1,  -1,   1,   1,   1,  -1,            
      1,   -1,   1,   1,   1,   1,   1,
      1,    1,  -1,  -1,   1,   1,   1,
      -1,    1,   1,   1,  -1,   1,   1,
      -1,   -1,   1,  -1,   1,  -1,   1,      
      1,   -1,  -1,   1,  -1,   1,   1,
      1,    1,  -1,   1,   1,  -1,  -1,
      -1,    1,   1,  -1,   1,   1,   1,
      -1,   -1,   1,  -1,  -1,   1,  -1,      
      1,   -1,  -1,   1,  -1,  -1,   1,
      -1,    1,  -1,   1,   1,  -1,   1,
      1,   -1,   1,  -1,   1,   1,  -1,
      -1,    1,  -1,  -1,  -1,   1,  -1,
      -1,   -1,   1,   1,  -1,  -1,   1,
      -1,   -1,  -1,  -1,   1,  -1,   1,
      -1,   -1,  -1,   1,  -1,   1,  -1,
      -1,   -1,  -1,  -1,  -1,  -1,  -1,      
      -a,    0,   0,   0,   0,   0,   0,
      a,    0,   0,   0,   0,   0,   0,
      0,   -a,   0,   0,   0,   0,   0,
      0,    a,   0,   0,   0,   0,   0,
      0,    0,  -a,   0,   0,   0,   0,
      0,    0,   a,   0,   0,   0,   0,
      0,    0,   0,  -a,   0,   0,   0,
      0,    0,   0,   a,   0,   0,   0,
      0,    0,   0,   0,  -a,   0,   0,
      0,    0,   0,   0,   a,   0,   0,
      0,    0,   0,   0,   0,  -a,   0,
      0,    0,   0,   0,   0,   a,   0,    
      0,    0,   0,   0,   0,   0,  -a,
      0,    0,   0,   0,   0,   0,   a   
    ),byrow = TRUE, ncol = k)
    PBCD7 <- as.data.frame(rbind(PBCD7,matrix(rep(0,k*n.center),ncol=k)))
    names(PBCD7) <- paste("X",1:k,sep="")
    return(PBCD7)
  }  
}

####################################################################
####         generate (Doehlert) Uniform Shell Designs (USD)     ###
####################################################################

gen.USD<- function(k, alpha = 1){
  if(!(k >=2 && k <=6) | !(k%%1==0)){
    print("This function now provides Uniform Shell (Doehlert) Designs for k = 2 to 6")
    stop
  }  
  if(k == 2){
    X1 = c(0,1,-1,.5,-.5,.5,-.5)
    X2 = c(0,0,0,.86602,-.86602,-.86602,.86602)
    Doehlert2 <- data.frame(X1,X2)
    return(alpha*Doehlert2)
  } 
  
  if(k == 3){
    X1 = c(0,1,-1,.5,-.5,.5,-.5,.5,-.5,.5,0,-.5,0)
    X2 = c(0,0,0,.86602,-.86602,-.86602,.86602,.28868,-.28868,-.28868,.57735,.28868,-.57735)
    X3 = c(rep(0,7),.81650,-.81650,-.81650,-.81650,.81650,.81650)
    Doehlert3 <- data.frame(X1,X2,X3)
    return(alpha*Doehlert3)
  }
  
  if(k == 4){
    X1 = c(0,1,-1,.5,-.5,.5,-.5,.5,-.5,.5,0,-.5,0,.5,-.5,.5,0,0,-.5,0,0)
    X2 = c(0,0,0,.86602,-.86602,-.86602,.86602,.28868,-.28868,-.28868,.57735,.28868,-.57735,
           .28868,-.28868,-.28868,.57735,0,.28868,-.57735,0)
    X3 = c(rep(0,7),.81650,-.81650,-.81650,-.81650,.81650,.81650,.20413,-.20413,-.20413,-.20413,
           .61238,.20413,.20413,-.61238)
    X4 = c(rep(0,13),.79057,-.79057,-.79057,-.79057,-.79057,.79057,.79057,.79057)
    Doehlert4 <- data.frame(X1,X2,X3,X4)
    return(alpha*Doehlert4)
  }
  
  if(k == 5){
    X1 = c(0,1,-1,.5,-.5,.5,-.5,.5,-.5,.5,0,-.5,0,.5,-.5,.5,0,0,-.5,0,0,
           .5,-.5,.5,0,0,0,-.5,0,0,0)
    X2 = c(0,0,0,.86602,-.86602,-.86602,.86602,.28868,-.28868,-.28868,.57735,.28868,-.57735,
           .28868,-.28868,-.28868,.57735,0,.28868,-.57735,0,.28868,-.28868,-.28868,.57735,0,0,
           .28868,-.57735,0,0)
    X3 = c(rep(0,7),.81650,-.81650,-.81650,-.81650,.81650,.81650,.20413,-.20413,-.20413,-.20413,
           .61238,.20413,.20413,-.61238,.20413,-.20413,-.20413,-.20413,.61238,0,.20413,.20413,
           -.61238,0)
    X4 = c(rep(0,13),.79057,-.79057,-.79057,-.79057,-.79057,.79057,.79057,.79057,
           .15812,-.15812,-.15812,-.15812,-.15812,.63246,.15812,.15812,.15812,-.63246)
    X5 = c(rep(0,21),.77460,rep(-.77460,5),rep(.77460,4))
    Doehlert5 <- data.frame(X1,X2,X3,X4,X5)
    return(alpha*Doehlert5)
  }
  
  if(k == 6){
    X1 = c(0,1,-1,.5,-.5,.5,-.5,.5,-.5,.5,0,-.5,0,.5,-.5,.5,0,0,-.5,0,0,
           .5,-.5,.5,0,0,0,-.5,0,0,0,.5,-.5,.5,0,0,0,0,-.5,0,0,0,0)
    X2 = c(0,0,0,.86602,-.86602,-.86602,.86602,.28868,-.28868,-.28868,.57735,.28868,-.57735,
           .28868,-.28868,-.28868,.57735,0,.28868,-.57735,0,.28868,-.28868,-.28868,.57735,0,0,
           .28868,-.57735,0,0,.28868,-.28868,-.28868,.57735,0,0,0,.28868,-.57735,0,0,0)
    X3 = c(rep(0,7),.81650,-.81650,-.81650,-.81650,.81650,.81650,.20413,-.20413,-.20413,-.20413,
           .61238,.20413,.20413,-.61238,.20413,-.20413,-.20413,-.20413,.61238,0,.20413,.20413,
           -.61238,0,.20413,-.20413,-.20413,-.20413,.61238,0,0,.20413,.20413,-.61238,0,0)
    X4 = c(rep(0,13),.79057,-.79057,-.79057,-.79057,-.79057,.79057,.79057,.79057,
           .15812,-.15812,-.15812,-.15812,-.15812,.63246,.15812,.15812,.15812,-.63246,
           .15812,-.15812,-.15812,-.15812,-.15812,.63246,0,.15812,.15812,.15812,-.63246,0)
    X5 = c(rep(0,21),.77460,rep(-.77460,5),rep(.77460,4),.1291,-.1291,-.1291,-.1291,-.1291,
           -.1291,.6455,.1291,.1291,.1291,.1291,-.6455)
    X6 = c(rep(0,31),.76376,rep(-.76376,6),rep(.76376,5))
    Doehlert6 <- data.frame(X1,X2,X3,X4,X5,X6)
    return(alpha*Doehlert6)
  } 
}

####################################################################
####            generate Box-Behnken Designs (BBD)               ###
####################################################################
gen.BBD <- function(k, n.center = 1){
  if(!(k >=3 && k <=7) | !(k%%1==0)){
    print("This function now provides Box-Behnken Designs for k = 3 to 7")
    stop
  } 
  if(k == 3 | k == 4 | k == 5){
    n.vars <- k
    Two.col<- as.matrix(gen.Factr(n.vars = 2, n.levels = 2, scale = TRUE))
    N.row <- 4*choose(n.vars,2) + n.center
    BBD.matrix  <- matrix(0,N.row,n.vars)
    count <- 0
    for(i in 1:(n.vars-1)){
      for(j in (i+1):n.vars){
        count <- count + 1
        BBD.matrix[(4*(count-1)+1):(4*count),i] <- Two.col[,1]
        BBD.matrix[(4*(count-1)+1):(4*count),j] <- Two.col[,2]
      }
    }
    BBD.matrix <- as.data.frame(BBD.matrix)
    names(BBD.matrix) <- paste("X", 1:k, sep = "")
    return(BBD.matrix)
  }
  if(k == 6){
    Th3<- as.matrix(gen.Factr(n.vars = 3, n.levels = 2, scale = TRUE))
    bch1<- cbind(Th3[,1], Th3[,2], rep(0,8), Th3[,3], rep(0,8), rep(0,8))
    bch2<- cbind(rep(0,8) , Th3[,1], Th3[,2], rep(0,8),  Th3[,3],rep(0,8))
    bch3<- cbind(rep(0,8) , rep(0,8), Th3[,1], Th3[,2],  rep(0,8), Th3[,3])
    bch4<- cbind(Th3[,1] , rep(0,8),  rep(0,8), Th3[,2],  Th3[,3], rep(0,8))
    bch5<- cbind(rep(0,8) , Th3[,1],  rep(0,8), rep(0,8),  Th3[,2], Th3[,3])
    bch6<- cbind(Th3[,1] , rep(0,8), Th3[,2], rep(0,8),  rep(0,8), Th3[,3])
    BBD.matrix <- rbind(bch1, bch2, bch3, bch4, bch5, bch6, matrix(rep(0,k*n.center),ncol=k))
    BBD.matrix <- as.data.frame(BBD.matrix)
    names(BBD.matrix) <- paste("X", 1:k, sep = "")
    return(BBD.matrix)      
  }
  if(k == 7){
    Th3<- as.matrix(gen.Factr(n.vars = 3, n.levels = 2, scale = TRUE))
    bch1<- cbind(Th3[,1] , Th3[,2] , rep(0,8) , Th3[,3],  rep(0,8), rep(0,8), rep(0,8))
    bch2<- cbind(rep(0,8) , Th3[,1] , Th3[,2] , rep(0,8),  Th3[,3], rep(0,8), rep(0,8))
    bch3<- cbind(rep(0,8) , rep(0,8) , Th3[,1] , Th3[,2],  rep(0,8), Th3[,3], rep(0,8))
    bch4<- cbind(rep(0,8) , rep(0,8) , rep(0,8) , Th3[,1],  Th3[,2], rep(0,8), Th3[,3])
    bch5<- cbind(Th3[,1] , rep(0,8) , rep(0,8) , rep(0,8),  Th3[,2], Th3[,3], rep(0,8))
    bch6<- cbind(rep(0,8) , Th3[,1] , rep(0,8) , rep(0,8),  rep(0,8), Th3[,2], Th3[,3])
    bch7<- cbind(Th3[,1] , rep(0,8) , Th3[,2] , rep(0,8),  rep(0,8), rep(0,8), Th3[,3])
    BBD.matrix <- rbind(bch1, bch2, bch3, bch4, bch5, bch6, bch7, matrix(rep(0,k*n.center),ncol=k))
    BBD.matrix <- as.data.frame(BBD.matrix)
    names(BBD.matrix) <- paste("X", 1:k, sep = "")
    return(BBD.matrix)      
  }
}

####################################################################
####            generate Exact D-optimal (Borkowski)             ###
####################################################################
Borkowski2003 <- function(criterion, k, N){
  if(!(k == 2 | k == 3) | !(k%%1==0) | !(criterion == "A" | criterion == "D" | criterion == "G"| criterion == "IV")){
    print("This function now provides only exact D-, A-, G-, and IV-optimal designs for k = 2 and 3")
    stop
  } 
  if(criterion == "D"){
    if(k == 2){
      if(N == 6){
        a = .394449; b = .131483;
        D.exact.k2.N6 <- data.frame(X1=c(1,-1, -1, 1, a,-b),
                                    X2=c(1, 1,-1, -a,-1, b))
        return(D.exact.k2.N6) 
      }
      if(N == 7){
        a = 0.067476; b = 0.091516;
        D.exact.k2.N7 <- data.frame(X1=c(-1,-1, 1, 1, 1, a,-b),
                                    X2=c(-1, 1,-1, 1,-a,-1, b))
        return(D.exact.k2.N7)
      }
      if(N == 8){
        a = 0.082078; b = -0.215160;
        D.exact.k2.N8 <- data.frame(X1=c(-1,-1, 1, 1, 0, 1,-1, 0),
                                    X2=c(-1, 1,-1, 1, 1, a, a, b))
        return(D.exact.k2.N8)
      }
      if(N == 9){
        D.exact.k2.N9 <- data.frame(X1=c(-1,-1, 1, 1, 0, 0, 1,-1, 0),
                                    X2=c(-1, 1,-1, 1, 1,-1, 0, 0, 0))
        return(D.exact.k2.N9)
      }
      if(N == 10){
        a = 0.099329; b = 0.016983; c = 0.024346;
        D.exact.k2.N10 <- data.frame(X1=c(-1,-1, 1, 1, 1, 1,-a,-b,-1, c),
                                     X2=c(-1, 1,-1, 1,-1, a,-1, 1, b,-c))
        return(D.exact.k2.N10)
      }
      if(N == 11){
        a = 0.107871; b = -0.047248; 
        D.exact.k2.N11 <- data.frame(X1=c(-1,-1, 1, 1,-1, 1, 0, 0, 1,-1, 0),
                                     X2=c(-1, 1,-1, 1,-1,-1, 1,-1, a, a, b))
        return(D.exact.k2.N11)
      }
      if(N == 12){
        a = 0.006125; b =  0.101286;  c = 0.023820; 
        D.exact.k2.N12 <- data.frame(X1=c(-1,-1, 1, 1,-1, 1, 1, 1,-a,-1,-b, c),
                                     X2=c(-1, 1,-1, 1,-1, 1,-1, a,-1, b, 1,-c))
        return(D.exact.k2.N12)
      }
      if((N <= 5) | (N >= 13) | !(N%%1 == 0)){
        print("Only N = 6 to 12")
        stop
      }
    }
    if(k == 3){
      if(N == 10){
        a = .2912; b = .1925;
        D.exact.k3.N10 <- matrix(c(
          a, -1, -1, -1, a, -1, -1, -1, a,
          -b, 1, -b, -b, -b, 1, 1, -b, -b, -1,1,1,1,1,-1,1,-1,1,1,1,1
        ),byrow=TRUE, ncol = 3)
        D.exact.k3.N10 <- as.data.frame(D.exact.k3.N10)
        names(D.exact.k3.N10) <- paste("X",1:k,sep="")
        return(D.exact.k3.N10)
      }
      if(N == 11){
        part1<- gen.CCD(n.vars = 3, n.center = 0, alpha = 0)[1:8,]
        D.exact.k3.N11 <- rbind(part1,c(1,0,0),c(0,1,0),c(0,0,1))
        return(D.exact.k3.N11)
      }
      if(N == 12){
        CCD<- gen.CCD(n.vars = 3, n.center = 0, alpha = 0)
        CCD <- CCD[-(9:14),]
        a = -0.0150; b = 0.0426;
        D.exact.k3.N12 <- rbind(CCD, c(1,a,a),c(a,1,a),c(a,a,1),c(b,-1,-1))
        return(D.exact.k3.N12)
      }
      if(N == 13){
        CCD<- gen.CCD(n.vars = 3, n.center = 0, alpha = 0)
        CCD <- CCD[-(9:14),]
        CCD <- CCD[-8,]
        a = 0.1463; b = 0.0644;
        D.exact.k3.N13 <- rbind(CCD, c(a,1,1),c(1,a,1),c(1,1,a),c(b,b,-1),c(b,-1,b),c(-1,b,b))
        return(D.exact.k3.N13)
      }
      if(N == 14){
        CCD<- gen.CCD(n.vars = 3, n.center = 0, alpha = 0)
        CCD <- CCD[-(9:14),]
        a = -0.0187; b = 0.0583;
        D.exact.k3.N14 <- rbind(CCD, c(a,1,1),c(1,a,1),c(1,1,a),c(b,b,-1),c(b,-1,b),c(-1,b,b))
        return(D.exact.k3.N14)
      }
      if(N == 15){
        CCD<- gen.CCD(n.vars = 3, n.center = 0, alpha = 0)
        CCD <- CCD[-(9:14),]
        D.exact.k3.N15 <- rbind(CCD, c(-.076,-1,1),c(-.0273,-1,-1),c(-1,-0.0905,-0.0282),
                                c(1,.0305,-1),c(1,-1,.0193),c(.073,1,-.0417),c(.1205,.0431,1))
        return(D.exact.k3.N15)
      }
      if(N == 16){
        D.exact.k3.N16 <- matrix(c(
          .0513,1,-1,1,1,-.0513,.029,1,1,-1,1,-.029,-.0819,-.0258,-1,1,-.0258,-1,1,-.0258,.0819,-1,-.0426,1
        ), byrow = TRUE, ncol = 3)
        D.exact.k3.N16 <- as.data.frame(D.exact.k3.N16)
        names(D.exact.k3.N16) <- paste("X",1:k,sep="")
        D.exact.k3.N16 <- rbind(D.exact.k3.N16, gen.CCD(3,0,1)[1:8,])      
        return(D.exact.k3.N16)
      }
      if((N <= 9) | (N >= 17) | !(N%%1 == 0)){
        print("Only N = 10 to 16")
        stop
      }
    }
  }
  if(criterion == "A"){
    if(k == 2){
      if(N == 6){
        a = .503816; b = -.220484;
        A.exact.k2.N6 <- data.frame(X1=c(1,-1, 0, -1, 1,0),
                                    X2=c(-1, -1,1, a,a, b))
        return(A.exact.k2.N6) 
      }
      if(N == 7){
        a = .478014; b = -.218806;
        A.exact.k2.N7 <- data.frame(X1=c(-1,1,0,-1,1,0,0),
                                    X2=c(-1,-1,1,a,a,b,b))
        return(A.exact.k2.N7)
      }
      if(N == 8){
        a = .078929; b = .571386; c = -.029071; d = -.179098; e = .100290;
        A.exact.k2.N8 <- data.frame(X1=c(-1,-1,1,a,1,c,c,-1),
                                    X2=c(-1,1,-1,1,b,d,d,e))
        return(A.exact.k2.N8)
      }
      if(N == 9){
        A.exact.k2.N9 <- data.frame(X1=c(-1,-1, 1, 1, 0, 0, 1,-1, 0),
                                    X2=c(-1, 1,-1, 1, 1,-1, 0, 0, 0))
        return(A.exact.k2.N9)
      }
      if(N == 10){
        A.exact.k2.N10 <- data.frame(X1=c(-1,-1, 1, 1,0,0,1,-1,0,0),
                                     X2=c(-1, 1,-1, 1,1,-1,0,0,0,0))
        A.exact.k2.N10<- list(Design = A.exact.k2.N10, Note = "This is also an exact IV-optimal design")
        return(A.exact.k2.N10)
      }
      if(N == 11){
        A.exact.k2.N11 <- data.frame(X1=c(-1,-1, 1, 1,0,0,1,-1,0,0,0),
                                     X2=c(-1, 1,-1, 1,1,-1,0,0,0,0,0))
        A.exact.k2.N11<- list(Design = A.exact.k2.N11, Note = "This is also an exact IV-optimal design")        
        return(A.exact.k2.N11)
      }
      if(N == 12){
        a = .059358; b =  -.021281;
        A.exact.k2.N12 <- data.frame(X1=c(-1,-1, 1, 1,1,1,-1,a,a,b,b,b),
                                     X2=c(-1, 1,-1, 1,0,0,0,-1,1,0,0,0))
        return(A.exact.k2.N12)
      }
      if((N <= 5) | (N >= 13) | !(N%%1 == 0) ){
        print("Only N = 6 to 12")
        stop
      }
    }
    if(k == 3){
      if(N == 10){
        a = .1749; b = .1072;
        A.exact.k3.N10 <- matrix(c(
          -a, -1,  1, 1, a, 1, 1, -1, -a,
          -1, -b, b, b, 1, b, b, -b, -1, -1,1,-1,-1,1,1,1,1,-1,-1,-1,-1
        ),byrow=TRUE, ncol = 3)
        A.exact.k3.N10 <- as.data.frame(A.exact.k3.N10)
        names(A.exact.k3.N10) <- paste("X",1:k,sep="")
        return(A.exact.k3.N10)
      }
      if(N == 11){
        a = .1223; b = .2166; c = .0914;
        A.exact.k3.N11 <- matrix(c(
          a,  -1,  1,
          -1,   a,  1,
          -1,  -1, -a,
          -b,  -b, -1,
          -b,   1,  b,
          1,  -b,  b,
          c,   c, -c,
          1,  -1, -1,
          1,   1, -1,
          -1,   1, -1,
          1,   1,  1
        ), byrow= TRUE, ncol = 3)
        A.exact.k3.N11 <- as.data.frame(A.exact.k3.N11)
        names(A.exact.k3.N11) <- paste("X",1:k,sep="")
        return(A.exact.k3.N11)
      }
      if(N == 12){
        A.exact.k3.N12 <- matrix(c(
          .1185,      1,      1,
          -1, -.1185,      1,
          -1,      1,  .0814,
          .0232, -.0232, -.1403,
          -.1811,     -1,  .2417,
          1,  .1811,  .2417,
          -.0474,  .0474,     -1,
          -1,      1,     -1,
          -1,     -1,     -1,
          1,      1,     -1,
          1,     -1,     -1,
          1,     -1,      1
        ), byrow= TRUE, ncol = 3)
        A.exact.k3.N12 <- as.data.frame(A.exact.k3.N12)
        names(A.exact.k3.N12) <- paste("X",1:k,sep="")
        return(A.exact.k3.N12)
      }
      if(N == 13){
        A.exact.k3.N13 <- matrix(c(
          0, -.0537,      0,
          0, -.0537,      0,
          -.1243,      1,      1,
          .1243,      1,     -1,
          1, -.6353,     -1,
          -1, -.6353,      1,
          1,  .2472,      1,
          -1,  .2472,     -1,
          1,      1, -.1243,
          -1,      1,  .1243, 
          0,     -1,      0,
          -1,     -1,     -1,
          1,     -1,      1
        ), byrow= TRUE, ncol = 3)
        A.exact.k3.N13 <- as.data.frame(A.exact.k3.N13)
        names(A.exact.k3.N13) <- paste("X",1:k,sep="")
        return(A.exact.k3.N13)
      }
      if(N == 14){
        A.exact.k3.N14<- gen.CCD(n.vars = 3, n.center = 0, alpha = 1)
        A.exact.k3.N14 = list(Design = A.exact.k3.N14, note = "This is also an exact G-optimal design")
        return(A.exact.k3.N14)
      }
      if(N == 15){
        A.exact.k3.N15<- gen.CCD(n.vars = 3, n.center = 1, alpha = 1)
        A.exact.k3.N15 = list(Design = A.exact.k3.N15, note = "This is also an exact IV-optimal design")
        return(A.exact.k3.N15)
      }
      if(N == 16){
        A.exact.k3.N16<- gen.CCD(n.vars = 3, n.center = 0, alpha = 1)
        A.exact.k3.N16[13,] <- c(0,0,-1); A.exact.k3.N16[14,] <- c(0,0,-1); 
        A.exact.k3.N16<- rbind(A.exact.k3.N16,c(0,-1,0),c(1,0,0))
        return(A.exact.k3.N16)
      }
      if((N <= 9)  | (N >= 17) | !(N%%1 == 0) ){
        print("Only N = 10 to 16")
      }
    }
  }
  if(criterion == "G"){
    if(k == 2){
      if(N == 6){
        a = -.193256; b = .522; c = .864001
        G.exact.k2.N6 <- matrix(c(-1,-1,a,a,1,b,b,1,c,-1,-1,c), byrow = TRUE, ncol=2)
        G.exact.k2.N6 <- as.data.frame(G.exact.k2.N6)
        names(G.exact.k2.N6) <- paste("X",1:k,sep="")
        return(G.exact.k2.N6)
      }
      if(N == 7){
        a = .937817; b = -.603132; c = .773037; d = -.050674;
        G.exact.k2.N7 <- matrix(c(1,1,a,-1,-1,a,b,-1,-1,b,c,d,d,c), byrow = TRUE, ncol=2)
        G.exact.k2.N7 <- as.data.frame(G.exact.k2.N7)
        names(G.exact.k2.N7) <- paste("X",1:k,sep="")
        return(G.exact.k2.N7)
      }
      if(N == 8){
        a = .052190; b = .063291; c = .824605;
        G.exact.k2.N8 <- matrix(c(1,1,1,-1,-1,1,-1,-1,1,-a,-1,a,-b,-c,b,c), byrow = TRUE, ncol=2)
        G.exact.k2.N8 <- as.data.frame(G.exact.k2.N8)
        names(G.exact.k2.N8) <- paste("X",1:k,sep="")
        return(G.exact.k2.N8)
      }
      if(N == 9){
        a = .420229
        G.exact.k2.N9 <- matrix(c(1,1,1,-1,-1,1,-1,-1,0,0,a,1,-a,-1,1,-a,-1,a), byrow = TRUE, ncol=2)
        G.exact.k2.N9 <- as.data.frame(G.exact.k2.N9)
        names(G.exact.k2.N9) <- paste("X",1:k,sep="")
        return(G.exact.k2.N9)
      }
      if(N == 10){
        a = -.430; b = .564028; c=.176604;
        G.exact.k2.N10 <- matrix(c(1,1,1,-1,-1,1,-1,-1,0,-1,1,a,-1,a,b,1,-b,1,0,c), byrow = TRUE, ncol=2)
        G.exact.k2.N10 <- as.data.frame(G.exact.k2.N10)
        names(G.exact.k2.N10) <- paste("X",1:k,sep="")
        return(G.exact.k2.N10)
      }
      if(N == 11){
        a = -.561111; b = .158889; c = .819313; d = .000843;
        G.exact.k2.N11 <- matrix(c(1,1,1,-1,-1,1,-1,-1,-1,-1,1,a,a,1,b,-1,-1,b,c,c,d,d), byrow = TRUE, ncol=2)
        G.exact.k2.N11 <- as.data.frame(G.exact.k2.N11)
        names(G.exact.k2.N11) <- paste("X",1:k,sep="")
        return(G.exact.k2.N11)
      }
      if(N == 12){
        a = -.681289; b = .061178; c = .519483; d = -.035116; e = .895400;
        G.exact.k2.N12 <- matrix(c(1,-1,-1,1,1,a,a,1,b,-1,-1,b,c,d,d,c,1,e,e,1,-1,-e,-e,-1), byrow = TRUE, ncol=2)
        G.exact.k2.N12 <- as.data.frame(G.exact.k2.N12)
        names(G.exact.k2.N12) <- paste("X",1:k,sep="")
        return(G.exact.k2.N12)
      }
      if((N <= 5) | (N >= 13) | !(N%%1 == 0) ){
        print("Only N = 6 to 12")
        stop
      }
    }
    if(k == 3){
      if(N == 10){
        G.exact.k3.N10 <- matrix(c(
          -.9998, .4518, .8786,
          1 ,     1, .3059,
          1 ,    -1,-.9110,
          -.8810,    -1, .9900,
          -1, -.7809, -1  ,
          .0615,    1  ,   1 ,
          1 , -.6490,   1 ,
          -.9558,  1    , -.7944,
          .7592, .6918, -1,
          -.0413, -.5508, -.0283
        ),byrow = TRUE, ncol=3)
        G.exact.k3.N10 <- as.data.frame(G.exact.k3.N10)
        names(G.exact.k3.N10) <- paste("X",1:k,sep="")
        return(G.exact.k3.N10)
      }
      if(N == 11){
        G.exact.k3.N11 <- matrix(c(
          -.1140, -1, 1,-1,-.1140,1,-1,-1,.1141,-.7550,-.7550,-1,-.7550,1,.7550,1,-.7550,.7550,.1858,.1858
          ,-.1858,1,-1,-1,1,1,-1,-1,1,-1,1,1,1
        ),byrow = TRUE, ncol=3)
        G.exact.k3.N11 <- as.data.frame(G.exact.k3.N11)
        names(G.exact.k3.N11) <- paste("X",1:k,sep="")
        return(G.exact.k3.N11)
      }
      if(N == 12){
        G.exact.k3.N12 <- matrix(c(
          1,-.7188,1,1,-1,-.9960,-1,-1,-.9900,.9930,.9841,.8816,-1,1,-.8898,
          -1,-.5815,1,-.4144,-1,.9600,-.7357,1,.9983,.9998,1,-.8238,-1,.0100,
          .0940,.3524,-1,.0469,.0928,.0913,-1
        ),byrow = TRUE, ncol=3)
        G.exact.k3.N12 <- as.data.frame(G.exact.k3.N12)
        names(G.exact.k3.N12) <- paste("X",1:k,sep="")
        return(G.exact.k3.N12)
      }
      if(N == 13){
        G.exact.k3.N13 <- matrix(c(
          .3938,-1,-1,1,-.3938,-1,1,-1,-.3938,-1,0,0,0,1,0,0,0,1,1,-1,1,1,1,-1,1,1,1,-1,-1,-1,-1,-1,1,
          -1,1,-1,-1,1,1
        ),byrow = TRUE, ncol=3)
        G.exact.k3.N13 <- as.data.frame(G.exact.k3.N13)
        names(G.exact.k3.N13) <- paste("X",1:k,sep="")
        return(G.exact.k3.N13)
      }
      if(N == 14){
        G.exact.k3.N14 <-  gen.CCD(3, n.center=0, alpha=1)
        G.exact.k3.N14 <- as.data.frame(G.exact.k3.N14)
        return(G.exact.k3.N14)
      }
      if(N == 15){
        G.exact.k3.N15 <-  gen.CCD(3, n.center=0, alpha=1)
        G.exact.k3.N15 <- as.data.frame(rbind(G.exact.k3.N15,c(1,-1,-1)))
        return(G.exact.k3.N15)
      }
      if(N == 16){
        G.exact.k3.N16 <-  gen.CCD(3, n.center=0, alpha=1)[-c(13,14),]
        G.exact.k3.N16 <- rbind(G.exact.k3.N16,c(0,0,-1),c(0,0,-1),c(0,1,0),c(1,0,0))
        return(G.exact.k3.N16)
      }
      if((N <= 9)  | (N >= 17) | !(N%%1 == 0) ){
        print("Only N = 10 to 16")
      }
    }
  }
  if(criterion == "IV"){
    if(k == 2){
      if(N == 6){
        a = .707479; b = .276367; c = .144868;
        IV.exact.k2.N6 <- matrix(c(-1,1,a,1,-1,-a,1,-b,b,-1,-c,c), byrow = TRUE, ncol=2)
        IV.exact.k2.N6 <- as.data.frame(IV.exact.k2.N6)
        names(IV.exact.k2.N6) <- paste("X",1:k,sep="")
        return(IV.exact.k2.N6)
      }
      if(N == 7){
        a = .310497; b = .796856; c = .147869;
        IV.exact.k2.N7 <- matrix(c(1,-1,-1,a,-a,1,1,b,-b,-1,c,-c,c,-c), byrow = TRUE, ncol=2)
        IV.exact.k2.N7 <- as.data.frame(IV.exact.k2.N7)
        names(IV.exact.k2.N7) <- paste("X",1:k,sep="")
        return(IV.exact.k2.N7)
      }
      if(N == 8){
        a = .003859; b = .768301;  c = .094936;
        IV.exact.k2.N8 <- matrix(c(1,-1,-1,1,-1,-a,a,1,1,b,-b,-1,c,-c,c,-c), byrow = TRUE, ncol=2)
        IV.exact.k2.N8 <- as.data.frame(IV.exact.k2.N8)
        names(IV.exact.k2.N8) <- paste("X",1:k,sep="")
        return(IV.exact.k2.N8)
      }
      if(N == 9){
        a = .044158; b = .836727; c = .044687;
        IV.exact.k2.N9 <- matrix(c(-1,1,-1,-1,-1,0,a,1,a,-1,1,b,1,-b,c,0,c,0), byrow = TRUE, ncol=2)
        IV.exact.k2.N9 <- as.data.frame(IV.exact.k2.N9)
        names(IV.exact.k2.N9) <- paste("X",1:k,sep="")
        return(IV.exact.k2.N9)
      }
      if(N == 10){
        IV.exact.k2.N10 <- data.frame(X1=c(-1,-1, 1, 1,0,0,1,-1,0,0),
                                      X2=c(-1, 1,-1, 1,1,-1,0,0,0,0))
        return(IV.exact.k2.N10)
      }
      if(N == 11){
        IV.exact.k2.N11 <- data.frame(X1=c(-1,-1, 1, 1,0,0,1,-1,0,0,0),
                                      X2=c(-1, 1,-1, 1,1,-1,0,0,0,0,0))
        return(IV.exact.k2.N11)
      }
      if(N == 12){
        IV.exact.k2.N12 <- matrix(c(-1,-1,-1,1,1,-1,1,1,0,1,0,-1,1,0,-1,0,0,0,0,0,0,0,0,0), byrow = TRUE, ncol=2)
        IV.exact.k2.N12 <- as.data.frame(IV.exact.k2.N12)
        names(IV.exact.k2.N12) <- paste("X",1:k,sep="")
        return(IV.exact.k2.N12)
      }
      if((N <= 5) | (N >= 13) | !(N%%1 == 0) ){
        print("Only N = 6 to 12")
        stop
      }
    }
    if(k == 3){
      if(N == 10){
        IV.exact.k3.N10 <- matrix(c(
          -.9605,-.1025,-.1025,.1025,.9605,-.1025,.1025,-.1025,.9605,-.2553,-1,-1,1,.2553,-1
          ,1,-1,.2553,-1,1,1,1,1,1,-1,1,-1,-1,-1,1
        ),byrow = TRUE, ncol=3)
        IV.exact.k3.N10 <- as.data.frame(IV.exact.k3.N10)
        names(IV.exact.k3.N10) <- paste("X",1:k,sep="")
        return(IV.exact.k3.N10)
      }
      if(N == 11){
        IV.exact.k3.N11 <- matrix(c(
          -.0589,.0589,-.0589,-1,-.1905,.1905,.1905,1,.1905,.1905,-.1905,-1,-.2349,-1,1,1,.2349,1,1,-1,-.2349,
          -1,-1,-1,-1,1,-1,-1,1,1,1,1,-1
        ),byrow = TRUE, ncol=3)
        IV.exact.k3.N11 <- as.data.frame(IV.exact.k3.N11)
        names(IV.exact.k3.N11) <- paste("X",1:k,sep="")
        return(IV.exact.k3.N11)
      }
      if(N == 12){
        a = .0925; b = .3120; c = .1685;
        IV.exact.k3.N12 <- matrix(c(
          -a,-a,a,-a,-a,a,-1,b,-b,b,-1,-b,b,b,1,-c,1,-1,1,-c,-1,1,1,c,-1,-1,1,1,-1,1,-1,-1,-1,-1,1,1
        ),byrow = TRUE, ncol=3)
        IV.exact.k3.N12 <- as.data.frame(IV.exact.k3.N12)
        names(IV.exact.k3.N12) <- paste("X",1:k,sep="")
        return(IV.exact.k3.N12)
      }
      if(N == 13){
        a = .0537; b = .1243; c = .6353; d = .2472;
        IV.exact.k3.N13 <- matrix(c(
          0,-a,0,0,-a,0,-b,1,1,b,1,-1,1,-c,-1,-1,-c,1,1,d,1,-1,d,-1,1,1,-b,-1,1,b,0,-1,0,-1,-1,-1,1,-1,1
        ),byrow = TRUE, ncol=3)
        IV.exact.k3.N13 <- as.data.frame(IV.exact.k3.N13)
        names(IV.exact.k3.N13) <- paste("X",1:k,sep="")
        return(IV.exact.k3.N13)
      }
      if(N == 14){
        IV.exact.k3.N14 <- matrix(c(
          -.8052,1,1,-.5653,-1,-1,-1,-.1528,1,-1,-1,.0475,.3061,-1,1,1,-1,-.7339,1,.8613,1,-1,.9168,-1,1,1,-1,1,-.1531,.1993,-.0174,1,.0152
          ,.0524,-.0482,-1,-.0240,.0628,-.0056,-.0240,.0628,-.0056       
        ),byrow = TRUE, ncol=3)
        IV.exact.k3.N14 <- as.data.frame(IV.exact.k3.N14)
        names(IV.exact.k3.N14) <- paste("X",1:k,sep="")
        return(IV.exact.k3.N14)
      }
      if(N == 15){
        IV.exact.k3.N15 <- gen.CCD(3,1,1)[-c(13,14),]
        IV.exact.k3.N15 <- rbind(IV.exact.k3.N15,c(0,0,-1),c(0,0,-1))
        names(IV.exact.k3.N15) <- paste("X",1:k,sep="")
        return(IV.exact.k3.N15)
      }
      if(N == 16){
        IV.exact.k3.N16 <- gen.CCD(3,2,1)[-c(13,14),]
        IV.exact.k3.N16 <- rbind(IV.exact.k3.N16,c(0,0,-1),c(0,0,-1))
        names(IV.exact.k3.N16) <- paste("X",1:k,sep="")
        return(IV.exact.k3.N16)
      }
      if((N <= 9)  | (N >= 17) | !(N%%1 == 0) ){
        print("Only N = 10 to 16")
      }
    }
  }
}

################# Contour plots ####################

spvcontour <- function(design.matrix, shape, max.radius = sqrt(2), length = 100, nlevels = 10,
                        title = "Contour of SPVs"){
  norm2 <- function(x){return(sqrt(sum(x^2)))}
  shuffle.fun <- function(row.vec){
    num.var <- length(row.vec)
    row.vec[shuffle(num.var)]
  }
  n.var <- ncol(design.matrix)
  if(!(n.var == 2)){
    return(print("The contour plot is only 2-factor design"))
  } 
  if( (shape == "circle") && !is.na(max.radius) ){
    design.matrix<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix))),byrow=FALSE,ncol=n.var))  
    N.generated <- 150*n.var*(2^(n.var))
    model.X <- model.matrix( ~quad(.) , design.matrix)
    nrow <- nrow(model.X)
    ncol <- ncol(model.X)
    M.inv<- nrow*solve(t(model.X)%*%model.X)
    
    N.of.pts.gen <- 10*N.generated
    set.seed(1972264)
    
    Urand<- runif(n=(10*n.var*N.generated),min=-max.radius,max=max.radius)
    Urand <- ifelse(Urand < (-max.radius+.01), -max.radius, Urand)
    Urand <- ifelse(Urand > ( max.radius-.01),  max.radius, Urand)
    Urand <- ifelse(-.01 < Urand & Urand < .01 ,  0, Urand)  
    
    rand.sphere.2times <- matrix(Urand ,byrow= TRUE, ncol=n.var)
    norm2.rand.2times <- apply(rand.sphere.2times, MARGIN = 1, norm2)
    norm2.rand.rep.2times <- matrix(rep(norm2.rand.2times, each = n.var),byrow= TRUE, ncol=n.var) 
    
    Urand.2 <- runif(n = N.of.pts.gen, min = 0, max = max.radius)
    Urand.2 <- ifelse(Urand.2 < .001, 0, Urand.2)
    Urand.2 <- ifelse(Urand.2 > ( max.radius-.001),  max.radius, Urand.2)
    radius.random  <- Urand.2
    
    Surface.pts.random <- radius.random*as.data.frame(rand.sphere.2times/norm2.rand.rep.2times)
    pred.model.random  <- model.matrix( ~quad(.), Surface.pts.random)
    
    Var.pred.random<- numeric(N.of.pts.gen)
    for(kk in 1:N.of.pts.gen){
      each.obs <- as.vector(pred.model.random[kk, ])
      Var.pred.random[kk] <- as.numeric(t(each.obs)%*%M.inv%*%each.obs)
    }
    par(mfrow=c(1,1),
        mai = c(0.75, 0.75, 0.75, 0.375),
        omi = c(0.075, 0.0375, 0.0375, 0.0375),
        mgp = c(2, 1, 0),
        xpd = FALSE)
    contour(akima::interp(Surface.pts.random$V1, Surface.pts.random$V2, Var.pred.random, duplicate = "mean",
                   xo=seq(min(Surface.pts.random$V1), max(Surface.pts.random$V1), length = length),
                   yo=seq(min(Surface.pts.random$V2), max(Surface.pts.random$V2), length = length),
                   linear = TRUE, extrap=FALSE), nlevels = nlevels, labcex = 1, col="black",main = title)
    r=max.radius
    nseg=360
    x.cent <- 0
    y.cent <- 0
    
    xx <- x.cent + r*cos( seq(0,2*pi, length.out=nseg) )
    yy <- y.cent + r*sin( seq(0,2*pi, length.out=nseg) )   
    lines(xx,yy, col='red')
    points(design.matrix, pch = 19, cex = 1.5)
  }else{
    if( (shape == "square")){
      N.generated <- 15000
      model.X <- model.matrix( ~quad(.) , design.matrix)
      nrow <- nrow(model.X)
      ncol <- ncol(model.X)
      M.inv<- nrow*solve(t(model.X)%*%model.X)
      
      Urand <- runif(n = N.generated*n.var, min = -1, max = 1)
      Urand <- ifelse(Urand < (-1+.01), -1, Urand)
      Urand <- ifelse(Urand > ( 1-.01),  1, Urand)
      Urand <- ifelse(Urand > -.01 & Urand < .01,  0, Urand)
      
      Cuboidal.pts.random  <- matrix(Urand, byrow = TRUE, ncol = n.var)      
      Cuboidal.pts.random  <- as.data.frame(Cuboidal.pts.random)
      pred.model.random  <- model.matrix( ~quad(.), Cuboidal.pts.random)
      radius.random<- apply(Cuboidal.pts.random, MARGIN = 1, norm2)
      
      size.pred.rand <- dim(pred.model.random)[1]
      Var.pred.random <- numeric(size.pred.rand)
      for(kk in 1:size.pred.rand){
        each.obs <- as.vector(pred.model.random[kk, ])
        Var.pred.random[kk] <- as.numeric(t(each.obs)%*%M.inv%*%each.obs)
      }
      par(mfrow=c(1,1),
          mai = c(0.75, 0.75, 0.75, 0.375),
          omi = c(0.075, 0.0375, 0.0375, 0.0375),
          mgp = c(2, 1, 0),
          xpd = FALSE)
      contour(akima::interp(Cuboidal.pts.random$V1, Cuboidal.pts.random$V2, Var.pred.random, duplicate = "mean",
                     xo=seq(min(Cuboidal.pts.random$V1), max(Cuboidal.pts.random$V1), length = length),
                     yo=seq(min(Cuboidal.pts.random$V2), max(Cuboidal.pts.random$V2), length = length),
                     linear = TRUE, extrap=FALSE), nlevels = nlevels, labcex = 1,main = title)
      points(design.matrix, pch = 19, cex = 1.5)
      
      segments(-1, -1, -1, 1, col = "red")
      segments(-1,  1,  1, 1, col = "red")
      segments( 1,  1,  1, -1, col = "red")
      segments( 1, -1, -1, -1, col = "red")
      
    }else{ return(print("Arguments are not valid"))}
  }  
}
