ml.est <- function (y, x=NULL, model = "LN", lambda=3,  w=0.05, lambda.fix=FALSE, w.fix=FALSE,
                     eps=1e-7, max.iter=500, t.outl=0.5, graph=FALSE)
{
#------------------------------------------------------------------------------
#           Individuazione degli outlier basata su un modello mistura di 2 gaussiane
#------------------------------------------------------------------------------
#         PARAMETRI 
#  y  = matrice ( o data.frame) -  Variabili dipendenti (con possibili errori)
#  x  = matrice ( o data.frame) - Variabili indipendenti (dati esatti. P.e. da archivio amministrativo)
#  model = Indica se i dati osservati hanno distribuzione log-normale (LN) o normale (N).
#  w = proporzione dei dati contaminati (peso a priori) 
#  max.iter = numero massimo di iterazioni per la convergenza EM
#  eps = soglia di accettazione
#  lambda = fattore di inflazione della varianza
#  graph = visualizzazione dei grafici durante l'elaborazione
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Copio i dati di input su aree di appoggio
#------------------------------------------------------------------------------
  memo.y <- y <-as.matrix(y)         
  memo.x <- x
#------------------------------------------------------------------------------
#        CONTROLLI SUI PARAMETRI 
#  Eliminazione dei record contenenti missing per la stima dei parametri
#------------------------------------------------------------------------------  
  ind.NA<- which(rowSums(is.na(y)) >0 )

  if (length(ind.NA) > 0 )    {
      warning(paste("Input matrix y contains", length(ind.NA), " (%",length(ind.NA)*100/nrow(y) ,
      ") rows with missing values not included in parameter's estimation\n" ))
      y <- y[-ind.NA,,drop=FALSE]
      if (length(x) > 0) {
	      x <-as.matrix(x) 
          x <- x[-ind.NA,,drop=FALSE]
	  }	  
   }
#------------------------------------------------------------------------------
#        CONTROLLI SUI PARAMETRI 
#------------------------------------------------------------------------------
  vars <- check.vars(y,x,model,parent="ml.est")

  if (vars$ret == -9) {
     stop(vars$msg.err) 
  }
  if (vars$ret != 0) {
     warning(vars$msg.err) 
  } 
  y <- as.matrix(vars$y)
  x <- as.matrix(vars$x) 


  
  #y <- as.matrix(y, dimnames = NULL)
   
  p <- ncol(y)
  n <- nrow(y)
  omega <- rep(1,n)  
   
#  covar <- x
     
#  x <- as.matrix(cbind (rep(1,n),x))     
  q <- ncol(x) 
#------------------------------------------------------------------------------  
 
 
#---------------------       DEFINIZIONE VARIABILI       ---------------------
  
  conv <- TRUE                                                                                     
  lik0 <- 10
  oldlik <- 0
  iter <- 0
  sing <- FALSE
  if (ncol(x)+ ncol(y) < 3)   # INSERIRE BOXPLOT
      graph=FALSE 
  if (graph )   {
   lambda_all <- lambda

   if (ncol(y) >= 2)  { 
      lab  <- colnames(y)[1:2]
      Var <- y[,1:2]  
   }   
   else if (ncol(y) == 1 & ncol(x) > 1) {
      lab  <- c(colnames(x)[2],colnames(y)[1])
      Var <- cbind(x[,2],y[,1])  
   }     
   
 #  windows()
   par(mfrow=c(2,1))
  }
  
#  B <- solve(t(x) %*% x) %*% t(x) %*% y      # B ha q (ncol(x)) righe e p (ncol(y)) variabili
  B <- solve((t(x) %*% x) + (10e-8* diag(rep(1,q)))) %*% t(x) %*% y      # B ha q (ncol(x)) righe e p (ncol(y)) variabili
  B0 <- B
 sigma <- (t(y - x%*%B) %*% (y - x%*%B)) / (n-1)
  sigma <- sigma + (10e-8* diag(rep(1,p)))
  sigma0 <- sigma
  sigma2 <- (1 + lambda) * sigma  
  w1 <- 1-w
#************************  INIZIO CICLO EM  ************************************
  while (iter < max.iter & conv == TRUE)
  {
     iter <- iter + 1
#***********************    E - STEP         ************************************      
    tau1 <- post.prob(y, x, B, sigma, w1, lambda)   
    tau2 <- 1 - tau1 
#***********************    M - STEP         ************************************    

#***********************    calcolo dei pesi  **********************************
     if (!w.fix)
         w1 <- sum(tau1)/n;

#***********************        omega          **********************************
  
     omega <- as.vector(tau1 + tau2 / (1+lambda))
    
#***********************        B         **********************************
     appo <- t(x) %*% (omega * x)
     appo <- solve(appo)
     B <- appo %*% t(x) %*% (omega * y)
     gc()
#***********************        sigma         **********************************

     dif <- y - x%*%B
     
     sigma <- (t(dif) %*% (omega * dif)) / n
     if (det(sigma)  <   10e-10)  {
         sing <- TRUE
         lambda.fix <- TRUE
         stop("Estimation has been stopped at current values of parameters because 
                 the determinant of covariance matrix is less than 10e-10")
      }
     gc()
#***********************        lambda        **********************************

#     q1 <- matrix(diag(dif %*% solve(sigma1) %*% t(dif)),n,1)        ##   DIM n,1
#     q2 <- matrix(diag(dif %*% solve(sigma2) %*% t(dif)),n,1)        ##   DIM n,1
     if (!lambda.fix)  {     
        s1 <- solve(sigma)
        appo <-  t(dif) %*% (as.vector(tau2) * dif) %*% s1
        lambda <- sum(diag(as.matrix(appo))) / (p * sum(tau2)) -1
     }
     gc()
     if (graph)   {     
        plot(Var,  col = "lightgrey", main= "EM IN ACTION...\n Identifiyng outliers",  xlab=lab[1], ylab=lab[2] )
        points(Var[tau2 > 0.5, ],pch=21,col="blue",bg=paste("cyan",sample(1:4,1),sep="")) 
        
        lambda_all <- c (lambda_all, lambda)
        plot( lambda_all, xlab="n. iterations", ylab="lambda")
     }
     
#***********************    CONVERGENZA      **********************************


     s1 <- solve(sigma)
     s2 <- s1 / (1 + lambda) 
     sigma2 <- (1+lambda)* sigma 
     q1 <- matrix(tensorizza (dif, s1),n,1) 
     q2 <- matrix(tensorizza (dif, s2),n,1) 
     rm (s1,s2)
     gc()
     
     q1 <- -0.5*q1
     q2 <- -0.5*q2
     
     ll <- w1 * exp(q1)  / sqrt(2*pi*det(sigma)) + (1-w1) * (exp(q2)) / sqrt(2*pi*det(sigma2))
     lik <- sum(log(ll))
     
     conv <- (abs(lik-oldlik) > eps*abs(lik-lik0) )

     #alpha <- sqrt((lambda+1) )
     oldlik <- lik
     if (iter == 1)  
        lik0 <- lik
   
  }
 #************************  FINE CICLO EM  ************************************
    gc() 

#   CALCOLO DEI VALORI PREVISTI
    yprev <-  pred.y(y=memo.y, x=memo.x, B, sigma, lambda, w=1-w1, model = model)    
#   CALCOLO DEL BIC
# N. parametri per il modello normale
    k1 <- ncol(x) * p + (p*(p+1))/2 # p=ncol(y) 
# N. parametri per il modello di contaminazione
    k2 <- k1 + 2 - w.fix - lambda.fix 
# Calcolo della verisimiglianza normale

    dati<-cbind(x,y)
    q <- ncol(x)
    norm.mv<-function(u){dmvnorm(u[q+1:p], t(B0)%*%u[1:q], sigma0, log=TRUE)}
    lik.n <- sum(apply(dati,1,norm.mv))
   ###############  calcolo di BIC e AIC per i due modelli #############
   
   BIC.n <- -2*lik.n + k1*log(n)
   BIC.mix <- -2*lik + k2*log(n)
   AIC.n <- 2*k1 - lik.n   
   AIC.mix <- 2*k2 - lik   
    
    ris <- list(
                 ypred = as.matrix(yprev[,1:(ncol(yprev)-3)]),  
                 B=B, 
                 sigma=sigma, 
                 lambda=lambda,
                 w=1-w1, 
                 tau=yprev$tau, 
                 outlier = yprev$outlier,
                 pattern= yprev$pattern,
                 is.conv = (iter < max.iter),
                 n.iter =iter,
                 bic.aic = c(BIC.norm=BIC.n, BIC.mix=BIC.mix, AIC.norm=AIC.n, AIC.mix=AIC.mix)
               )

    attr(ris, "model") <- model
  
    class(ris) <- c(class(ris), "mlest" )
    if (iter == max.iter)
      warning (paste("EM algorithm failed to converge: stop after", max.iter, "iterations"))
    ris
 }

