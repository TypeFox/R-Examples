check.vars <- function (y, x, model="LN", parent="pred.y") {
#print(parent.env(env) )  vedere come posso conoscere il nome della funzione chiamante (altrimenti lo passo)
   EPS <- 0.0000001
   lista  <- list(ret=0, msg.err=NULL, y=NULL, x=NULL)
#   ERRORI BLOCCANTI
   if (!is.numeric(as.matrix(y))  )  {
      lista$msg.err <- "Variables must be numeric "
      lista$ret <- -9
   }
   else if (!inherits(y,c("data.frame", "matrix","numeric", "integer"))  )  {
      lista$msg.err <- "Variables must be supplied as a matrix or as a data frame "
      lista$ret <- -9
   }
    else if ( model == 'LN' & sum(y<0, na.rm=TRUE) )  {
        lista$msg.err <-"Negative values are not allowed if model is LN\n"
        lista$ret <- -9
    }

   y <- as.matrix (y)
   n <- nrow(y)
   if (!is.null(x) & length(x) > 0) {
      x<- as.matrix(x)
      if (!is.numeric(x) )  {
       lista$msg.err <- "Covariates must be numeric "
       lista$ret <- -9
     }   
     else if (!inherits(x,c("data.frame", "matrix","numeric", "integer")) )  {
      lista$msg.err <- "Covariates must be supplied as a matrix or as a data frame "
      lista$ret <- -9
      }
     else if (sum(is.na(x)) > 0)  {
        lista$msg.err <-"Covariates can not have missing values "
        lista$ret <- -9
     } else if (nrow(x) != nrow(y))  {
        lista$msg.err <-"Variables y and x must have the same number of rows"
        lista$ret <- -9
     }
     else if ( model == 'LN' &  sum(x<0) != 0 )  {
        lista$msg.err <-"Negative values are not allowed if model is LN\n"
        lista$ret <- -9
     }
   }
    
   
   if (lista$ret == -9)
       return (lista)
#    WARNING




   ind0<-NULL
#   PREPARAZIONE DATI IN BASE AL MODELLO
  if (model == "LN") {
      ind0 <- which(y == 0)
      if (length(ind0)> 0)  {
        y[ind0] <- EPS
        if (parent!="ml.est")  {
          lista$ret <- lista$ret+1
          lista$msg.err <-  rbind( lista$msg.err,                 
                paste(length(ind0)," response variable values (%", round((length(ind0)*100)/n, 2),
                       ") equal 0 are substituted by ",EPS, "\n",sep="") )
        }
      }
      y <- as.matrix(log(y))
      if (length(x)>0)  {
#        x <- as.matrix (x)
        if ( sum( x == 0) > 0)   {
            ind <- which(x==0)
            x<- as.matrix(x )
            x[ind] <- EPS  
            if (parent!="ml.est")   {                    
              lista$ret <- lista$ret+1
              lista$msg.err <- rbind( lista$msg.err,                  
                paste(length(ind)," covariate values (%", round((length(ind)*100)/length(x), 2),
                     ") equal 0 are substituted by ",EPS, "\n", sep="") )
            }
        }        
        x <- log(x)
      }
      else
        x <- NULL
  } else {
     y <- as.matrix(y)
     x <- x
  }
   lista$y <- y
   lista$x <- as.matrix(cbind(rep(1,n),x))
   return (lista)
 }
 
 
  