adjCov <- function(x, id){


if(x$VC$triv == TRUE){


if( !is.null(x$X3) ) {mul2 <- x$X4; mul3 <- x$X5; mul4 <- x$X6} 
if(  is.null(x$X3) )  mul2 <- mul3 <- mul4 <- 1 
                                       
scores <- cbind( c(x$fit$dl.de1)*x$X1, 
                 c(x$fit$dl.de2)*x$X2, 
                 c(x$fit$dl.de3)*x$X3,
                 c(x$fit$dl.dtheta12.st)*mul2,
                 c(x$fit$dl.dtheta13.st)*mul3,
                 c(x$fit$dl.dtheta23.st)*mul4  ) 

}




if(x$VC$triv == FALSE){


cont2par <- x$VC$m2   
cont3par <- x$VC$m3  
bin.link <- x$VC$bl  

Vb <- x$Vb 

if(x$Cont == "NO"){ ###


if( x$margins[2] %in% bin.link && x$Model != "BPO0"){

if(is.null(x$X3) )  mul <- 1
if(!is.null(x$X3) ) mul <- x$X3

scores <- cbind( c(x$fit$dl.dbe1)*x$X1, c(x$fit$dl.dbe2)*x$X2, c(x$fit$dl.drho)*mul )

                                                    }

if( x$Model == "BPO0" ){

scores <- cbind( c(x$fit$dl.dbe1)*x$X1, c(x$fit$dl.dbe2)*x$X2 )

                       }



if( x$margins[2] %in% cont2par ){

if( !is.null(x$X3) && !is.null(x$X4) ) {mul1 <- x$X3; mul2 <- x$X4} 
if(  is.null(x$X3) &&  is.null(x$X4) )  mul1 <- mul2 <- 1 
                                       
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*x$X2, 
                 c(x$fit$dl.dsigma.st)*mul1,
                 c(x$fit$dl.dteta.st)*mul2       )                                           
                                 }

if( x$margins[2] %in% cont3par ){

if( !is.null(x$X3) && !is.null(x$X4) && !is.null(x$X5)) {mul1 <- x$X3; mul2 <- x$X4; mul3 <- x$X5} 
if(  is.null(x$X3) &&  is.null(x$X4) &&  is.null(x$X5))  mul1 <- mul2 <- mul3 <- 1 
                                       
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*x$X2, 
                 c(x$fit$dl.dsigma.st)*mul1,
                 c(x$fit$dl.dnu.st)*mul2,
                 c(x$fit$dl.dteta.st)*mul3       )     
                                       
                                }




} ###



if(x$Cont == "YES"){ ###




if( x$margins[1] %in% cont2par && x$margins[2] %in% cont2par ){

if( !is.null(x$X3) ) {mul1 <- x$X3; mul2 <- x$X4; mul3 <- x$X5} 
if(  is.null(x$X3) )  mul1 <- mul2 <- mul3 <- 1 
                                       
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*x$X2, 
                 c(x$fit$dl.dsigma21.st)*mul1,
                 c(x$fit$dl.dsigma22.st)*mul2,
                 c(x$fit$dl.dteta.st)*mul3       ) 
                 
                 
}






if( x$margins[1] %in% cont3par && x$margins[2] %in% cont3par ){

if( !is.null(x$X3) ) {mul1 <- x$X3; mul2 <- x$X4; mul3 <- x$X5; mul4 <- x$X6; mul5 <- x$X7} 
if(  is.null(x$X3) )  mul1 <- mul2 <- mul3 <- mul4 <- mul5 <- 1 
                                       
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*x$X2, 
                 c(x$fit$dl.dsigma21.st)*mul1,
                 c(x$fit$dl.dsigma22.st)*mul2,
                 c(x$fit$dl.dnu1.st)*mul3,
                 c(x$fit$dl.dnu2.st)*mul4,
                 c(x$fit$dl.dteta.st)*mul5       ) 
                 
                 
}






if( x$margins[1] %in% cont2par && x$margins[2] %in% cont3par ){

if( !is.null(x$X3) ) {mul1 <- x$X3; mul2 <- x$X4; mul3 <- x$X5; mul4 <- x$X6} 
if(  is.null(x$X3) )  mul1 <- mul2 <- mul3 <- mul4 <- 1 
                                       
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*x$X2, 
                 c(x$fit$dl.dsigma21.st)*mul1,
                 c(x$fit$dl.dsigma22.st)*mul2,
                 c(x$fit$dl.dnu2.st)*mul3,
                 c(x$fit$dl.dteta.st)*mul4       ) 
                 
                 
}




if( x$margins[1] %in% cont3par && x$margins[2] %in% cont2par ){

if( !is.null(x$X3) ) {mul1 <- x$X3; mul2 <- x$X4; mul3 <- x$X5; mul4 <- x$X6} 
if(  is.null(x$X3) )  mul1 <- mul2 <- mul3 <- mul4 <- 1 
                                       
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*x$X2, 
                 c(x$fit$dl.dsigma21.st)*mul1,
                 c(x$fit$dl.dsigma22.st)*mul2,
                 c(x$fit$dl.dnu1.st)*mul3,
                 c(x$fit$dl.dteta.st)*mul4       ) 
                 
                 
}

 

}



} # end TRIV










scores <- aggregate.data.frame(scores,by=list(id),FUN=sum)[,-1]
nclusters <- dim(scores)[1]
meat   <- (nclusters-1)*var(scores)
covsan <- Vb %*% meat %*% Vb
x$Vb <- covsan

rm(scores, nclusters, meat, covsan, Vb)

x

}


