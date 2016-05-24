logLikguilds <- function(theta_x,theta_y,alpha_x,alpha_y,J,Sx,Sy,Nx,Ny,KDA_X,KDA_Y,prefactor1,prefactor2,verbose=TRUE) {
 thrs <- 10;
 
 f <- function(x) {
   return( -1 * evaluateLogLik(x,theta_x,theta_y,alpha_x,alpha_y,J,Sx,Sy,Nx,Ny,KDA_X,KDA_Y) );
 }

 maxes <- fminbnd(f,0,1,maxiter=500,tol=1e-4);

 ymax = -1 * maxes$fmin;
 xmax = maxes$xmin;
 xlft = 0;
 xrgt = 1;
 eps = .Machine$double.eps; 

 checkLeftX  = evaluateLogLik(eps,theta_x,theta_y, alpha_x,alpha_y,J,Sx,Sy,Nx,Ny,KDA_X,KDA_Y)  - ymax + thrs;
 checkRightX = evaluateLogLik(1-eps,theta_x,theta_y, alpha_x,alpha_y,J,Sx,Sy,Nx,Ny,KDA_X,KDA_Y)- ymax + thrs;

 if(checkLeftX < 0) {
     g <- function(x) {
        return(evaluateLogLik(x,theta_x,theta_y,alpha_x,alpha_y,J,Sx,Sy,Nx,Ny,KDA_X,KDA_Y) - ymax + thrs);
     }
     xlft <- (uniroot(g,c(eps,xmax)))$root;
 }

 if(checkRightX < 0) {
     h <- function(x) {
        return(evaluateLogLik(x,theta_x,theta_y,alpha_x,alpha_y,J,Sx,Sy,Nx,Ny,KDA_X,KDA_Y) - ymax + thrs);
     }
     xrgt <- (uniroot(h,c(xmax,1-eps)))$root;
 }

 calcLLEXP <- function(x) {
   out <- exp(evaluateLogLik(x,theta_x,theta_y,alpha_x,alpha_y,J,Sx,Sy,Nx,Ny,KDA_X,KDA_Y)-ymax);
   return(out);
 }

 aux <- integrate(f=calcLLEXP,lower=xlft,upper=xrgt,abs.tol=1e-9);

 y = ymax + log(aux$value) + prefactor1 + Sx * log(theta_x/2) + prefactor2 + Sy * log(theta_y/2) + lgamma((theta_x/2) + (theta_y/2)) - (lgamma(theta_x/2) + lgamma(theta_y/2));
 
 if(verbose==TRUE) cat(sprintf("%4.3f       %4.4f       %4.4f       %4.4f        %4.3f\n",theta_x,theta_y,alpha_x,alpha_y,y));flush.console(); 
 
 return(y);
}


maxLikelihood.Guilds <- function(initVals,model,method,SADX,SADY,verbose=TRUE) {
  KDA_X <- calcKDA(SADX);
  KDA_Y <- calcKDA(SADY);
    
  x <- c(table(SADX));
  freq_x <- c();
  for(i in 1:length(x)) freq_x[i] <- x[[i]];
  prefactor1 = -( sum(log(SADX)) + sum(lgamma(1+freq_x)) );
  
  x2 <- c(table(SADY));
  freq_y <- c();
  for(i in 1:length(x2)) freq_y[i] <- x2[[i]];
  prefactor2 = -( sum(log(SADY)) + sum(lgamma(1+freq_y)) );

  Sx = length(SADX);
  Sy = length(SADY);
  Nx = sum(SADX);
  Ny = sum(SADY);
  J = Nx + Ny;
  
  loglikverbose <- verbose;
  if(method == "simplex") loglikverbose <- FALSE;

  g <- function(v) {  
  	incorrectLength <- 0;
  	if(model == "D0" && length(v) != 2) incorrectLength <- 1;
  	if(model == "D1" && length(v) != 3) incorrectLength <- 1;
  
  	if(incorrectLength == 1) { cat("Input vector is of incorrect length\n"); return; };

  	if(model == "D0") { theta_x = v[1]; theta_y = v[1]; alpha_x = v[2]; alpha_y = v[2];}
  	if(model == "D1") { theta_x = v[1]; theta_y = v[1]; alpha_x = v[2]; alpha_y = v[3];}
	
  	if(alpha_x<0||alpha_y<0||theta_x<1||theta_y<1) return(-Inf)
  	if(alpha_x>(1-(1e-8))||alpha_y>(1-(1e-8))) return(-Inf)
  
  	y = logLikguilds(theta_x,theta_y,alpha_x,alpha_y,J,Sx,Sy,Nx,Ny,KDA_X,KDA_Y,prefactor1,prefactor2,verbose);
  	return(y);
  }


  x;
  if(method == "simplex") {
    	x <- simplex(initVals,g,verbose);
  }
  if(method == "subplex") {
	x <- subplex(initVals,g);
  }
  return(x);
}

logLikelihood.Guilds <- function(parameters,model,SADX,SADY,verbose=TRUE) {

  KDA_X <- calcKDA(SADX);
  KDA_Y <- calcKDA(SADY);
  
  x <- c(table(SADX));
  freq_x <- c();
  for(i in 1:length(x)) freq_x[i] <- x[[i]];
  prefactor1 = -( sum(log(SADX)) + sum(lgamma(1+freq_x)) );
  
  x2 <- c(table(SADY));
  freq_y <- c();
  for(i in 1:length(x2)) freq_y[i] <- x2[[i]];
  prefactor2 = -( sum(log(SADY)) + sum(lgamma(1+freq_y)) );

  Sx = length(SADX);
  Sy = length(SADY);
  Nx = sum(SADX);
  Ny = sum(SADY);
  J = Nx + Ny;
  
  
  if(model == "D0") { theta_x = parameters[1]; theta_y = parameters[1]; alpha_x = parameters[2]; alpha_y = parameters[2];}
  if(model == "D1") { theta_x = parameters[1]; theta_y = parameters[1]; alpha_x = parameters[2]; alpha_y = parameters[3];}

  LL <- -Inf;
  
  if(verbose==TRUE) {
   cat("Chosen model: ",model,"\n");
   cat("Now starting to calculate likelihood of: \n"); x2 <- parameters;
   if(model == "D0") cat("Theta X =",theta_x," Theta Y =","Theta X","\t Alpha X =", alpha_x," Alpha Y =","Alpha X" ,"\n");
   if(model == "D1") cat("Theta X =",theta_x," Theta Y =","Theta X","\t Alpha X =", alpha_x," Alpha Y =", alpha_y    ,"\n"); 

   flush.console();
  }
  


  LL <- logLikguilds(theta_x,theta_y,alpha_x,alpha_y,J,Sx,Sy,Nx,Ny,KDA_X,KDA_Y,prefactor1,prefactor2,verbose=);
  
  if(verbose==TRUE) cat("Likelihood is ",LL,"\n");
  
  return(LL);
}
