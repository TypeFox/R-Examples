evaluateCondLik <- function(v,theta_x,theta_y,alpha_x,alpha_y,Nx,Ny) {
  nx = v;
  ny = 1-nx;
  J = Nx+Ny;
  I_X = alpha_x * nx * (J-1) / (1 - alpha_x * nx - alpha_y * ny);
  I_Y = alpha_y * ny * (J-1) / (1 - alpha_x * nx - alpha_y * ny);

  a = lgamma(J+1);
  b = rep(0,length(I_X));
  poch_X = rep(0,length(I_X));
  poch_Y = rep(0,length(I_X));
  
  for(cnt in 1:length(I_X)) {
    b[cnt] = lgamma(I_X[cnt] + I_Y[cnt] + J) - lgamma(I_X[cnt] + I_Y[cnt]);
	poch_X[cnt] = lgamma(I_X[cnt] + Nx) - lgamma(I_X[cnt]);
	poch_Y[cnt] = lgamma(I_Y[cnt] + Ny) - lgamma(I_Y[cnt]);
  }
  
  c = a - b;
  
  h = poch_X + poch_Y - (lgamma(Nx+1) + lgamma(Ny+1));
   
  k  =  lgamma((theta_x/2) + (theta_y/2)) - (lgamma(theta_x/2) + lgamma(theta_y/2));
  l  = ((theta_x/2) - 1) * log(nx) + ((theta_y/2) - 1) * log(ny);
  
  result = c + h + k + l;
  return(result);
}


calcConditional <- function(v,model,Nx,Ny) {
  incorrectLength <- 0;
  if(model == "D0" && length(v) != 2) incorrectLength <- 1;
  if(model == "D1" && length(v) != 3) incorrectLength <- 1;
  
  if(incorrectLength == 1) { cat("Input vector is of incorrect length\n"); return; };
  
  if(model == "D0") { theta_x = v[1]; theta_y = v[1]; alpha_x = v[2]; alpha_y = v[2];}
  if(model == "D1") { theta_x = v[1]; theta_y = v[1]; alpha_x = v[2]; alpha_y = v[3];}
 
  if(alpha_x<0||alpha_y<0||theta_x<1||theta_y<1) return(-Inf)
  if(alpha_x>(1-(1e-8))||alpha_y>(1-(1e-8))) return(-Inf)
  
  thrs <- 10;
 
  f <- function(x) {
    return( -1 * evaluateCondLik(x,theta_x,theta_y,alpha_x,alpha_y,Nx,Ny) );
  }
  
  maxes <- fminbnd(f,0,1,maxiter=500,tol=1e-4);

 ymax = -1 * maxes$fmin;
 xmax = maxes$xmin;
 xlft = 0;
 xrgt = 1;
 eps = .Machine$double.eps; 

 checkLeftX  = evaluateCondLik(eps,theta_x,theta_y,alpha_x,alpha_y,Nx,Ny)  - ymax + thrs;
 checkRightX = evaluateCondLik(1-eps,theta_x,theta_y,alpha_x,alpha_y,Nx,Ny)  - ymax + thrs;


 if(checkLeftX < 0) {
     g <- function(x) {
        return(evaluateCondLik(x,theta_x,theta_y,alpha_x,alpha_y,Nx,Ny) - ymax + thrs);
     }
     xlft <- (uniroot(g,c(eps,xmax)))$root;
 }

 if(checkRightX < 0) {
     h <- function(x) {
        return(evaluateCondLik(x,theta_x,theta_y,alpha_x,alpha_y,Nx,Ny) - ymax + thrs);
     }
     xrgt <- (uniroot(h,c(xmax,1-eps)))$root;
 }

 calcLLEXP <- function(x) {
   out <- exp(evaluateCondLik(x,theta_x,theta_y,alpha_x,alpha_y,Nx,Ny)-ymax);
   return(out);
 }

  aux <- integrate(f=calcLLEXP,lower=xlft,upper=xrgt,abs.tol=1e-9);

  LL <- log(aux$value) + ymax;

  return(LL);  
}


logLikelihood.Guilds.Conditional <- function(parameters,model,SADX,SADY,verbose=TRUE) {
  Sx = length(SADX);
  Sy = length(SADY);
  Nx = sum(SADX);
  Ny = sum(SADY);
  J = Nx + Ny;
  
  LL <- -Inf;

  if(verbose==TRUE) {
   cat("Chosen model: ",model,"\n");
   cat("Now starting to calculate likelihood of: \n"); x2 <- parameters;
   if(model == "D0") cat("Theta X =",x2[1]," Theta Y =","Theta X","\t Alpha X =", x2[2]," Alpha Y =","Alpha X" ,"\n");
   if(model == "D1") cat("Theta X =",x2[1]," Theta Y =","Theta X","\t Alpha X =", x2[2]," Alpha Y =", x2[3]    ,"\n"); 

   flush.console();
  }
   
  LL <- logLikelihood.Guilds(parameters,model,SADX,SADY,verbose);
  
  #conditional part:
  conditional_part <- calcConditional(parameters,model,Nx,Ny);
   
  output <- LL - conditional_part;
  return(output);
}

conditional.LogLik <- function(v,model,J,Sx,Sy,Nx,Ny,KDA_X,KDA_Y,prefactor1,prefactor2,verbose=TRUE) {  
  if(model == "D0") { theta_x = v[1]; theta_y = v[1]; alpha_x = v[2]; alpha_y = v[2];}
  if(model == "D1") { theta_x = v[1]; theta_y = v[1]; alpha_x = v[2]; alpha_y = v[3];}
  
  if(alpha_x<0||alpha_y<0||theta_x<1||theta_y<1) return(-Inf)
  if(alpha_x>(1-(1e-8))||alpha_y>(1-(1e-8))) return(-Inf)
  
  
  LL <- logLikguilds(theta_x,theta_y,alpha_x,alpha_y,J,Sx,Sy,Nx,Ny,KDA_X,KDA_Y,prefactor1,prefactor2,verbose);
  cond_LL <- calcConditional(v,model,Nx,Ny);
  out <- LL - cond_LL;
  return(out);
}




maxLikelihood.Guilds.Conditional <- function(initVals,model,method,SADX,SADY,verbose=TRUE) {
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

  g <- function(x) {
    out <- -1 * conditional.LogLik(x,model,J,Sx,Sy,Nx,Ny,KDA_X,KDA_Y,prefactor1,prefactor2,verbose);
    return(out);
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



