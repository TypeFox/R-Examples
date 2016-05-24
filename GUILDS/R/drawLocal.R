
pm_sad <- function(th,I,j) {

	k1 = unique(round(10^seq(0,log10(j-1),len=1000)))     #k1=unique(round(logspace(0,log10(j-1),1000)));
	res=rep(0,length(k1))  								  #res=zeros(1,length(k1));
	for(cnt in 1:length(k1)) {
		k = k1[cnt];
		f <- function(x) {
			return(-1 * pm_sadaux(x,I,th,j,k));
		}
		xmax = optimize(f,c(0,1))$minimum
		ymax = pm_sadaux(xmax,I,th,j,k)
		y0 = pm_sadaux(0,I,th,j,k);
		if(ymax<(y0+1)) {
			xmin = 0;
		} else {
			f <- function(x) {
				return(pm_sadaux(x,I,th,j,k) - ymax + 1);
			}
			#xL = f(1e-10);
			#xR = f(xmax);
			#cat(cnt,xL,xR,"\n"); flush.console();
			xmin = uniroot(f,c(1e-10,xmax))$root;
		}
		x1 = (xmax+1)/2;
		y1 = pm_sadaux(x1,I,th,j,k);
		while(ymax < (y1+1)) {
			x1 = (x1 + 1) / 2;
			y1 = pm_sadaux(x1,I,th,j,k);
		}
		f <- function(x) {
				return(pm_sadaux(x,I,th,j,k)-ymax+1);
		}
		
		xplu = uniroot(f,c(xmax,x1),tol=eps())$root;
		xlft = xmax - 10 * (xmax - xmin); xlft = max(0,xlft);
		xrgt = xmax + 10 * (xplu - xmax); xrgt = min(1,xrgt);
		
		f <- function(x) {
		    return( exp(pm_sadaux(x,I,th,j,k)));
		}
		
		res[cnt] = integrate(f,xlft,xrgt,abs.tol=1e-9)$val  #checked!!!!!1
	}

	k2 = 1:(j-1);  #the original algorithm had 1:j, but at the interp1 algorithm does not evaluate at the ends.
	kesk = interp1(x=log2(k1),y=k1*res,xi=log2(k2),method="spline"); 	#kesk = spline(log2(k1),k1.*res,log2(k2));
	k2 <- c(k2,j);
	kesk <- c(kesk,0);
	esk = kesk / k2;
	esk[esk<eps()] = eps();
	return(esk);
}

pm_sadaux <- function(x,I,th,j,k) {
    #cat(x,I,th,j,k,"\n"); flush.console();
	y = matrix(0,nrow=size(x)[1],ncol=size(x)[2]);
	idx0=which(x==0);
	if(length(idx0)) {
		aux = log(th) + lgamma(I) - lgamma(I+j) +
			  lgamma(j+1) - lgamma(k+1) - lgamma(j-k+1) +
			  lgamma(k) +   lgamma(I+j-k) - lgamma(I) + log(I);
		y[idx0] = aux;
	}
	idx1=which(x==1);
	if(length(idx1)) {
		  if ((k<j)||((k==j)&&(th>1))) {
			  y[idx1] = -Inf;
		  } else if ((k==j)&&(th==1)) {
		      y[idx1] = 0;
		  } else if ((k==j)&&(th< 1)) {
			  y[idx1] = Inf;
		   }
	}
	idxp = which((x>0)&(x<1));
	xx=x[idxp];
	aux = log(th) + lgamma(I)-lgamma(I+j) +
		  lgamma(j+1) - lgamma(k+1) - lgamma(j-k+1) +
		  lgamma(I*xx+k) - lgamma(I*xx) +
		  lgamma(I*(1-xx)+j-k)-lgamma(I*(1-xx)) +
		  (th-1)*log(1-xx)-log(xx);
	y[idxp] = aux;
    return(y);
}

	
	




drawLocalCond <- function(theta,alpha_x,alpha_y,JX,JY) {
    #make two SADs
  J = JX + JY;
  SADX <- c();
  SADY <- c();

  nx = getPX(theta,alpha_x,alpha_y,JX,JY);
  ny = 1 - nx; 

  I_X = alpha_x * nx * (J-1) / (1 - alpha_x*nx-alpha_y*ny); #update I_X and I_Y accordingly
  I_Y = alpha_y * ny * (J-1) / (1 - alpha_x*nx-alpha_y*ny);

  SADX = pm_sad(theta, I_X, JX);
  SADY = pm_sad(theta, I_Y, JY);

  output <- list( guildX = SADX,guildY = SADY);
  
  return(output);
}

drawLocal <- function(theta,alpha_x,alpha_y,J) {
    #make two SADs
  SADX <- c();
  SADY <- c();

  nx = rbeta(1,theta,theta);
  ny = 1 - nx; 

  I_X = alpha_x * nx * (J-1) / (1 - alpha_x*nx-alpha_y*ny); #update I_X and I_Y accordingly
  I_Y = alpha_y * ny * (J-1) / (1 - alpha_x*nx-alpha_y*ny);

  probs <- c();
  allN <- 0:J;
  if(is.infinite(I_X) && is.infinite(I_Y)) {
     probs = exp( lgamma(J+1) - (lgamma(allN + 1) + lgamma(J - allN + 1)) + allN * log(nx) + (J-allN) * log(ny));
  } else {
    probs = PolyaEggenberger(I_X,I_Y,J,allN); #set up a probability vector
  }

  NX = sample(0:J,1,replace=TRUE,prob=probs);
  NY = J - NX;

  SADX = pm_sad(theta, I_X, NX);
  SADY = pm_sad(theta, I_Y, NY);

  output <- list( guildX = SADX,guildY = SADY);
  
  return(output);
}