generate.Guilds <- function(theta,alpha_x, alpha_y, J) {

  #make two SADs
  SADX <- c();
  SADY <- c();

  #first draw nx and ny from a beta distribution
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

  SADX = generate.ZSM(theta,I_X,NX);
  SADY = generate.ZSM(theta,I_Y,NY);

  output <- list( guildX = SADX,guildY = SADY);
  
  return(output);
}


PolyaEggenberger <- function(theta_x,theta_y,J,N) {
  a = lgamma(J+1); #J!
  b = lgamma(theta_x + theta_y + J) - lgamma(theta_x + theta_y); #(theta_x + theta_y)_J pochhammer
  
  c1 = lgamma(theta_x + N) - lgamma(theta_x); #(theta_x)_Nx  pochhammer
  c2 = lgamma(theta_y + J - N) - lgamma(theta_y); #(theta_y)_Ny  pochhammer
  d = lgamma(N + 1) + lgamma(J-N+1);
  
  return( exp(a - b + c1 + c2 - d));
}

localComm <- function(alpha_x,alpha_y,Jx,Jy,px) {
  J = Jx + Jy;

  nx = px;
  ny = 1 - nx;

  I_X = alpha_x * nx * (J-1) / (1 - alpha_x*nx-alpha_y*ny); #update I_X and I_Y accordingly
  I_Y = alpha_y * ny * (J-1) / (1 - alpha_x*nx-alpha_y*ny);

 if(is.infinite(I_X) && is.infinite(I_Y)) {
     output = exp( lgamma(J+1) - (lgamma(Jx + 1) + lgamma(J - Jx + 1)) + Jx * log(nx) + (J-Jx) * log(ny));
 } else {
  output <- PolyaEggenberger(I_X,I_Y,J,Jx);
 }
  return(output);
}

rho <- function(theta,px) {
  a <- lgamma(2 * theta) - 2 * lgamma(theta)
  b <- (theta - 1) * log(px) + (theta - 1) * log((1-px));
  output <- exp(a+b);
  return(output);
}

getPX <- function(theta,alpha_x,alpha_y,JX,JY) {
   
   px = (1:(1000-1))/1000
   calcLocal <- function(x) {
	a <- localComm(alpha_x,alpha_y,JX,JY,x) * rho(theta,x)
	return(a);
   }
   divider <- integrate(f=calcLocal,lower=0,upper=1,abs.tol=1e-9)$value;


   probs <- calcLocal(px) / divider;

   return(sample(px,1,prob=probs));
}

generate.Guilds.Cond <- function(theta,alpha_x,alpha_y,JX,JY) {

  J = JX + JY;
  SADX <- c();
  SADY <- c();

  nx = getPX(theta,alpha_x,alpha_y,JX,JY);
  ny = 1 - nx; 

  I_X = alpha_x * nx * (J-1) / (1 - alpha_x*nx-alpha_y*ny); #update I_X and I_Y accordingly
  I_Y = alpha_y * ny * (J-1) / (1 - alpha_x*nx-alpha_y*ny);

  SADX = generate.ZSM(theta,I_X,JX);
  SADY = generate.ZSM(theta,I_Y,JY);

  output <- list( guildX = SADX,guildY = SADY);
  
  return(output);

}


generate.ZSM <- function(theta,I,J) {

  #based on urn.gp code by Rampal Etienne available as a supplementary
  #online appendix for his 2005 paper in Ecology Letters

  if(is.infinite(I)) I = .Machine$double.xmax

  anc_ct=0; spec_ct=0;
  ancestor=rep(NaN,J);
  species=rep(NaN,J);
  ancestor_species=rep(NaN,J);

  for(j in 1:J) { #loop over each invidual 
     if( runif(1,0,1) <= I/(I+j-1) ) { #ancestor not previously encountered
        anc_ct = anc_ct + 1;
        ancestor[j] = anc_ct;
        if( runif(1,0,1) <= theta / (theta+anc_ct - 1)) {#new ancestor is a new species
            spec_ct = spec_ct + 1;
            species[j] = spec_ct;
            ancestor_species[anc_ct] = spec_ct;
         } else { #new ancestor is an existing species but new migration into local community
            prior_lineage = ceil(runif(1,0,1) * (anc_ct-1));
            s = ancestor_species[prior_lineage];
            species[j] = s;
            ancestor_species[anc_ct] = s;
         }
     } else { #descendant of existing individual
        descend_from = ceil(runif(1,0,1) * (j-1)); #select one individual at random
        ancestor[j] = ancestor[descend_from];
        species[j] = species[descend_from];
     }
  }
  
  x <- c(table(species));
  abund <- c();
  for(i in 1:length(x))  abund[i] <- x[[i]];
  abund=sort(abund,decreasing=T);
  return(abund);
}

generate.ESF <- function(theta,I,J) {
  return(generate.ZSM(theta,I,J));
}