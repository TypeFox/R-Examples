prestonSort <- function(A) {
	output <- rep(0,20);

	for(k in 1:length(output)) {
	  start <- 2^(k-1);
	  end <- -1 + 2^(k);
	  if(end > length(A)) {
		end <- length(A);
		X <- sum(A[start:end]);
		output[k] <- X; 
		break;
	  }
	  X <- sum(A[start:end]);
	  output[k] <- X;
	}
	toRemove <- c();
	for(k in length(output):1) {
		if(output[k] < 1e-6) toRemove <- c(toRemove,k);
		if(output[k] > 1e-6) break;
	}
	output <- output[-toRemove];	

	return(output);
}

expected.SAD <- function(theta,m,J) {
  I = (J-1)* m / (1-m)
  aux <- pm_sad(theta,I,J);
  SAD <- prestonSort(aux);
  return(SAD);  
}

expected.SAD.Guilds <- function(theta,alpha_x,alpha_y,J,n_replicates=100) {
  meanX <- rep(0,J);
  meanY <- rep(0,J);

  for(r in 1:n_replicates) {
		M <- drawLocal(theta,alpha_x,alpha_y,J);
		for(m in 1:length(M$guildX)) {
			meanX[m] <- meanX[m] + M$guildX[m];
		}
		for(m in 1:length(M$guildY)) {
			meanY[m] <- meanY[m] + M$guildY[m];
		}
  }
  meanX <- meanX / n_replicates;
  meanY <- meanY / n_replicates;

  gX <- prestonSort(meanX);
  gY <- prestonSort(meanY);
  
  output <- list( guildX = gX,guildY = gY);
  return(output);
}

expected.SAD.Guilds.Conditional <- function(theta,alpha_x,alpha_y,Jx,Jy,n_replicates=100) {
  meanX <- rep(0,Jx);
  meanY <- rep(0,Jy);

  for(r in 1:n_replicates) {
		M <- drawLocalCond(theta,alpha_x,alpha_y,Jx,Jy);
		for(m in 1:length(M$guildX)) {
			meanX[m] <- meanX[m] + M$guildX[m];
		}
		for(m in 1:length(M$guildY)) {
			meanY[m] <- meanY[m] + M$guildY[m];
		}
  }
  meanX <- meanX / n_replicates;
  meanY <- meanY / n_replicates;

  gX <- prestonSort(meanX);
  gY <- prestonSort(meanY);
  
  output <- list( guildX = gX,guildY = gY);
  return(output);
}