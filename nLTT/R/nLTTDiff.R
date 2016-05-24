nLTTstat <- function(tree1, tree2, distanceMethod = "abs")
{
 	diff <- -10;
 	if(distanceMethod == "abs") diff = normLTTdiffABS(tree1,tree2);
 	if(distanceMethod == "squ") diff = normLTTdiffSQ(tree1,tree2);
 	 	
 	if(diff < 0) {cat("chosen unknown distance method!\n"); flush.console();}
 	return(diff);	
}

nLTTstat_exact <- function(tree1, tree2, distanceMethod = "abs")
{
 	diff <- -10;
 	if(distanceMethod == "abs") diff = normLTTdiffexactABS(tree1,tree2);
 	if(distanceMethod == "squ") diff = normLTTdiffexactSQ(tree1,tree2);
 	 	
 	if(diff < 0) {cat("chosen unknown distance method!\n"); flush.console();}
 	return(diff);	
}



normLTTdiffexactABS <- function(tree1,tree2)
{
	b_times <- c(-1 * rev(sort(branching.times(tree1))),0);
  	lineages <- c(2:length(b_times),length(b_times));
  	b_times_N <- 1 - b_times / min(b_times); #normalize branching times  
  	lineages_N <- lineages / max(lineages);  #normalize lineages  	
	
	b_times2 <- c(-1 * rev(sort(branching.times(tree2))),0);
    lineages2 <- c(2:length(b_times2),length(b_times2));  
    b_times2_N <- 1 - b_times2 / min(b_times2); #normalize branching times  
    lineages2_N <- lineages2 / max(lineages2);  #normalize lineages  
	
	allBtimes <- unique(sort(c(b_times_N,b_times2_N)));
	diff <- 0;
	for(k in 2:length(allBtimes))
	{
			tim <- allBtimes[k];
			index1 <- max(which(b_times_N <= tim));
			index2 <- max(which(b_times2_N <= tim));
			lins1 <- lineages_N[index1];
			lins2 <- lineages2_N[index2];
			dt <- allBtimes[k] - allBtimes[k-1]
			diff <- diff + dt * abs(lins1-lins2);
	}
	return(diff);
}

normLTTdiffexactSQ <- function(tree1,tree2)
{
	b_times <- c(-1 * rev(sort(branching.times(tree1))),0);
  	lineages <- c(2:length(b_times),length(b_times));
  	b_times_N <- 1 - b_times / min(b_times); #normalize branching times  
  	lineages_N <- lineages / max(lineages);  #normalize lineages  	
	
	b_times2 <- c(-1 * rev(sort(branching.times(tree2))),0);
    lineages2 <- c(2:length(b_times2),length(b_times2));  
    b_times2_N <- 1 - b_times2 / min(b_times2); #normalize branching times  
    lineages2_N <- lineages2 / max(lineages2);  #normalize lineages  
	
	allBtimes <- unique(sort(c(b_times_N,b_times2_N)));
	diff <- 0;
	for(k in 2:length(allBtimes))
	{
			tim <- allBtimes[k];
			index1 <- max(which(b_times_N <= tim));
			index2 <- max(which(b_times2_N <= tim));
			lins1 <- lineages_N[index1];
			lins2 <- lineages2_N[index2];
			dt <- allBtimes[k] - allBtimes[k-1]
			diff <- diff + dt * (lins1-lins2)*(lins1-lins2);
	}
	return(diff);
}






normLTTdiffABS <- function(tree1, tree2) {

  b_times <- c(-1 * rev(sort(branching.times(tree1))),0);
  lineages <- c(2:length(b_times),length(b_times));
  b_times_N <- 1 - b_times / min(b_times); #normalize branching times  
  lineages_N <- lineages / max(lineages);  #normalize lineages  
  ltt1 <- approxfun(b_times_N,lineages_N,method="constant");

  b_times2 <- c(-1 * rev(sort(branching.times(tree2))),0);
  lineages2 <- c(2:length(b_times2),length(b_times2));  
  b_times2_N <- 1 - b_times2 / min(b_times2); #normalize branching times  
  lineages2_N <- lineages2 / max(lineages2);  #normalize lineages  
  ltt2 <- approxfun(b_times2_N,lineages2_N,method="constant");

  f <- function(t,x,p) { #function f is the absolute difference in time t: 0 <= t < 1
       output <- abs( ltt1(t) - ltt2(t));
       return(list(output));
  }
  
  times <- (0:100)/100;
  int_1 <- lsoda(0,times,func=f,tcrit=c(1)); #integrate over t: 0 < t < 1, notice tcrit=0 indicating t should never be larger than 0.
  total_area <- int_1[length(times),2]
  
  return(total_area);
}

normLTTdiffSQ <- function(tree1, tree2) {

  b_times <- c(-1 * rev(sort(branching.times(tree1))),0);
  lineages <- c(2:length(b_times),length(b_times));
  b_times_N <- 1 - b_times / min(b_times); #normalize branching times  
  lineages_N <- lineages / max(lineages);  #normalize lineages  
  ltt1 <- approxfun(b_times_N,lineages_N,method="constant");

  b_times2 <- c(-1 * rev(sort(branching.times(tree2))),0);
  lineages2 <- c(2:length(b_times2),length(b_times2));  
  b_times2_N <- 1 - b_times2 / min(b_times2); #normalize branching times  
  lineages2_N <- lineages2 / max(lineages2);  #normalize lineages  
  ltt2 <- approxfun(b_times2_N,lineages2_N,method="constant");

  f <- function(t,x,p) { #function f is the absolute difference in time t: 0 <= t < 1
       output <- ltt1(t) - ltt2(t);
	   output <- output * output;
       return(list(output));
  }
  
  times <- (0:100)/100;
  int_1 <- lsoda(0,times,func=f,tcrit=c(1)); #integrate over t: 0 < t < 1, notice tcrit=0 indicating t should never be larger than 0.
  total_area <- int_1[length(times),2]
  
  return(total_area);
}


