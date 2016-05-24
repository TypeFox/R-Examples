gen8ped <-
function(nfunnels=1, nperfam=50, nssdgen=6, nseeds=1, iripgen=0, seed=1)
{
  # could allow nseeds to be a vector, in which case it should have length ?  - one entry for each stage. if length 1, assume that  
 
  ped <- cbind(c(1:8), rep(0, 8), rep(0, 8))

  # special circumstance is just the ABCDEFGH cross
  if (nfunnels==1)
  {
    # create 2-way crosses
    ped <- rbind(ped, c(9, 1, 2), c(10, 3, 4), c(11, 5, 6), c(12, 7, 8))

    # create 4-way crosses
    ped <- rbind(ped, cbind(c(13:(12+nperfam*2)), rep(c(9,11), nperfam), rep(c(10, 12), nperfam)))
  
    n4 <- nrow(ped)+1
    # create 8-way crosses
    ped <- rbind(ped, cbind(c((nrow(ped)+1):(nrow(ped)+nperfam)), seq(13, nrow(ped), 2), seq(14, nrow(ped), 2)))

    n8 <- nrow(ped) 
 
    # irip (note: not really the correct formulation of irip, assumes infinite
    # population from random funnels, which with one funnel is inappropriate)
    if (iripgen>0)
    for (i in 1:iripgen)
    {
      for (j in n4:n8)
	ped <- rbind(ped, c(nrow(ped)+1, j, sample(setdiff(n4:n8, j), 1)))
      n4 <- n8+1
      n8 <- nrow(ped)
    }

    obs <- rep(0, nrow(ped))

    # self
    for (i in 1:nperfam)
    {
     	index <- i+n4-1
	for (j in 1:nseeds)
	{
	  obs <- c(obs, rep(0, nssdgen-1), 1)
	  ped <- rbind(ped, c(nrow(ped)+1, index, index))
	  ped <- rbind(ped, cbind(c((nrow(ped)+1):(nrow(ped)+nssdgen-1)), c(nrow(ped):(nrow(ped)+nssdgen-2)), c(nrow(ped):(nrow(ped)+nssdgen-2))))  
       	}
    }
  } # end of nfunnels==1
	

  # otherwise take a sample of all possible funnels
  if (nfunnels>1)
  {
    # Generate matrix of all possible funnels
    pairs <- combn(1:8, 2)

    quads <- vector()
    for (i in 1:(ncol(pairs)-3))
    	quads <- cbind(quads, rbind(matrix(rep(pairs[,i], choose(6-min(pairs[,i])+1, 2)), nrow=2, ncol=choose(6-min(pairs[,i])+1, 2)), combn(setdiff((min(pairs[,i])+1):8, pairs[,i]), 2)))
 
    eights <- vector()
    for (i in 1:(ncol(quads)-1))
    {
	addons <- which(lapply(apply(quads[,i:ncol(quads)], 2, function(x) return(intersect(quads[,i], x))), length)==0)
 	eights <- cbind(eights, rbind(matrix(rep(quads[,i], length(addons)), nrow=4, ncol=length(addons)), quads[,(i:ncol(quads))[addons]]))
    }

    whichfun <- sample(ncol(eights), nfunnels)

    n1 <- 8
    # create 2-way crosses
    for (i in 1:nfunnels)
      	ped <- rbind(ped, cbind(c((nrow(ped)+1):(nrow(ped)+4)), eights[seq(1, 8, 2) ,whichfun[i]], eights[seq(2, 8, 2), whichfun[i]]))
    n2 <- nrow(ped)
 
    # create 4-way crosses
    ped <- rbind(ped, cbind(c((nrow(ped)+1):(nrow(ped)+2*nfunnels*nperfam)), rep(seq(n1+1, nrow(ped), 2), each=nperfam), rep(seq(n1+2, nrow(ped), 2), each=nperfam)))

    n4 <- nrow(ped)+1
    # create 8-way crosses
    ped <- rbind(ped, cbind(c((nrow(ped)+1):(nrow(ped)+nfunnels*nperfam)), n2+rep(1:nperfam, nfunnels)+rep(2*nperfam*(0:(nfunnels-1)), each=nperfam), n2+nperfam+rep(1:nperfam, nfunnels)+rep(2*nperfam*(0:(nfunnels-1)), each=nperfam)))

    n8 <- nrow(ped)

    # irip
    if (iripgen>0)
    for (i in 1:iripgen)
    {
      for (j in n4:n8)
	ped <- rbind(ped, c(nrow(ped)+1, j, sample(setdiff(n4:n8, j), 1)))
      n4 <- n8+1
      n8 <- nrow(ped)
    }

    obs <- rep(0, nrow(ped))
    
    # self
    for (i in 1:(nperfam*nfunnels))
    {
     	index <- i+n4-1
	for (j in 1:nseeds)
	{
	  obs <- c(obs, rep(0, nssdgen-1), 1)
	  ped <- rbind(ped, c(nrow(ped)+1, index, index))
	  ped <- rbind(ped, cbind(c((nrow(ped)+1):(nrow(ped)+nssdgen-1)), c(nrow(ped):(nrow(ped)+nssdgen-2)), c(nrow(ped):(nrow(ped)+nssdgen-2))))  
       	}
    }
  }

  # fourth column is whether individual was genotyped
  ped <- cbind(ped, obs)
  ped <- as.data.frame(ped)
  names(ped) <- c("id", "m", "f", "obs")

  return(ped)
}
