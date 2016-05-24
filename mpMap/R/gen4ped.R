gen4ped <-
function(nfunnels=1, nperfam=50, nssdgen=6, nseeds=1, iripgen=0)
{
  if (!(nfunnels %in% c(1,3))) 
	stop("Only intended to generate pedigrees with 1 or 3 funnels")

  obs <- vector()
  # start with founders
  ped <- rbind(c(1,0,0), c(2,0,0), c(3,0,0), c(4,0,0))
  ped <- rbind(ped, c(5,1,2), c(6,3,4))
  if (nfunnels==3)
 	ped <- rbind(ped, c(7,1,3), c(8,2,4), c(9,1,4), c(10,2,3))
  
  n1 <- nrow(ped)+1
 
  for (j in 1:nfunnels)
	ped <- rbind(ped, cbind(c((nrow(ped)+1):(nrow(ped)+nperfam)), rep((j-1)*2+5, nperfam), rep((j-1)*2+6, nperfam)))
 
  n2 <- nrow(ped)  

  # at this point have done all the mixing, need to do AI and then selfing
  if (iripgen>0)
  for (i in 1:iripgen)
  {
    for (j in n1:n2)
    	ped <- rbind(ped, c(nrow(ped)+1, j, sample(setdiff(n1:n2, j), 1)))
    n1 <- n2+1
    n2 <- nrow(ped)
  }

  obs <- rep(0, nrow(ped))
  for (i in 1:(nperfam*nfunnels))
  { 
    index <- i+n1-1

    for (j in 1:nseeds)
    {
	obs <- c(obs, rep(0, nssdgen-1), 1)
    	ped <- rbind(ped, c(nrow(ped)+1, index, index))
	ped <- rbind(ped, cbind(c((nrow(ped)+1):(nrow(ped)+nssdgen-1)), c(nrow(ped):(nrow(ped)+nssdgen-2)), c(nrow(ped):(nrow(ped)+nssdgen-2))))
    }
  }

  # fourth column is whether individual was genotyped
  ped <- cbind(ped, obs)
  ped <- as.data.frame(ped)
  names(ped) <- c("id", "m", "f", "obs")

  return(ped)
}

