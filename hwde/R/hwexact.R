`hwexact` <-
function(obs.hom1=68, obs.hets=28, obs.hom2=4, data=NULL, loci="locus1",
           observed="Observed")
{
  if(is.data.frame(obs.hom1))data=obs.hom1
  if(!is.null(data)){
    loc <- data[, loci]
    types <- decode.genotypes(loc)$types
    obs.hom1 <- data[loc==types[1], observed]
    obs.hom2 <- data[loc==types[3], observed]
    obs.hets <- data[loc==types[2], observed]
  }
  if(length(obs.hom1) != length(obs.hets) | length(obs.hets) != length(obs.hom2))
    stop('length of obs.hom1, obs.hets, and obs.hom2 must be the same')

  return(unlist(lapply(1:length(obs.hom1), function(i){snphwe(obs.hom1[i], obs.hets[i], obs.hom2[i])})))
}

