quadIT <-
function(nodes=14, mu=0, sigma=1,absrange=5, ngr=1, ergE=NULL) 
{
  # riqv_quer = resulting object of the Estep function. NULL if a parametric prior should be estimated
  # ergE != NULL if nonpar estimation is required
  
  if(all(is.null(ergE)))
  {
    
  nwpgru <- mapply(function(gruu,mug,si)
  {
    
    if(length(nodes)==1)
    {
      
      if(nodes < 7) 
      {
        nodes <- 7
        cat("#nodes == too small --> set to 7\n")
        absrange <- 4
      }
      
      quadP <- seq(absrange*(-1), absrange, length.out = nodes) 
      
      quadweight  <- dnorm(quadP)
      quadweight1 <- quadweight/sum(quadweight)
      quadP_shift <- quadP*si + mug
      
    } else  {
      quadP <- nodes - mean(nodes)
      
      quadweight  <-  dnorm(quadP)
      quadweight1 <- quadweight/sum(quadweight)
      quadP_shift <- quadP*si + mug
    }
    
    list(nodes=quadP_shift,weights=quadweight1)
  },gruu=1:ngr,mug=mu,si=sigma,SIMPLIFY=FALSE)
  return(nwpgru)
    
  } else { # nonparametric things
    
  # loop over GROUPS 
    nwpgru <- mapply( function(pg,gruu)
    {
      At <- rowSums(pg) # estimate empirical weights
      quadnodes <- seq(absrange*(-1), absrange, length.out = nodes) # create nodes 
      
      if(any(At< 0.0001))
        { # if the values get too small (metric = Expected number of persons) set a small number
        At[At < 0.0001] <- 0.0001
        }
      
      if(gruu == 1) # for reference group or in case there is only one group
        {
        QN <- stdrdize_hist(At,quadnodes) # first group gets an standardized histogram
        } else 
          {
            AtdAt <- At / sum(At)
            QN <- list(nodes=quadnodes,weights=AtdAt)
          }
        
      
  return(QN)
      
    }, pg=ergE$fiqG,gruu=1:ngr, SIMPLIFY=FALSE)
  

    
  }

  
  
  
  
}
