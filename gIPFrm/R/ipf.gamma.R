ipf.gamma <-
function(ModelMatrix, ObsTable, gamma, tol, estimand)  
{
    nr <- nrow(ModelMatrix);
    nc <- ncol(ModelMatrix);
    if( length(ObsTable) != nc) stop("Dimensions of the model matrix and the data vector do not match");
    est <- (estimand == "probabilities") + (estimand == "intensities");
    if( est != 1) stop("The estimand is not specified correctly.")
    Total <- sum(ObsTable);
    m1 <- rep(1,nc);    ### Start with all 1's
    m2 <- m1;
    pa <- rep(1,nr);
   ## Subset Sums
    SuffStat <-(estimand == "probabilities")*suff.stat(ModelMatrix, ObsTable/Total)+
              (estimand == "intensities")*suff.stat(ModelMatrix, ObsTable);
       
   ## Compute suff stat (subset sums) for m2:

    SuffStat_m2 <- suff.stat(ModelMatrix, m2);
   
   ## Check whether all subset sums are positive: 
    if(sum(SuffStat > 0) < nr) stop("Not all subset sums are positive. Some model parameters may be non-estimable.");
    
    SingleCells <- single.cells(ModelMatrix);
    idx  <- 0;
    exponent <- 1;
    while( max(abs(SuffStat_m2 - gamma*SuffStat)) > tol)
    {
      for ( j in 1:nr)  ## for each subset
      {
        for(i in 1:nc) ## look for all cels on the subset
        {
           if((j %in% SingleCells[,1])&&(i %in% SingleCells[,2]))
           {
              exponent <- 1;
           }
           else
           {
              exponent <- ModelMatrix[j,i];
           } 
           m2[i] <- m1[i]* 
                   (gamma*SuffStat[j]/(m1 %*% ModelMatrix[j,]))^exponent;
           if (ModelMatrix[j,i] > 0) {idx=i}
        }
        pa[j]  <-  pa[j] *     
        (m2[idx] / m1[idx])^(1/ModelMatrix[j,idx])

        m1 <- m2; 
      }
      SuffStat_m2 <- suff.stat(ModelMatrix, m2); 
     }
     result <- list( model.matrix = ModelMatrix,
                     observed.data = ObsTable,
                     fitted.values = (estimand =="probabilities")*Total*m2 + 
                    (estimand =="intensities")*m2, 
                     model.parameters = pa);   
      return(result)
}
