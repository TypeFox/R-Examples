discordance.delta <-
function(x, names = NULL, theta = 0.05, corrected = TRUE, 
                              printToTex = FALSE, directory = NULL, file.name = NULL){
## internal functions ##  
  delta.prime <-
    function(
      q,
      theta
    ) 
    {
      sum(q<theta) / length(q)
    }
  
  delta <-
    function(
      rA,
      rB,
      q,
      theta
    ) 
    {
      sum((rA+rB)*(q<theta)) / (sum(rA) + sum(rB))
    }
##############
  
  
  if(is.null(names)){
    samples <- names(x)
  }
  
  n.s <- length(x)
  out.df <- data.frame(samples = samples, 
                       Delta = 1:n.s, 
                       Delta.prime = 1:n.s)
  
  theta = theta
  
  if (corrected == T){
    ind <- 7
  }else{
    ind <- 3
  }
  
  for (i in 1:n.s) {
    in.df <- as.data.frame(x[[i]])
    out.df$Delta[out.df$samples == samples[i]] <- 
      delta(
        in.df$freqA,
        in.df$freqB,
        in.df[ind],
        theta
      )
    out.df$Delta.prime[out.df$samples == samples[[i]]] <- 
      delta.prime(
        in.df[[ind]],
        theta
      )
  }
  if(printToTex ==T){
    if(is.null(file.name)){
      file.name = paste("discordanceDelta_", Sys.Date(), ".txt", sep = "" )
    }
    if(!is.null(directory)){
      file.name = paste(directory, "/", file.name, sep = "")
    }
    
    table = xtable(out.df, label="tab:Delta", digits = c(3), caption = "$\\Delta$ and $\\Delta'$.")
    print(table, file = file.name)
  }
  return(out.df)
}
