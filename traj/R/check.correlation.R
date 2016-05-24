check.correlation <-
function(output, verbose = TRUE, is.return = FALSE)
{
  cor.mat = cor(output)

  mes.names = names(output)
  
  is.corr = FALSE
  
  corr.var = NULL
  
  for(i_row in mes.names[-length(mes.names)]){

    i_pos = which( mes.names == i_row)
    res.names = mes.names[(i_pos + 1) : length(mes.names)]
    
    for(i_col in res.names){

      if(cor.mat[i_row, i_col] > 0.999){
        
        corr.var = rbind(corr.var, c(i_col, i_row))
        
        if(verbose){ 
          print(paste("Correlation of ",i_row, " and ", i_col, " : ", round(cor.mat[i_row, i_col],3), sep = ""))
          is.corr = TRUE
        }
      }
    }
  }
  if(!is.corr && verbose)
    print("No correlations found. That is good.")
  if(is.return)
    return(corr.var)
}
