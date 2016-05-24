PePNRM <-
function(quads,mueERG=mueERG)
{
  
  thetas <- mapply(function(eachG,ql)
  {
    PPzae  <- ql$nodes * eachG$LjmalA
    thet <- colSums(PPzae) / eachG$Pji_schlange
    thet
  },eachG=mueERG$mue_hat_g,ql=quads,SIMPLIFY=FALSE)
  
  
  return(thetas=thetas) 

}
