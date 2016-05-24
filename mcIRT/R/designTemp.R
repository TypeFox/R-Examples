designTemp <-
function(ngru,nit,TYPE="NLM")
  {
  if(TYPE=="NLM")
      {
        designL <- lapply(1:4,function(x)
                    {
                    matrix(1,nrow=ngru,ncol=nit)  
                    })
        
        names(designL) <- c("alpha","beta","zeta","lambda")
        
      }
  
  if(TYPE=="NRM")
      {
        designL <- lapply(1:2,function(x)
                    {
                      matrix(1,nrow=ngru,ncol=nit)  
                    })
        
        names(designL) <- c("zeta","lambda")
        
      }
  
  
  designL  
}
