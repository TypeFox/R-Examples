single.cells <-
function(ModelMatrix)
{
   setR <- NULL;
   setC <- NULL;
   all.single.cells <- NULL;
   nr <- nrow(ModelMatrix);
   nc <- ncol(ModelMatrix);
   for(j in 1:nr)
   {
      for(i in 1:nc)
      {
         if((sum(ModelMatrix[j,]>0)==1)&&(sum(ModelMatrix[,i]>0)==1)&&
            (ModelMatrix[j,i] !=0)) 
         {
            all.single.cells <- rbind(all.single.cells, c(j,i))
         }
      }
    }
   return (all.single.cells)
}
