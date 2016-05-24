genGLCM <-
function(direction, distance, rawmat)
{
   #direction: 1=east(right)  2=south(down)  3=west(left)  4=north(up)
  
   #generate new GLCM-Matrix after calculation of size of GLCM 
   GLCM<-matrix(0, ncol=255, nrow=255) 
   (GLCM)

   #count occurrences and fill the matrix
   if(direction==1)   ### east-calculation
   {
     for(i in 1:dim(rawmat)[2]-distance)   #for all coloumns: -distance because of direction east 
     {
         for(a in 1:dim(rawmat)[1])        #for all rows
         {        
            GLCM[rawmat[a,i],rawmat[a,i+distance]]<-GLCM[rawmat[a,i],rawmat[a,i+distance]]+1
         } 
     }#eo for 
   }#eo east calculation 

   if(direction==2)   ### south-calculation
   {
     for(i in 1:dim(rawmat)[1]-distance)   #for all row: -distance because of direction south 
     {
         for(a in 1:dim(rawmat)[2])        #for all coloumns
         {        
            GLCM[rawmat[a,i+distance],rawmat[a,i]]<-GLCM[rawmat[a,i+distance],rawmat[a,i]]+1
         } 
     }#eo for 
   }#eo south calculation 

   #add the matrix to its transponse to make it symmetrical
   transGLCM<-t(GLCM)
   print("GLCM generation succesful")
   GLCM<-GLCM+transGLCM

   #normalize the matrix to turn it into probabilities
   GLCMprob<-round(GLCM/sum(GLCM), digits=4)   

   return(GLCMprob)
}
