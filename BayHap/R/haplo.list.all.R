`haplo.list.all` <-
function(nlocus){

               nhetero<-nlocus
               haplo1<-rep(1,nlocus)
               haplo2<-rep(0,nlocus)
               hetero<-rep(1,nlocus)
               rows<-2^(nlocus-1)

               cols<-nlocus
               H<-NULL
              
              

               llista<-haplo.list(nlocus,haplo1,haplo2,hetero,nhetero,rows,cols)
               Hmat1<-llista[[1]]
               Hmat2<-llista[[2]]
               exponent<-nhetero-1

               files<-2^exponent
               H<-matrix(rep(9,2*nlocus*files),nrow=2*files)
              
               i<-1
               j<-1

               for(i in 1:files){

                   for(j in 1:nlocus){
                       H[2*i-1,j]<-Hmat1[i+rows*(j-1)]

                       H[2*i,j]<-Hmat2[i+rows*(j-1)]

                   }

               }

return(H)

}

