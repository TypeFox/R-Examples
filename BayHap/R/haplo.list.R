`haplo.list` <-
function(nlocus,haplo1,haplo2,hetero,nhetero,rows,cols){

             
               h1mat<-rep(9,cols*rows)
               h2mat<-h1mat
               hetero2<-NULL
               l<-1
               for(l in 1:nlocus){

                  h1mat[(1)+(rows)*(l-1)]<-haplo1[l]

                  h2mat[(1)+(rows)*(l-1)]<-haplo2[l]

               }
           
               n<-1
               if(nhetero>=2){
                 hetero2<-hetero
                 while(hetero2[n]==0){

                    n=n+1

                 }
                 hetero2[n]<-0
                 fila<-1
                 q<-1
                 for(q in 1:nlocus){
                     if(hetero2[q]==1) {

                        a<-0

                        zero<-0

                        exponent<-nhetero-1

                        files<-2^exponent
                        i<-1

                        for(i in 1:files){

                           if(h1mat[(i)+rows*(zero)]!=9){

                              a<-a+1

                           }

                        }

                        p<-1
                        for(p in 1:a){

                            fila<-fila+1

                            if(q<(nlocus)){

                              b<-1
                              for(b in 1:nlocus){

                                 if(b==q){

                                    h1mat[(fila)+(rows)*(b-1)]<-h2mat[(p)+(rows)*(b-1)]

                                    h2mat[(fila)+(rows)*(b-1)]<-h1mat[(p)+(rows)*(b-1)]

                                 }else{

                                    h1mat[(fila)+(rows)*(b-1)]<-h1mat[(p)+(rows)*(b-1)];

                                    h2mat[(fila)+(rows)*(b-1)]<-h2mat[(p)+(rows)*(b-1)];

                                 }

                              }

                            }else{ 
                              r<-1

                              for(r in 1:nlocus){

                                  if(r==(nlocus)){

                                     h1mat[(fila)+(rows)*(r-1)]<-h2mat[(p)+(rows)*(r-1)];

                                     h2mat[(fila)+(rows)*(r-1)]<-h1mat[(p)+(rows)*(r-1)];

                                  }else{

                                     h1mat[(fila)+(rows)*(r-1)]<-h1mat[(p)+(rows)*(r-1)];

                                     h2mat[(fila)+(rows)*(r-1)]<-h2mat[(p)+(rows)*(r-1)];

                                  }

                              }

                           }

                        }
                      }
                  }               

              }
                  
return(list(h1mat=h1mat,h2mat=h2mat))
            
}

