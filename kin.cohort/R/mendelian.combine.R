`mendelian.combine` <-
function(p,inp,out){
#
# p:vector with allele prequencies
# inp: vector with respective number op genotypes in proband
#
three2two<- function(p,inp,out){
   q<- 1-p
   if(inp==2){ # dominant model
      tpm0<-mendelian(p)
      tpm2<-array(data=0, dim=c(2, 2, 3))
      for (k in 1:3) { # rel type
         for (i in 1:3) {# A
            tpm0[1,i,k]<-tpm0[1,i,k]* q^2 
            tpm0[2,i,k]<-tpm0[2,i,k]* 2*p*q 
            tpm0[3,i,k]<-tpm0[3,i,k]* p^2
         }
         tpm2[1,1,k]<- tpm0[1,1,k] / q^2
         tpm2[1,2,k]<-(tpm0[2,1,k] + tpm0[3,1,k]) / q^2
         tpm2[2,1,k]<-(tpm0[1,2,k] + tpm0[1,3,k]) / ( p^2 + 2*p*q )
         tpm2[2,2,k]<-(tpm0[2,2,k] + tpm0[3,2,k] + tpm0[2,3,k] + tpm0[3,3,k]) / ( p^2 + 2*p*q )
      }  
      tpm<-tpm2
   } else { # codominant model

      tpm<-mendelian(p)
      if(out==2){ # collapse output to carrier / non-carrier
         tpm2<-array(data=0, dim=c(3, out, 3))
         for (i in 1:3) {# proband geno
             for (k in 1:3) {  # rel type
                tpm2[i,1,k]<-tpm[i,1,k]                # aa
                tpm2[i,2,k]<-tpm[i,2,k] + tpm[i,3,k]   # Aa + AA           
             }
         }
         tpm<-tpm2
      }
   } 
tpm             
}

      if (length(p)==1)
         tpm<-three2two(p,inp[1],out[1])
   	else {    # 2 genotypes
         cp<-list(NULL,NULL)
         for (l in 1:2){
             if (inp[l]==2) out[l]<-2
             cp[[l]]<-three2two(p[l],inp[l],out[l])
         }
   
         ngeno.pro<-inp[1]*inp[2]
      	tpm = array(data=0, dim=c(ngeno.pro, out[1]*out[2], 3))
         gen1<- as.numeric(gl(inp[1],inp[2])) # 11 22 33
         gen2<- rep(1:inp[2],inp[1])          # 12 12 12
         ngeno.rel<- out[1]*out[2] # 4
         car1<- as.numeric(gl(out[1],out[2])) # 11 22
         car2<- rep(1:out[1],out[2])          # 12 12
         for (k in 1:3) {
            for (i in 1:ngeno.pro) {
                 for (j in 1:ngeno.rel) {
                      tpm[i,j,k] <- cp[[1]][gen1[i],car1[j] ,k]*cp[[2]][gen2[i],car2[j] ,k]
   
         } } }

   }
   tpm
}

