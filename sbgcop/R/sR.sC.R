"sR.sC" <-
function(sC) {
p<-dim(sC)[1]
s<-dim(sC)[3]
sR<-array(dim=c(p,p,s) )
dimnames(sR)<-dimnames(sC)

for(l in 1:s) {

C<-sC[,,l]
R<-C*NA
for(j in 1:p){
R[j,-j]<- C[j,-j]%*%solve(C[-j,-j])
             }
sR[,,l]<-R     
              }
sR                    }

