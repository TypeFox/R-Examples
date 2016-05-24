fake_row_make<-function(modsA,modsB,n.modsA,n.modsB,is.exa=T){
# modsA is the list containing the modalities of each variable
# n.modsA is the number of modalities of each variable  
fake.row=matrix(1,1,length(n.modsA))

catDiff=(n.modsA-n.modsB)
chkDiff=which((n.modsA-n.modsB)!=0)

zero.cat=list()
for(e in 1:length(catDiff)){
  zero.cat[[e]]=setdiff(modsA[[e]],modsB[[e]])
  }

#print(zero.cat[[chkDiff[2]]])

for(f in 1:length(chkDiff)){
    for(g in 1:length(zero.cat[[chkDiff[f]]])){
      #print(length(zero.cat[[chkDiff[f]]]))
      if(is.exa==T){
        frow=matrix(1,1,length(n.modsA)) 
      }
      else{
        frow=matrix(0,1,length(n.modsA))
      }
        frow[chkDiff[f]]=zero.cat[[chkDiff[f]]][g]
        fake.row=rbind(fake.row,frow)
        }
    }
fake.row=fake.row[-1,]
fake.row
}