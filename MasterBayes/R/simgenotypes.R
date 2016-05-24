"simgenotypes" <-function(A, E1=0, E2=0, ped, no_dup=1, prop.missing=0, marker.type="MSW"){

  id<-unique(ped[,1]) 
   
  pedN<-matrix(1, nrow=length(id), ncol=3)
  pedN[,3]<-as.numeric(match(ped[,3], id, nomatch=NA))
  pedN[,2]<-as.numeric(match(ped[,2], id, nomatch=NA))
  pedN[,1]<-1:length(id)

  if(length(names(A))==0 | (marker.type!="MSC" & marker.type!="MSW")){
    names(A)<-paste("l", 1:length(A), sep="")
  }
  if(length(names(A[[1]]))==0){
    alleles<-lapply(A, function(x){1:length(x)-1+as.numeric(marker.type=="MSC" || marker.type=="MSW")})
  }else{
    alleles<-lapply(A, function(x){names(x)})
  }

  true_genotypes<-lapply(names(A), function(x){x=list()})
  obs_genotypes<-lapply(names(A), function(x){x=list()})
  Gobs<-lapply(names(A), function(x){x=list()})

  base_ind<-which(is.na(pedN[,2])==TRUE & is.na(pedN[,3])==TRUE)
  descend_ind<-which(is.na(pedN[,2])==FALSE | is.na(pedN[,3])==FALSE)


for(l in 1:length(A)){

  true_genotypes[[l]]<-matrix(NA, length(pedN[,1]),2)
  true_genotypes[[l]][base_ind,]<-matrix(sample(alleles[[l]],length(base_ind)*2, replace=TRUE, prob=A[[l]]), length(base_ind),2)

  if(length(descend_ind)>0){
    for(off in descend_ind){

      dam_o<-pedN[,2][off]

      if(is.na(dam_o)==F){
        true_genotypes[[l]][,1][off]<-sample(true_genotypes[[l]][dam_o,], 1)
      }else{
        true_genotypes[[l]][,1][off]<-sample(alleles[[l]], 1, prob=A[[l]])
      }

      sire_o<-pedN[,3][off]

      if(is.na(sire_o)==F){
        true_genotypes[[l]][,2][off]<-sample(true_genotypes[[l]][sire_o,], 1)
      }else{
        true_genotypes[[l]][,2][off]<-sample(alleles[[l]], 1, prob=A[[l]])
      }
    }
  }

  if(no_dup>0){
    Gobs[[l]]<-matrix(NA, 0, 1+as.numeric(marker.type!="AFLP")) 
    for(i in 1:no_dup){

      if(marker.type=="MSW"){
        obs_genotypes[[l]]<-true_genotypes[[l]]
        hets<-which(obs_genotypes[[l]][,1]!=obs_genotypes[[l]][,2])
        drop_out<-which(rbinom(length(hets),1, prob=(2*E1*(1-E1))/(1-(E1^2)))==1) # these genotypes drop out
        if(length(drop_out)>0){
           which_allele<-rbinom(length(drop_out),1, prob=0.5)+1 # which allele drops out
           obs_genotypes[[l]][,1][hets[drop_out][which(which_allele==1)]]<-obs_genotypes[[l]][,2][hets[drop_out][which(which_allele==1)]]
           obs_genotypes[[l]][,2][hets[drop_out][which(which_allele==2)]]<-obs_genotypes[[l]][,1][hets[drop_out][which(which_allele==2)]]
        }
        stochastic<-which(rbinom(length(obs_genotypes[[l]]), 1, prob=E2)==1)
        if(length(stochastic)>0){
          for(s in stochastic){
            alt_al<-alleles[[l]][-which(alleles[[l]]==obs_genotypes[[l]][s])]
            if(length(alt_al)!=1){
              obs_genotypes[[l]][s]<-sample(alt_al, 1)
            }else{
              obs_genotypes[[l]][s]<-alt_al
            }
          }
        } 
      }
      if(marker.type=="MSC"){
         obs_genotypes[[l]]<-true_genotypes[[l]]
         stochastic<-which(rbinom(dim(obs_genotypes[[l]])[1], 1, prob=E2*(2-E2))==1)
         if(length(stochastic)>0){
           obs_genotypes[[l]][stochastic,]<-sample(names(A[[l]]), length(stochastic)*2, prob=A[[l]], replace=TRUE)
         }
      }
      if(marker.type=="AFLP"){
        obs_genotypes[[l]]<-true_genotypes[[l]]
        zeros<-which(obs_genotypes[[l]]==0)
        ones<-which(obs_genotypes[[l]]==1)
        ZtoO<-zeros[which(rbinom(length(zeros), 1, E2)==1)]
        OtoZ<-ones[which(rbinom(length(ones), 1, E1)==1)]
        obs_genotypes[[l]][ZtoO]<-1
        obs_genotypes[[l]][OtoZ]<-0
        obs_genotypes[[l]]<-as.matrix(apply(obs_genotypes[[l]], 1, function(x){as.numeric(1%in%x)}))
      }
      if(marker.type=="SNP"){
        obs_genotypes[[l]]<-true_genotypes[[l]]
        hets<-which(obs_genotypes[[l]][,1]!=obs_genotypes[[l]][,2])
        homs<-which(obs_genotypes[[l]][,1]!=obs_genotypes[[l]][,2])
        het2hom<-hets[which(rbinom(length(hets), 1, E1)==1)]
        hom2het<-homs[which(rbinom(length(homs), 1, E2)==1)]
        obs_genotypes[[l]][,1][het2hom]<-obs_genotypes[[l]][,2][het2hom]
        obs_genotypes[[l]][,1][hom2het]<-abs(obs_genotypes[[l]][,2][hom2het]-1)
      }
      Gobs[[l]]<-rbind(Gobs[[l]],obs_genotypes[[l]])
  }
  if(marker.type=="MSC" | marker.type=="SNP" | marker.type=="MSW"){
    Gobs[[l]]<-genotype(Gobs[[l]])
    true_genotypes[[l]]<-genotype(true_genotypes[[l]])
  }
  if(marker.type=="AFLP"){
    Gobs[[l]]<-genotypeD(Gobs[[l]])
    true_genotypes[[l]]<-genotype(true_genotypes[[l]], alleles=c("0","1"), reorder="no")
  }
  if(prop.missing!=0){
    Gobs[[l]][sample(1:length(Gobs[[l]]), floor(length(Gobs[[l]])*prop.missing))]<-NA
  }
}
}
names(true_genotypes)<-names(A)
if(no_dup>0){
names(Gobs)<-names(A)
}
list(G=true_genotypes, Gid=id[pedN[,1]], Gobs=Gobs, id=rep(id[pedN[,1]], no_dup))
}
