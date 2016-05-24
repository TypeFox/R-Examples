insertPed<-function(ped, founders=NULL){

   ped[,1]<-as.character(ped[,1])
   ped[,2]<-as.character(ped[,2])
   ped[,3]<-as.character(ped[,3])

   mmothers<-na.omit(ped[,2][which(ped[,2]%in%ped[,1]==FALSE)])
   mfathers<-na.omit(ped[,3][which(ped[,3]%in%ped[,1]==FALSE)])
   if(is.null(founders)==FALSE){
     founders<-na.omit(founders[which(founders%in%ped[,1]==FALSE)])
   }
   mparents<-unique(c(mmothers, mfathers, founders))

   nped<-ped[rep(1,length(mparents)),]
   nped[,1]<-mparents
   nped[,2]<-NA
   nped[,3]<-NA

   nped<-rbind(nped,ped)
   colnames(nped)<-colnames(ped)
   nped
}

