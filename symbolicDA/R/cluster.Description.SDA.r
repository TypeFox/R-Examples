cluster.Description.SDA<-function(table.Symbolic,clusters,precission=3){
#  attach(table.Symbolic)
  if(!.is.symbolic(table.Symbolic))stop("table.Symbolic should be symbolic data table")
  clusterNo<-max(clusters)
  resul<-as.list(1:clusterNo)
  for(i in 1:clusterNo){
    objectsInCluster=NULL
    for(j in 1:nrow(table.Symbolic$individuals)){
      if(clusters[j]==i) objectsInCluster<-c(objectsInCluster,j)
    }
    table.Symbolic$variablesDescription<-NULL
    for(k in 1:nrow(table.Symbolic$variables)){
      maxInterval<-NA
      minInterval<-NA
      categories<-NA
      minProbCategories<-NA
      maxProbCategories<-NA
      avgProbCategories<-NA
      sumProbCategories<-NA
      if(table.Symbolic$variables[k,"type"]=="IC"){
        maxInterval<-max(table.Symbolic$indivIC[objectsInCluster,k,2])
        minInterval<-min(table.Symbolic$indivIC[objectsInCluster,k,1])
      }
      if(table.Symbolic$variables[k,"type"]=="MN" || table.Symbolic$variables[k,"type"]=="N"){
         cat<-NULL
         for(j in 1:nrow(table.Symbolic$indivN)){
            if((sum(table.Symbolic$indivN[j,1]==objectsInCluster)!=0) && table.Symbolic$indivN[j,2]==k){
            cat<-c(cat,as.vector(table.Symbolic$indivN[j,3]))
            }
         }
         cat<-sort(unique(cat))
#         gcat<<-cat
#         if(k==3){
#          print(gcat)
#          stop()
#         }
         
         ccat<-NULL
         for(cc in cat){
            for(j in 1:nrow(table.Symbolic$detailsListNom)){
              if(table.Symbolic$detailsListNom[j,"details_no"]==as.integer(table.Symbolic$variables[k,"details"]) && table.Symbolic$detailsListNom[j,"num"]==cc){
                ccat<-c(ccat, as.character(table.Symbolic$detailsListNom[j,"label"]))
              }
            }
         }
         categories<-paste(ccat,collapse=";")
      }
      if(table.Symbolic$variables[k,"type"]=="NM"){
         cat<-NULL
         freq<-NULL
         for(j in 1:nrow(table.Symbolic$table.Symbolic$indivNM)){
            if((sum(table.Symbolic$table.Symbolic$indivNM[j,1]==objectsInCluster)!=0) && table.Symbolic$table.Symbolic$indivNM[j,2]==k){
              cat<-c(cat,table.Symbolic$table.Symbolic$indivNM[j,3])
              freq<-c(freq,table.Symbolic$table.Symbolic$indivNM[j,4])
            }
         }
         #print(k)
         #print(cat)
        minProbCategories<-""
        maxProbCategories<-""
        avgProbCategories<-""
        sumProbCategories<-""
         for(j in 1:max(cat)){
            #print(k)
            #print(cat)
            freqCat<-freq[cat==j]
            for(z in 1:nrow(table.Symbolic$table.Symbolic$detailsListNomModif)){
              if(table.Symbolic$table.Symbolic$detailsListNomModif[z,"details_no"]==as.matrix(table.Symbolic$variables)[k,"details"] && table.Symbolic$table.Symbolic$detailsListNomModif[z,"num"]==j){
                #print("znalazlem")
                catName<-as.character(table.Symbolic$table.Symbolic$detailsListNomModif[z,"label"])
              }
            }
            #print(freqCat) 
            #print(table.Symbolic$variables[k,"details"]))
            #print(catName)
            #if(k==2)stop("")
            minProbCategories<-paste(minProbCategories,round(min(freqCat),precission)," ",catName,";",sep="")
            maxProbCategories<-paste(maxProbCategories,round(max(freqCat),precission)," ",catName,";",sep="")
            avgProbCategories<-paste(avgProbCategories,round(mean(freqCat),precission)," ",catName,";",sep="")
            sumProbCategories<-paste(sumProbCategories,round(sum(freqCat),precission)," ",catName,";",sep="")
         }
         ccat<-NULL
         for(cc in cat){
         }
         categories<-paste(ccat,collapse=";")
      }
      vd<-c(as.character(table.Symbolic$variables[k,"name"]),as.character(table.Symbolic$variables[k,"label"]))
      if(any(as.matrix(table.Symbolic$variables)[,"type"]=="IC")){
        vd<-c(vd,c(minInterval,maxInterval))
      }
      if(any(as.matrix(table.Symbolic$variables)[,"type"]=="N") || any(as.matrix(table.Symbolic$variables)[,"type"]=="MN")){
        vd<-c(vd,c(categories))
      }
      if(any(as.matrix(table.Symbolic$variables)[,"type"]=="NM")){
        vd<-c(vd,c(minProbCategories,maxProbCategories,avgProbCategories,sumProbCategories))
      }
      table.Symbolic$variablesDescription<-rbind(table.Symbolic$variablesDescription,vd)
      dimnames(table.Symbolic$variablesDescription)[[1]]<-NULL
    }
    dn<-c("variable name","label")
    if(any(as.matrix(table.Symbolic$variables)[,"type"]=="IC")){
      dn<-c(dn,c("min value","max value"))
    }
    if(any(as.matrix(table.Symbolic$variables)[,"type"]=="N") || any(as.matrix(table.Symbolic$variables)[,"type"]=="MN")){
      dn<-c(dn,c("categories"))
    }
    if(any(as.matrix(table.Symbolic$variables)[,"type"]=="NM")){
      dn<-c(dn,c("min probabilities","max probabilities","avg probabilities","sum probabilities"))
    }
    print(table.Symbolic$variablesDescription)
    print(dn)
    dimnames(table.Symbolic$variablesDescription)[[2]]<-dn
    resul[[i]]<-table.Symbolic$variablesDescription
  }
#  detach(table.Symbolic)
  resul
}
