retrieve_MycoFlor<-function(species,TRAITS,rest,data_myco){

    res<-new("results")
    if(is.null(TRAITS)){
        results<-NULL
    }else{
        temp_df<-as.data.frame(species,row.names = species)
        res_df<-data_myco[data_myco$Plant.species%in%species,]
        res_df<-merge(temp_df,res_df,by.x=0,by.y="Plant.species",all.x=TRUE)
        res_df$Plant.species<-as.character(res_df$species)
        temp_df<-with(res_df,tapply(Mycorrhizal.status,species,FUN = function(x){paste(unique(x),collapse = "-")}))
        temp_df<-as.data.frame(temp_df)

        names(temp_df)<-"MycoFlor"
        results<-temp_df

    }
    res@results<-results
    ## remove(list=c("myco"),pos =".GlobalEnv")
    stringa<-"Hempel, S., G\303\266tzenberger, L., K\303\274hn. I., Michalski, S.G.,\n Rillig, M.C., Zobel, M., and Moora, M., 2013. Mycorrhizas in the Central European flora:\n relationships with plant life history traits and ecology. Ecology 94: 1389-1399.\n"
    Encoding(stringa)<-"unicode"
    res@bibliography<-stringa
    return(res)
}
