PLANTS <- function(species_list,TRAITS){

    
    res<-new("results")
    if(is.null(TRAITS)){
        res@results<-NULL
    }else{


        env<-new.env(parent = parent.frame())
        PLANTS<-get("PLANTS",envir=env)
        species<-data.frame(species=species_list)
        
        ##temp_df<-merge(species,PLANTS,by.x="species",by.y="acceptedname",all.x=T)
        temp_df<-merge(species,PLANTS,by.x="species",by.y="matchedname",all.x=T)
        row.names(temp_df)<-temp_df$species
        temp_df<-temp_df[,!names(temp_df)%in%c("species"),drop=F]
        temp_df<-temp_df[,TRAITS,drop=F]

                                        #   remove(list=c("BROT"),envir = env)
        res@results<-temp_df
    }
    stringa<-"Green, W. (2009) USDA PLANTS Compilation, version 1, 09-02-02."
    Encoding(stringa)<-"unicode"
    res@bibliography<-stringa

    return(res)

}
