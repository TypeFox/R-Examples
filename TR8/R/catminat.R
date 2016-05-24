catminat<-function(species_list,TRAITS,catminat_df){
    res<-new("results")
    ## env<-new.env(parent = parent.frame())
    ## load(local_Catminat,envir = env)
    ## catminat_df<-get("catminat_df",envir=env)

    if(is.null(TRAITS)){
        res@results<-NULL
    }else{

        
        DF<-catminat_df[catminat_df$species_name%in%species_list,c("species_name",TRAITS)]
        row.names(DF)<-DF$species_name
        DF<-DF[,TRAITS,drop=FALSE]

        res@results<-DF
    }
    stringa<-"Julve, P., 1998 ff. - Baseflor. Index botanique, \303\251cologique et chorologique de la flore de France. Version : 26 November 2014 . http://perso.wanadoo.fr/philippe.julve/catminat.htm"
    Encoding(stringa)<-"unicode"
    res@bibliography<-stringa
    return(res)
}



