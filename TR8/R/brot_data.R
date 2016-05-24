brot_data<-function(species_list,TRAITS){
        
    res<-new("results")

        if(is.null(TRAITS)){
        res@results<-NULL
    }else{
            
            env<-new.env(parent = parent.frame())
            BROT<-get("BROT",envir=env)
            species<-data.frame(species=species_list)
    
            temp_df<-merge(species,BROT,by.x="species",by.y="taxa",all.x=T)
            row.names(temp_df)<-temp_df$species
            temp_df<-temp_df[,!names(temp_df)%in%c("species"),drop=F]
            temp_df<-temp_df[,TRAITS,drop=F]
            
            ##   remove(list=c("BROT"),envir = env)
            res@results<-temp_df
        }
    stringa<-"Paula S, Arianoutsou M, Kazanis D, Tavsanoglu \303\207, Lloret F, Buhk C, Ojeda\n F, Luna B, Moreno JM, Rodrigo A, Espelta JM, Palacio S, Fern\303\241ndez-Santos\n B, Fernandes PM, and Pausas JG. 2009.\n Fire-related traits for plant species of the Mediterranean Basin. Ecology 90: 1420. \n AND \n Paula S. & Pausas J.G. 2013. BROT: a plant trait database for\n Mediterranean Basin species.\n Version 2013.06. URL: http://www.uv.es/jgpausas/brot.htm"
    Encoding(stringa)<-"unicode"
    res@bibliography<-stringa

    return(res)
}


brot_aggregate<-function(x){

    ## check if numeric data are available in the
    ## sub-vector passed to aggregate
    numeric_data<-x[grep("[0-9]+\\.[0-9]+",x)]

    ## if true, then use those numeric data
    ## instead of the qualitative ones
    
    if(length(numeric_data)>0){

        ## some data are in the form "0.6-0.9", thus we need to
        ## split the two numbers ...
        
        numeric_data<-strsplit(numeric_data,"-")

        ## .. convert both from string to numbers and
        ## calculate their mean

        numeric_data<-unlist(lapply(numeric_data,function(Z){mean(as.numeric(Z))}))
        res<-mean(numeric_data)

    }else{

        ## if we have only categorical variable, I propose to take only
        ## the most frequent levels
        tavola<-table(x)
        res<-names(which(tavola==max(tavola)))
        res<-paste(unlist(res),collapse=", ")

        res<-unlist(res)
    }

    return(unlist(res))
}




brot_download_to_local_directory<-function(directory){

    BROT<-read.delim("http://www.uv.es/jgpausas/brot/BROT_2013.06_data.txt",sep="\t",header=T)
    BROT$Data<-as.character(BROT$Data)

    BROT<-aggregate(Data~taxa+Trait,BROT,FUN=brot_aggregate)
    BROT$taxa<-as.character(BROT$taxa)
    BROT$Trait<-as.character(BROT$Trait)

    temp_BROT<-reshape::cast(data=BROT,taxa~Trait,value="Data",fun.aggregate=function(x){paste(x,collapse=", ")})
    temp_BROT<-as.data.frame(temp_BROT)
    temp_BROT[]<-lapply(temp_BROT,as.character)
    temp_BROT[temp_BROT==""]<-NA
    BROT <- temp_BROT
    save(file=file.path(directory,"BROT.Rda"),BROT)
    ##    return(temp_BROT)
}

