
codes_column<-c(
    "High elevation"="high_el"
   ,"Low elevation"="low_el"
   ,"Seed mass"="seed_mas_cal"
   ,"Annual seed production per plant (c)"="seed_pro_cal"
   ,"Blade area (c)"="blade_area"
   ,"Leaf area to sapwood area"="leaf_to_sap_area"
   ,"Leaf Narea"="leaf_N_area"
   ,"Leaf Nmass"="leaf_N_mass"
   ,"Maximum height (c)"="max_height_cal"
   ,"Specific leaf area"="sla_cal"
   ,"Wood density"="wood_dens"
   ,"Blade length (c)"="blade_length"
   ,"Blade width (c)"="blade_width"
   ,"Leaf thickness"="leaf_thick"
   ,"Leaf type"="leaf_type"
)


eflora<-function(species_list,TRAITS){

    res<-new("results")
    if(is.null(TRAITS)){
        res@results<-NULL
    }else{

        ## query of Electronic Flora of California is done at the genus level
        ## in the form of
        ## http://ucjeps.berkeley.edu/cgi-bin/DT.pl?wanted_genus=Acer&wanted_genus=Achillea
        ## thus I extract the genuses of the requested species and then add "&wanted_genus="
        ## before each species
        genus_list<-sapply(species_list,function(x){gsub("^(\\w+)\\s+.*","\\1",x,perl=T)})
        code_url<-paste(paste("wanted_genus=",genus_list,sep=""),collapse="&")
        ## then I paste the string to the base url to have the complete URL for the query
        base_url<-"http://ucjeps.berkeley.edu/cgi-bin/DT.pl?"
        long_url<-paste(base_url,code_url,sep="")
        

        lista<-list()
        ## read the data
        berk<-readHTMLTable(long_url)[[2]]
        ## remove units of measure from the data values column
        berk$V3<-gsub("\\t.*$","",berk$V3)
        berk<-berk[,c("V1","V2","V3")]
        ## recode the trait names in the short form used in the gui
        berk$V2<-revalue(berk$V2,codes_column,warn_missing=FALSE)

        ## create a species*trait dataframe
        ## dati<-cast(V1~V2,value="V3",data=berk,fun.aggregate=function(x){paste(x,sep="-")})
        berk<-cast(V1~V2,value="V3",data=berk)
        ## add empty columns for each traits which was requested but for which no data
        ## were found in the traitbase
        nomi<-TRAITS[!TRAITS%in%names(berk)]
        for(i in nomi){
            berk[[i]]<-"<NA>"
        }
        
        ## create a dataframe with the requested species' names
        TP<-data.frame(V1=species_list)
        berk<-merge(TP,berk,all.x=TRUE)
        row.names(berk)<-berk$V1
        berk<-berk[,TRAITS,drop=FALSE]
        berk<-droplevels.data.frame(berk)
        res@results<-berk
    }

    stringa<-paste("Jepson Flora Project, 2006. Ecological Flora of California, 23 July 2006, ",format(Sys.time(), "%b-%d-%Y"),sep="")
    Encoding(stringa)<-"unicode"
    res@bibliography<-stringa

    return(res)


}
