catminat_replace<-function(origin,change){
    
    origin<-as.character(origin)
    
    for(i in names(change)){
        origin<-gsub(i,change[i],origin)
    }
    return(origin)
}




catminat_download_to_local_directory<-function(directory){


    url<-"http://philippe.julve.pagesperso-orange.fr/baseflor.xlsx"
    ##baseflor<-read.xls (url, sheet = 1, header=T,method="tab")
    ##baseflor<-read.xls (url, sheet = 1, header=T,method="tab",fileEncoding="utf-8")
    temp_dest<-tempfile(fileext=".xlsx")
    download.file(url,temp_dest,mode="wb")

    catminat_df<-read_excel(temp_dest,sheet=1,col_names=T,col_types=rep("text",60))
    
    catminat_df<-catminat_df[grep("[0-9]+",catminat_df$rang_taxinomiqu,invert=T),]


    ##### TODO correggere nomi delle specie\


    ## remove entries for which CHOROLOGIE=="?"
    catminat_df<-catminat_df[catminat_df$CHOROLOGIE!="?",]
    
    ## change columns' name in order to
    ## avoid possible conflicts with non-
    ## ascii
    ## french chars
    recode_catminat_values<-c("Lumi.{1}re"="ell_L_fr",
                              "Temp.{1}rature"="elle_T_fr",
                              "Continentalit.{1}"="ell_C_fr",
                              "Humidit.{1}_.{1}daphique"="ell_U_fr",
                              "R.{1}action_du_sol_.pH."="ell_R_fr",
                              "Niveau_trophique"="ell_N_fr",
                              "Salinit.*"="ell_S_fr",
                              "Texture"="Soil_texture_fr",
                              "Mati.*re_organique"="organic_matter_fr",
                              "Humidit.*_atmosph.*rique"="ell_U_atm_fr",
                              "pollinisation"="poll_vect_fr",
                              "diss.*mination"="dissemination_fr",
                              "couleur_fleur"="flower_colour_fr",
                              "fruit"="fruit_type_fr",
                              "sexualit.*"="sex_reprod_fr",
                              "ordre_maturation"="order_of_maturation",
                              "inflorescence"="inflorescence_fr",
                              "Nom.Phytobase"="species_name",
                              "TYPE_BIOLOGIQUE"="li_form_fr"
                              )

    for(i in names(recode_catminat_values)){
        names(catminat_df)<-gsub(i,recode_catminat_values[i],names(catminat_df))
    }

    catminat_df$species_name<-as.character(catminat_df$species_name)
    catminat_df<-catminat_df[!is.na(catminat_df$species_name),]


    
    catminat_df$species_name<-gsub("&amp","&",catminat_df$species_name,perl=TRUE)
    catminat_df$species_name<-gsub(";","",catminat_df$species_name,perl=TRUE)
    catminat_df$species_name<-gsub("\\s+\\*$","",catminat_df$species_name,perl=TRUE)
    catminat_df$species_name<-gsub("\\s+A$","",catminat_df$species_name,perl=TRUE)
    catminat_df$species_name<-gsub("\\s+B$","",catminat_df$species_name,perl=TRUE)


    

    ## I split flowering dates into 2 columns,
    ## flower begin and end
    flowering<-as.character(catminat_df$floraison)
    beg_flow_fr<-as.numeric(gsub("^([0-9]+)+\\-([0-9])+$","\\1",flowering))
    end_flow_fr<-as.numeric(gsub("^([0-9]+)+\\-([0-9])+$","\\2",flowering))
    catminat_df<-data.frame(catminat_df,beg_flow_fr,end_flow_fr)

    
    ## I choose only a subset of the available columns in baseflor dataset
    selected_columns_catminat<-c(
        "species_name",
        "CHOROLOGIE",
        "inflorescence_fr",
        "sex_reprod_fr",
        "order_of_maturation",
        "poll_vect_fr",
        "fruit_type_fr",
        "dissemination_fr",
        "flower_colour_fr",
        "macule",
        "type_ligneux",
        "li_form_fr",
        "ell_L_fr",
        "elle_T_fr",
        "ell_C_fr",
        "ell_U_atm_fr",
        "ell_U_fr",
        "ell_R_fr",
        "ell_N_fr",
        "ell_S_fr",
        "Soil_texture_fr",
        "organic_matter_fr",
        "beg_flow_fr",
        "end_flow_fr",
        "PhytobaseID"
        )

    catminat_df<-catminat_df[,selected_columns_catminat]


    ## recode inflorescenses types
    inflorescences<-c(
        "capitule de capitules$"="Compound capitulum",
        "capitule simple$"="Capitulum",
        "c.*ne$"="Cone",                                
        "corymbe$"="Corymb",
        "corymbe de capitules$"="Corymb of capitula",
        "cyathe$"="Cyathium",
        "cyme bipare$"="Dichasia cyme",                         
        "cyme biscorpio.*de$"="Scorpioid cyme",
        "cyme capituliforme$"="Capituliform cyme ",
        "cyme d'.+pis$"="Cyme of spikes",
        "cyme d'ombelles$"="Cyme of umbels",                     
        "cyme de capitules$"="Cyme of capitula",
        "cyme de glom.+rules$"="Cyme of glomerula",
        "cyme multipare$"="Pleiochasium",
        "cyme unipare h.*lico.*de$"="Helicoid cyme",              
        "cyme unipare scorpio.*de$"="Scorpioid cyme",
        ".*pi d'.+pillets$"="Spike of spikelets",
        ".*pi de capitules$"="Spike of capitula",
        ".+pi de cymes triflores$"="Spike of three-flowers cymes",              
        ".+pi simple$"="Simple spike",
        "fleur solitaire lat.+rale$"="Solitary lateral flower",
        "fleur solitaire terminale$"="Solitary terminal flower",
        "glom.+rules$"="Glomerula",                          
        "glom.*rules spiciformes$"="Spike-like glomerula",
        "ombelle d'ombellules$"="Umbel of umbels",
        "ombelle simple$"="Simple umbel",
        "ombelle simple d'.+pis$"="Simple umbel of spikes",               
        "ombelle simple de capitules$"="Simple umbel of capitula",
        "panicule d'.+pillets$"="Panicle of spikelets",
        "panicule spiciforme$"="Spike-like panicle",
        "rac.{1}me capituliforme$"="Capitulum-like raceme",                 
        "rac.{1}me d'.+pis$"="Raceme of spikes",
        "rac.{1}me d'ombelles$"="Raceme of umbels",
        "rac.{1}me de capitules$"="Raceme of capitula",
        "ra.{1}me de cymes bipares$"="Raceme of dichasia cymes",             
        "rac.{1}me de cymes unipares h.{1}lico.{1}des$"="Raceme of helicoid cymes",
        "rac.{1}me de cymes unipares scorpio.{1}des$"="Raceme of scorpioid cymes",
        "rac.{1}me de rac.{1}mes$"="Raceme of racemes",
        "rac.{1}me simple$"="Simple raceme",                       
        "rac.{1}me de cymes bipares$"="Raceme of dichasia cymes", 
        "spadice$"="Spadix",
        "verticille d'ombelles$"="Verticil of umbels"               
        )
    catminat_df$inflorescence_fr<-catminat_replace(catminat_df$inflorescence_fr,inflorescences)

    ## recode fruit types
    fruit_types=c("ak.{1}ne$"="achene",
        "baie$"="berry",
        "capsule$"="capsule",
        "caryopse$"="caryopsis",
        "c.{1}ne$"="cone",
        "drupe$"="drupe",
        "follicule$"="follicle",
        "catminat_dfusse$"="legume",
        "pyxide$"="pyxid",
        "samare$"="samara",
        "silique$"="silique"  
        )
    catminat_df$fruit_type_fr<-catminat_replace(catminat_df$fruit_type_fr,fruit_types)

    ## recode flower colours
    flower_colours<-c("blanc$"="white",
                      "jaune$"="yellow",
                      "vert$"="green",
                      "marron$"="brown",
                      "bleu$"="blue",
                      "jaune$"="yellow",
                      "jauna$"="yellow",
                      "noir$"="black"
                      )
    catminat_df$flower_colour_fr<-catminat_replace(catminat_df$flower_colour_fr,flower_colours)

    ## recode dissemination types
    dissemination<-c(
        "an.*mochore$"="anemochores",
        "myrm.*cochore$"="myrmecochores",
        "myrm.*cochore$"="myrmecochores",
        "autochore$"="autochores",
        "barochore$"="barochores",
        "endozoochore$"="endozoochores",
        "endozoochorie$"="endozoochores",
        ".+pizoochore$"="epizoochores",              
        "dyszoochore$"="dyszoochores", 
        "hydrochore$"="hydrochores"                             
        )
    catminat_df$dissemination_fr<-catminat_replace(catminat_df$dissemination_fr,dissemination)

    ## recode sexual reproduction types
    sex_reprod<-c(
        "androdio.{1}que$"="Androdioecy",
        "gynomono.{1}que" ="Gynomonoecious",
        "gynodio.{1}que$"="Gynodioecious",
        "polygame$"="Polygamous",
        "mono.{1}que$"="Monoecious",
        "dio.{1}que$"="Dioecious",
        "hermaphrodite$"="Hermaphroditic",
        "polygame$"="Polygamous"
        )
    catminat_df$sex_reprod_fr<-catminat_replace(catminat_df$sex_reprod_fr,sex_reprod)

    ## remove species names where the phrase
    ## "sans nom" is found
    catminat_df<-catminat_df[grep("sans nom",catminat_df$species_name,invert=TRUE),]

    ## recode pollen vector
    poll_vec<-c(
        "an.{1}mogame"="wind",
        "autogame"="self",
        "apogame"="apogamy",
        "entomogame"="insect",
        "hydrogame"="water"
        )
    catminat_df$poll_vect_fr<-catminat_replace(as.character(catminat_df$poll_vect_fr),poll_vec)

    ## recode life_form_fr

    ## -- for the moment leave the original values
    
    
    ## remove entries without species names
    catminat_df<-catminat_df[catminat_df$species_name!="",]


    ## Remove double entries

    ## beware: catminat is now read with readxl package and empty cells are coded as <NA>. not as empty strings
    ##    catminat_df<-catminat_df[!(catminat_df$species_name=="Juncus articulatus" & catminat_df$sex_reprod_fr==""),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Juncus articulatus" & is.na(catminat_df$sex_reprod_fr)),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Capparis spinosa" & catminat_df$PhytobaseID=="13450"),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Centaurium ery. erythraea" & catminat_df$PhytobaseID=="2361"),]
##    catminat_df<-catminat_df[!(catminat_df$species_name=="Centaurium ery. erythraea" & catminat_df$PhytobaseID=="2361"),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Chenopodium ambrosioides" & is.na(catminat_df$flower_colour_fr)),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Chenopodium opulifolium" & catminat_df$CHOROLOGIE=="cosmopolite"),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Cornus sanguinea" & catminat_df$PhytobaseID=="2258"),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Crataegus mon. monogyna" & catminat_df$PhytobaseID=="1657"),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Daphne laureola" & catminat_df$PhytobaseID=="2084"),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Eryngium bourgatii" & catminat_df$PhytobaseID=="11913"),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Festuca ovi. ovina" & catminat_df$PhytobaseID=="352"),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Gaudinia fragilis" & catminat_df$li_form_fr=="test(hbis)"),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Genista salzmannii" & catminat_df$PhytobaseID=="12481"),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Medicago turbinata" & catminat_df$PhytobaseID=="12456"),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Ophrys bertolonii" & catminat_df$PhytobaseID=="725"),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Populus nigra" & catminat_df$PhytobaseID=="7038"),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Rhamnus saxatilis" & catminat_df$PhytobaseID=="16181"),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Rosmarinus officinalis" & catminat_df$PhytobaseID=="16080"),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Salix rosmarinifolia" & catminat_df$PhytobaseID=="795"),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Salvia officinalis" & catminat_df$PhytobaseID=="14526"),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Taraxacum sagittilobum" & catminat_df$PhytobaseID=="15326"),]
    catminat_df<-catminat_df[!(catminat_df$species_name=="Thlaspi rot. rotundifolium" & catminat_df$PhytobaseID=="1274"),]
    
    catminat_df<-catminat_df[!(catminat_df$species_name=="Linum bienne" & catminat_df$li_form_fr=="test(hbis)"),]

    ## remove double entries
    catminat_df<-catminat_df[!(duplicated(catminat_df)),]

    ## fix flowering dates for Medicago turbinata
    catminat_df[catminat_df$PhytobaseID==12456,c("beg_flow_fr","end_flow_fr")]<-c(4,6)
    ## and remove double entry
    catminat_df<-catminat_df[catminat_df$PhytobaseID!=11083,]

    ## fix CHOROLOGIE for Rubus cinerascens and  Rubus conspicuus and then remove double entries
    catminat_df<-catminat_df[(1:nrow(catminat_df))!=which(catminat_df$species_name=="Rubus cinerascens")[1],]
    catminat_df<-catminat_df[(1:nrow(catminat_df))!=which(catminat_df$species_name=="Rubus conspicuus")[1],]

    ##catminat_df<-catminat_df[row.names(catminat_df)!=(which(catminat_df$species_name=="Erodium glandulosum")[1]),]
    ##catminat_df<-catminat_df[row.names(catminat_df)!=(which(catminat_df$species_name=="Erodium rupicola")[1]),]
    ##catminat_df<-catminat_df[row.names(catminat_df)!=(which(catminat_df$species_name=="Onosma echioides")[1]),]

    ## remove CHOROLOGIE column
    catminat_df<-catminat_df[,names(catminat_df)!="CHOROLOGIE"]
    
    save(file=file.path(directory,"catminat.Rda"),catminat_df)
}

