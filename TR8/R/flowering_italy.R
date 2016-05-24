

luirig<-function(species){
    ## The website provides flowering periods with month expressed
    ## as roman numbers, thus a lookup table is needed for
    ## conversion purposes

    ##lookup_month<-data.frame(code=c(1:12,"NA"),roman=c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII"," "),stringsAsFactors = T)
    lookup_month<-data.frame(code=c(1:12),roman=c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII"),stringsAsFactors = T)
    base_url<-"http://luirig.altervista.org/flora/taxa/floraspecie.php?genere="
    ## a single page contains all the data for all the species belonging to the same genus

    ## This version adopt a different strategy: in luirig web pages there are href in the form "genus+species",
    ## thus these are used as target for our search: in the html structure, these nodes are precedeed by a <br>
    ## which contains the flowering dates which are extracted by means of regexp
    
    genus<-gsub("(^[a-zA-Z]+) .*","\\1",species,useBytes = TRUE)
    url<-paste(base_url,genus,sep = "")
    RES<-list()

    if(url.exists(url)){

        ## in the href species names are lowercase and species and genus are linked by a plus sign
        specie<-tolower(species)
        specie<-gsub("\\s+","+",specie)
        
        temp_pag<-htmlParse(url)
        ## search for the target node
        query=paste('//*[@href="index1.php?scientific-name=',specie,'"]',sep="")
        value=xpathApply(temp_pag,query,xmlValue)

        ## is our species found in the html parsed page?
        if(length(value)>0){
            ## extract our node
            uno<-getNodeSet(temp_pag,query)
            ## go back of 2 parents node (where the flowering data are found)
            due<-xmlParent(xmlParent(xmlParent(uno[[1]])))
            ## extract the text of the node
            vai<-xmlValue(due,"//*[text()[contains(., 'Fiorit')]]")
            ## extract flowering dates
            beg_fl<-gsub(".*Fiorit:([IVX]+)-([IVX]*)Tipo.*","\\1",vai)
            end_fl<-gsub(".*Fiorit:([IVX]+)-([IVX]*)Tipo.*","\\2",vai)        
            
            beg_fl<-mapvalues(beg_fl,lookup_month$roman,lookup_month$cod,warn_missing=FALSE)
            end_fl<-mapvalues(end_fl,lookup_month$roman,lookup_month$cod,warn_missing=FALSE)
            
            if(beg_fl%in%1:12&end_fl%in%1:12){
                RES[[species]]["IT_beg_flow"]=beg_fl
                RES[[species]]["IT_end_flow"]=end_fl
            }else{
                RES[[species]]["IT_beg_flow"]=NA
                RES[[species]]["IT_end_flow"]=NA
            }
        }else{
            RES[[species]]["IT_beg_flow"]=NA
            RES[[species]]["IT_end_flow"]=NA
        }
        }else{
            RES[[species]]["IT_beg_flow"]=NA
            RES[[species]]["IT_end_flow"]=NA
            
    }
    return(RES[[species]])
}


    


    

##' ##' get_italian_flowering get the beginning and the end of the flowering
##' phase for the italian flora. Values are based on Pignatti and retrieved
##' from the \samp{http://luirig.altervista.org/}
##'
##' 
##' 
##' @title get_italian_flowering 
##' @param species_list : a vector containing species names
##' @param TRAITS a vector containing the traits to be downloaded (used as a check for
##' tr8_gui() created variables)
##' @return a dataframe with two columns, the beginning and the end month (expressed as a number from 1 to 12)
##' of the flowering phase
##' @references \samp{http://luirig.altervista.org/}
##' @author Gionata Bocci <boccigionata@@gmail.com>
get_italian_flowering<-function(species_list,TRAITS,rest){
    res<-new("results")
    if(is.null(TRAITS)){
        res@results<-NULL
    }else{
        test<-("IT_beg_flow"%in%TRAITS||"IT_end_flow"%in%TRAITS)
        if(length(TRAITS)==0||test){
            ## names_month<-c("Gennaio","Febbraio","Marzo","Aprile","Maggio","Giugno","Luglio","Agosto","Settembre","Ottobre","Novembre","Dicembre")
            ## base_url="http://luirig.altervista.org/flora/taxa/index1.php?scientific-name="##vicia+sativa
            ## base_url_bis<-"http://luirig.altervista.org/schede/ae/"
            ## ##flower_date<-as.data.frame(names_month)
            flower_date<-list()
            for(species in species_list){
                Sys.sleep(rest)
                #species<-species_name
                #cur<-tolower(species_name)
                ## some species are found in url pattern base_url/genus+species
                ##cur<-gsub("^([A-Za-z]+) ([A-Za-z-]+).*","\\1+\\2",cur)
                ## other species are found in url pattern base_urlbis/genus_species.htm
                ## cur_bis<-gsub("^([A-Za-z]+) ([A-Za-z-]+).*","\\1_\\2",tolower(species))
                ## url<-paste(base_url,cur,sep="")
                ## url_bis<-paste(base_url_bis,cur_bis,".htm",sep="")
                ## check firt pattern
                

            ##     if(url.exists(url)){
            ##         ##flower_date<-merge(flower_date,luirig(url,species_name))
            ##             tp<-luirig(url,species_name)
            ##             for(i in 1:nrow(tp)){
            ##                 flower_date[[species]][[as.character(tp$names_month[i])]]<-tp[i,species]
            ##             }
            ##     }else{
            ##         for(i in names_month){
            ##                 flower_date[[species]][[as.character(i)]]<-NA
            ##             }
            ##     }
            ## }

                flower_date[[species]]<-luirig(species)
            }
            go<-ldply(flower_date)
            row.names(go)<-go$.id
            go<-go[,names(go)!=".id",drop=FALSE]
            go<-data.frame(go)
            res@results<-data.frame(go)
        }else{
            res@results<-NULL
        }
    }
    res@bibliography<-"Pignatti Sandro, 1982 Flora d'Italia.\nEdagricole, Bologna.\ndata retrieved from http://luirig.altervista.org/ "
    return(res)
}





