##
## Class Tr8 is used as a "containter" for all other functions and classes
## Class Tr8
## needed to download traits data from various databases
## 
## @name Tr8
## @rdname Tr8-Class
## @exportClass Tr8
## @slot species_list a list of species for which traits data are to be searched
## @slot results dataframe containing scraped traits
## @slot not_valid species whose name were not present in the Ecoflora database
## @slot double_names species for which more than one name was found
setClass("Tr8",representation =list(species_list="vector",results="data.frame",not_valid="vector",double_names="vector",bibliography="list",reference="data.frame",issues="ANY"))

## Method issues
##
## @name issues
## @rdname Tr8-Class
## @exportMethod issues
setGeneric(name="issues",def=function(.Object){standardGeneric("issues")})

## @rdname Tr8-Class
## @param .Object an object of class Tr8
setMethod(f="issues",
          signature = "Tr8",
          definition = function(.Object){
    check=FALSE
    ## warning for double names
    if(length(.Object@double_names)>0){
        cat("\n")
        cat("\t WARNING\n")
        cat("\n")
        for(sp in .Object@double_names){
            cat(paste("\tFor species",sp," multiple matched names were found\n"))
        }
        cat("\n")
        check=TRUE
    }
    ## warning for missing species
    if(length(.Object@not_valid)>0){
        cat("\n")
        cat("\t WARNING\n")
        cat("\n")
        for(sp in .Object@not_valid){
            cat(paste("\tFor species",sp," no matched names were found\n"))
        }
        cat("\n")
        check=TRUE
    }
    ## advice about potential issues given by Ecoflora
    if(!is.null(.Object@issues)){
        cat(.Object@issues)
        check=TRUE
    }
    if(check){
        cat("\tPlease check that these results are consistent with your orginal dataset!\n")
        ##              return(.Object)
    }else{

        cat("No particular problems were faced in the data retrieval process.\n")

    }
    
}

)


## Method lookup
##
## @name lookup
## @rdname Tr8-Class
## @exportMethod lookup
setGeneric(name="lookup",def=function(.Object){standardGeneric("lookup")})

## @rdname Tr8-Class
## @aliases lookup, Tr8-Class
## @param .Object an object of class Tr8
setMethod(f="lookup",
          signature="Tr8",
          definition = function(.Object){
    REF<-.Object@reference
    RES<-.Object@results
    DF<-REF[REF$short_code%in%names(RES),]

    cat("\n")
    cat("\n")
    cat("*****************************************************************")
    cat("\n")
    cat("To interpret the traits data, please refer to the following table\n")
    cat("\n")
    cat(sprintf("%-30s\t%-40s\t%-40s\n"," code","description","reference database\n"))
    cat(sprintf("%-30s\t%-40s\t%-40s\n"," ----","-----------","------------------\n"))
    for(i in 1:nrow(DF)){
        cat(sprintf("%-30s\t%-40s\t%-30s\n",DF[i,2],DF[i,3],DF[i,4]))
    }
    cat("\n")
    cat(sprintf("%-30s\t%-40s\t%-40s\n"," ****","***********","******************\n"))
    ##tp<-.Object@reference
    ##tp<-tp[,c("short_code","description","db")]
    tp<-DF[,c("short_code","description","db")]
    names(tp)<-revalue(names(tp),c("short_code"="code","db"="reference database"))
    return(invisible(tp))
}
)



## @rdname Tr8-Class
## @aliases show, Tr8-Class
## @param object an object of class Tr8
setMethod(f="show",
          signature="Tr8",
          function(object){
    print(object@results)
}
)


## A method to extract results form a Tr8 object
##
## the slots @results will be returned by the function
setGeneric("extract_traits",def=function(object){standardGeneric("extract_traits")})

setMethod(f="extract_traits",
          signature="Tr8",
          function(object){
    return(object@results)
}
)


## Method bib
                                        #@name bib
                                        #@rdname Tr8-Class
                                        #@exportMethod bib
setGeneric(name="bib",def=function(.Object){standardGeneric("bib")})


setMethod(f="bib",
          signature="Tr8",
          definition = function(.Object){
    env<-new.env(parent = parent.frame())
    data(column_list,envir = env)
    column_list<-get("column_list",envir=env)

    cat("\n")
    cat("Please use the following references for the data you retrieved with tr8()\n")
    cat("\n")
    for(db in names(.Object@bibliography)){
        cat("************************************************\n")
        cat("\n\nFor the following traits:\n\n")
        for(trait in .Object@bibliography[[db]]){
            cat("\t * ",paste(trait),"\n")
        }
        cat("\nplease use:\n\n")
        cat(db,fill=TRUE)
        cat("\n")
    }
    cat("************************************************\n")
}

)

## @rdname Tr8-Class
## @aliases bib, Tr8-Class
## @param .Object an object of class Tr8
## setMethod(f="bib",
##           signature="Tr8",
##           definition = function(.Object){
##               cat("\n")
##               cat("Please use the following references when using data retrieved with tr8()\n")
##               cat("\n")
##               ## cat(sep="\n",strwrap("BIOLFLOR - Eine Datenbank zu biologisch-ökologischen Merkmalen zur Flora von Deutschland. Schriftenreihe für Vegetationskunde 38: 1-333.  (Bundesamt für. Bonn, Bundesamt für Naturschutz)"))
##               ## cat("\n")
##               ## cat(sep="\n",strwrap("Kleyer, M., Bekker, R.M., Knevel, I.C., Bakker, J.P, Thompson, K., Sonnenschein, M., Poschlod, P., Van Groenendael, J.M., Klimes, L., Klimesova, J., Klotz, S., Rusch, G.M., Hermy, M., Adriaens, D., Boedeltje, G., Bossuyt, B., Dannemann, A., Endels, P., G\xf6tzenberger, L., Hodgson, J.G., Jackel, A-K., Kuehn, I., Kunzmann, D., Ozinga, W.A., Römermann, C., Stadler, M., Schlegelmilch, J., Steendam, H.J., Tackenberg, O., Wilmann, B., Cornelissen, J.H.C., Eriksson, O., Garnier, E., Peco, B. (2008): The LEDA Traitbase: A database of life-history traits of Northwest European flora. Journal of Ecology 96: 1266-1274."))
##               ## cat("\n")
##               ## cat(sep="\n",strwrap("Fitter, A . H. and Peat , H. J., 1994, The Ecological Flora Database, J. Ecol., 82, 415-425."))
##               ## cat("\n")
##               ## cat(sep="\n",strwrap("Pignatti S., Menegoni P., Pietrosanti S., 2005, Biondicazione attraverso le piante vascolari. Valori di indicazione secondo Ellenberg (Zeigerwerte) per le specie della Flora d'Italia. Braun-Blanquetia 39, Camerino, pp.  97."))
##               ## cat("\n")
##               ## cat(sep="\n",strwrap("Asem A. Akhmetzhanova, Nadejda A. Soudzilovskaia, Vladimir G. Onipchenko, Will K. Cornwell, Vladimir A. Agafonov, Ivan A. Selivanov, and Johannes H. C. Cornelissen. 2012. A rediscovered treasure: mycorrhizal intensity database for 3000 vascular plants species across the former Soviet Union. Ecology 93:689. URL: http://esapubs.org/Archive/ecol/E093/059/default.htm"))
##               ## cat("\n")
##           }
##           )


#' \code{tr8}: a function for retrieving functional traits data from various
#' databases. 
#' 
#' \code{tr8} makes use of other function provided by the \code{TR8} package in
#' order to query various databases and provide the user with a dataframe
#' containing traits data for the species of interest. 
#' 
#' @param species_list a vector containing names of the plant species for which
#' traits data want to be extracted.
#' @param download_list a 
#' @param gui_config if set to TRUE a GUI for selecting traits of interest is shown (default is TRUE)
#' @return data.frame containing various traits data for the species of interest
#' @author Gionata Bocci <boccigionata@@gmail.com>
#' @seealso \code{\link{available_traits}}, \code{\link{ecoflora}}, \code{\link{leda}}, \code{\link{biolflor}},\code{\link{pignatti_f}}
#' @references Please always use the following citations any time you use trait
#' data retrieved with \code{tr8}
#' 
#' \bold{BiolFlor}
#'
#' BIOLFLOR - Eine Datenbank zu biologisch-ökologischen Merkmalen zur Flora von Deutschland. Schriftenreihe für Vegetationskunde 38: 1-333.  (Bundesamt für. Bonn, Bundesamt für Naturschutz)
#' 
#' \bold{Ecoflora}
#' 
#' Fitter, A . H. and Peat , H. J., 1994, The Ecological Flora Database, J.
#' Ecol., 82, 415-425.  \samp{http://www.ecoflora.co.uk}
#' 
#' \bold{LEDA traitbase}
#' Kleyer, M., Bekker, R.M., Knevel, I.C., Bakker, J.P, Thompson, K., Sonnenschein, M., Poschlod, P.,
#' Van Groenendael, J.M., Klimes, L., Klimesova, J., Klotz, S., Rusch, G.M., Hermy, M., Adriaens, D.,
#' Boedeltje, G., Bossuyt, B., Dannemann, A., Endels, P., Götzenberger, L., Hodgson, J.G., Jackel, A-K.,
#' Kühn, I., Kunzmann, D., Ozinga, W.A., Römermann, C., Stadler, M., Schlegelmilch, J., Steendam, H.J.,
#' Tackenberg, O., Wilmann, B., Cornelissen, J.H.C., Eriksson, O., Garnier, E., Peco, B. (2008):
#' The LEDA Traitbase: A database of life-history traits of Northwest European flora.
#' Journal of Ecology 96: 1266-1274.
#'
#' \bold{Akhmetzhanova et al, 2012}
#' 
#' Asem A. Akhmetzhanova, Nadejda A. Soudzilovskaia, Vladimir G. Onipchenko,
#' Will K. Cornwell, Vladimir A. Agafonov, Ivan A. Selivanov, and Johannes H. C. Cornelissen. 2012.
#' A rediscovered treasure: mycorrhizal intensity database for 3000 vascular plants
#' species across the former Soviet Union. Ecology 93:689.
#' URL: http://esapubs.org/Archive/ecol/E093/059/default.htm
#'
#' 
#' \bold{Pignatti et al., 2005}
#' 
#' Pignatti S., Menegoni P., Pietrosanti S., 2005, Biondicazione attraverso le piante vascolari.
#' Valori di indicazione secondo Ellenberg (Zeigerwerte) per le specie della Flora d'Italia.
#' Braun-Blanquetia 39, Camerino, pp.  97.
#'
#' #' @examples \dontrun{
#' #My_traits<-tr8(species_list=c("Abies alba"),download_traits=c("le_area","h_max","h_min"))
#' }
#' @export tr8
tr8<-function(species_list,download_list=NULL,gui_config=FALSE,synonyms=FALSE){

    ## if(tryCatch(nsl("www.cran.r-project.org"), error =function(e){return(FALSE)},warning=function(w){return(FALSE)})==FALSE){
    ##     stop("You need a working internet connection to use tr8()")
    ## }
    
    ## get column_list dataset
    env<-new.env(parent = parent.frame())
    data(column_list,envir = env)
    column_list<-get("column_list",envir=env)
    ## load lookup table
    ## convert it to a data frame
    temp_dframe<-ldply(column_list)
    names(temp_dframe)<-c("long_code","short_code","description","db")

    op<-options()
    options("guiToolkit"="tcltk")
    ## rest is used for Sys.sleep in all the retrieving functions
    rest=0.01

    ## dir.create does not seem to work under windows, thus
    ## TR8 will not try to create its own subdirectory any more
    ## but will simply use the standard user_data_dir
    ## appname <- "TR8"
    ## appauthor <- "GioBo"
    ## directory<-user_data_dir(appname, appauthor)
    directory<-user_data_dir()
    

    
    if(missing(species_list)||!is.character(species_list)){
        message("\ntr8() accepts only a list of plant species names \nplease read help(tr8)\n")
    }else{
        traits_list<-list()
        ## if the user wants to manually sets the parameters to download
        if(gui_config)
        {
            ##gmessage(title="TR8 reminder!","Please always use the appropriate citations for the downloaded data.\n
            ##\n Run the bib() function on the downloaded data to get the correct bibliographic citations to be used.\n")
            ## run the gui
            traits_list<-tr8_config()
        }else{
            for(db in c("BiolFlor","LEDA","Ecoflora","Pignatti","AMF","Catminat","BROT","PLANTS","EFlora_Cal")){
                                        #db<-temp_dframe$db[temp_dframe$short_code==i]
                data_db<-temp_dframe[temp_dframe$db==db,]
                if(sum(data_db$short_code%in%download_list)>0){
                    code<-data_db$long_code[data_db$short_code%in%download_list]
                }else{code<-NULL}
                traits_list[db]<-list(code)
            }
        }
        if(synonyms==TRUE){
            
            check_names<-tnrs(species_list)
            check_names<-check_names[,c("submittedname","acceptedname","matchedname")]
            
            reference_names<-lapply(species_list,function(x){
                
                sp_names<-check_names[check_names$submittedname==x,]
                sp_names<-unique(unlist(sp_names))
                sp_names<-sp_names[grep("^\\w+ \\w+.*$",sp_names)]
                return(sp_names)

            }
            )
            names(reference_names)<-species_list
            species_list<-unique(as.vector(unlist(reference_names)))
            
        }
        
        
        ## retrieve traits from ecolora function
        local_ecoflora<-file.path(directory,"ECOFLORA_df.Rda")
        if(file.exists(local_ecoflora)){
            load(local_ecoflora)}else{
            if(length(traits_list$Ecoflora)>0){
                local_storage(db="Ecoflora",directory)
                load(local_ecoflora)
            }
        }
        eco_traits<-ecoflora(species_list,TRAITS=traits_list$Ecoflora,rest=rest)

        ## check if an already downloaded version of the LEDA database
        ## exists and, if so, use it otherwise download a copy, but only
        ## if at least one LEDA trait is needed
        local_leda<-file.path(directory,"leda_database.Rda")
        if(file.exists(local_leda)){
            load(local_leda)}else{
            if(length(traits_list$LEDA)>0){

                ## unfortunately nls() does not work on Windows, thus I think it's better to remove that
                ## if(tryCatch(nsl("www.cran.r-project.org"), error =function(e){return(FALSE)},warning=function(w){return(FALSE)})==FALSE){
                ##     stop("You neither have a working internet connection nor locally stored LEDA files.\n  Please re-run tr8() function when your internet connection is working.")
                ## }
                url_leda="http://www.uni-oldenburg.de/en/landeco/research/projects/LEDA/Data%20Files/"
                if(tryCatch(url.exists(url_leda), error =function(e){return(FALSE)},warning=function(w){return(FALSE)})==FALSE){
                    stop("\n\n LEDA website is probably down.\n Please re-run tr8() without selecting LEDA as a source of data \n (or re-try later).\n\n")
                }
                
                local_storage(db="LEDA",directory)
                load(local_leda)
            }else{rearranged<-NULL}
        }
        leda_traits<-leda(species_list,TRAITS=traits_list$LEDA,rearranged=rearranged)
        
        ## retrieve data from BiolFlor
        biolflor_traits<-biolflor(species_list,TRAITS=traits_list$BiolFlor,rest=rest)
        
        ## retrieve data from Pignatti
        pignatti_traits<-pignatti_f(species_list,TRAITS=traits_list$Pignatti)

        ## retrieve flowering periods for Italy
        it_flowering<-get_italian_flowering(species_list,TRAITS=traits_list$Pignatti,rest=rest)
        
        ## add AMF
        ## if AMF is not already downloaded, local_storage is run (but only
        ## if this trait is requred
        TRAIT_AK="Akhmetzhanova"
        ##is the user interested in downloadin Akhmetzhanova?
        if("Myco_infection"%in%traits_list$AMF){
            ## then check if the dataset has already been downloaded
            local_amf<-file.path(directory,"myco.Rda")
            if(file.exists(local_amf)){
                load(local_amf)}else{
                ##  ## otherwise download it now
                ## if(tryCatch(nsl("www.cran.r-project.org"), error =function(e){return(FALSE)},warning=function(w){return(FALSE)})==FALSE){
                ##      stop("You neither have a working internet connection nor locally stored files from Akhmetzhanova et al.\n  Please re-run tr8() function when your internet connection is working.")
                ##  }

                local_storage(db="Akhmetzhanova",directory)
                load(local_amf)
            }
        }else{myco<-NULL
            TRAIT_AK<-NULL
        }
        
        amf_traits<-retrieve_amf(species_list,TRAITS=TRAIT_AK,rest=rest,data_myco=myco)

        ## check&download mycoflor

        TRAIT_MYC="MycoFlor"
        ##is the user interested in downloadin MycoFlor?
        if("MycoFlor"%in%traits_list$AMF){
            ## then check if the dataset has already been downloaded
            local_amf<-file.path(directory,"MycoFlor.Rda")
            if(file.exists(local_amf)){
                load(local_amf)}else{
                ## ## otherwise download it now
                ## if(tryCatch(nsl("www.cran.r-project.org"), error =function(e){return(FALSE)},warning=function(w){return(FALSE)})==FALSE){
                ##     stop("You neither have a working internet connection nor locally stored files from MycoFlor.\n  Please re-run tr8() function when your internet connection is working.")
                ## }

                local_storage(db="MycoFlor",directory)
                load(local_amf)
            }
        }else{
            MycoFlor<-NULL
            TRAIT_MYC<-NULL
        }

        amf_MycoFlor<-retrieve_MycoFlor(species_list,TRAITS=TRAIT_MYC,rest=rest,data_myco=MycoFlor)



        
        ## check if an already downloaded version of the Catminat database
        ## exists and, if so, use it otherwise download a copy, but only
        ## if at least one Catminat trait is needed
        local_Catminat<-file.path(directory,"catminat.Rda")
        if(file.exists(local_Catminat)){
            load(local_Catminat)}else{
            if(length(traits_list$Catminat)>0){
                local_storage(db="Catminat",directory)
                load(local_Catminat)
            }else{catminat_df<-NULL}
        }
        ##        leda_traits<-leda(species_list,TRAITS=traits_list$LEDA,rearranged=rearranged)
        catminat_traits<-catminat(species_list,TRAITS=traits_list$Catminat,catminat_df)
        
        ## check if an already downloaded version of the BROT database
        ## exists and, if so, use it otherwise download a copy, but only
        ## if at least one BROT trait is needed
        local_BROT<-file.path(directory,"BROT.Rda")
        if(file.exists(local_BROT)){
            load(local_BROT)}else{
            if(length(traits_list$BROT)>0){
                local_storage(db="BROT",directory)
                load(local_BROT)
            }else{BROT_df<-NULL}
        }
        
        brot_traits <- brot_data(species_list,TRAITS=traits_list$BROT)
        
        ## download traits from Electronic Flora of Californa
        efloracal_traits<-eflora(species_list,TRAITS=traits_list$EFlora_Cal)
        
        ## check if an already downloaded version of the PLANTS database
        ## exists and, if so, use it otherwise download a copy, but only
        ## if at least one BROT trait is needed
        local_PLANTS<-file.path(directory,"PLANTS.Rda")
        if(file.exists(local_PLANTS)){
            load(local_PLANTS)}else{
            if(length(traits_list$PLANTS)>0){
                local_storage(db="PLANTS",directory)
                load(local_PLANTS)
            }else{PLANTS_df<-NULL}
        }
        
        PLANT_traits <- PLANTS(species_list,TRAITS=traits_list$PLANTS)
        
        
        ## merge the results
        tr8_traits<-data.frame(species_list,row.names=species_list)
        bibliography=list()
        potential_issues<-c()
        for(i in c(eco_traits,biolflor_traits,leda_traits,pignatti_traits,it_flowering,amf_traits,amf_MycoFlor,catminat_traits,brot_traits,PLANT_traits,efloracal_traits)){
            ## merge the dataframes only if they contain data
            if(!is.null(i@results))
            {
                ## clean dataframe column names
                i@results<-column_conversion(i@results)
                ## update the bibliography (Adding the required sources
                
                bibliography[[i@bibliography]]=names(i@results)
                tr8_traits=merge(tr8_traits,i@results,by.x=0,by.y=0,all=TRUE)
                row.names(tr8_traits)<-tr8_traits$Row.names
                tr8_traits<-tr8_traits[,-1,drop=FALSE]
                potential_issues<-c(potential_issues,i@issues)
            }
        }

        ## remove column species_list
        row_names<-row.names(tr8_traits)
        ## names_columns<-names(tr8_traits)[!(names(tr8_traits)%in%c("Row.names","species_list"))]
        ## names(tr8_traits)<-names_columns
        ## tr8_traits<-as.data.frame(tr8_traits[,names_columns],row.names = row_names)
        tr8_traits<-tr8_traits[,!(names(tr8_traits)%in%c("Row.names","species_list")),drop=FALSE]

        obj<-new("Tr8")
        ##    obj@double_names<-unique(c(eco_traits@double_names,leda_traits@double_names))
        ##    obj@not_valid<-intersect(intersect(eco_traits@not_valid,leda_traits@not_valid),pignatti_traits@not_valid)

        ## biolflor_clean is not needed any more
        ## tr8_traits<-biolflor_clean(tr8_traits)
        tr8_traits<-column_conversion(tr8_traits)

        
        if(synonyms==TRUE){
            
            reference_names<-ldply(lapply(reference_names,ldply))
            names(reference_names)<-c("original_names","synonyms")
            tr8_traits<-merge(reference_names,tr8_traits,by.x="synonyms",by.y=0)
            ## in this case, where synonyms are required, then
            ## row names is left with "numbers" since many strange coincidences may
            ## happen (eg. two different species may have been found under the
            ## same synonym, eg. using them as row.names would rais an error (and orginal
            ## names cannot be used for the very same reason)
            ##row.names(tr8_traits)<-tr8_traits$synonyms
            ##tr8_traits<-tr8_traits[,names(tr8_traits)!="synonyms"]
        }


        obj@reference<-temp_dframe
        obj@results<-tr8_traits
        obj@bibliography<-bibliography
        obj@issues<-potential_issues
        
                                        #    issues(obj)
        ##return(obj)
                                        #    return(tr8_traits)

        options(op)
        remove(list=c("column_list"), envir = env)    
        return(obj)
    }
}


