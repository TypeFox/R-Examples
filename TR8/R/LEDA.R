leda_extract_from_local_df<-function(local_df,species_list,column_variables){
    as.data.frame(species_list)->spec
    TEMP<-as.data.frame(local_df[,column_variables])
    names(TEMP)<-column_variables
    row.names(TEMP)<-row.names(local_df)
    spec<-merge(spec,TEMP,by.x="species_list",by.y=0,all.x=TRUE)
    DF<-spec[,column_variables,drop=FALSE]
    row.names(DF)<-spec$species_list
    return(DF)
}



## load("data/LEDA_df.rda")

### * leda() function to extract data from local LEDA dataframe


setClass("Leda",representation = list(results="data.frame",not_valid="vector",double_names="vector"))


##' Used to shorten levels of Seed longevity
##' according to LEDA traitbase
##'
##' Short codes are obtained adopting the following rules:
##'
##' \describe{
##' \item{LT}{"long-term persistent"}
##' \item{ST}{"short-term persistent"}
##' \item{T}{"transient"}            
##' \item{P}{"present"}}
##' @title seed_simplify a function to shorten levels for seed longevity 
##' @param tp : a list obtained from the LEDA traitbase
##' @return a vector
##' @author Gionata Bocci <boccigionata@@gmail.com>
seed_simplify<-function(tp){
    tp<-lapply(tp,function(x)(gsub("long-term persistent","LT",x)))
    tp<-lapply(tp,function(x)(gsub("short-term persistent","PT",x)))
    tp<-lapply(tp,function(x)(gsub("transient","T",x)))
    tp<-lapply(tp,function(x)(gsub("present","P",x)))
    tp<-lapply(tp,function(x)(paste(sort(x),collapse="-")))
    tp<-unlist(tp)
    return(tp)
}

#' Extracts functional traits from the LEDA traitbase.
#' 
#' \code{leda} allows the user to extract data from \emph{LEDA_df} which is a
#' \samp{http://www.leda-traitbase.org/LEDAportal/}.
#' subset of the data available on the LEDA traitbase website
#' 
#' The function returns a data.frame with species as rows and LEDA functional
#' traits as columns.  \code{NA} will be used for those traits which do not
#' have values in the LEDA traitbase.  Species names are converted to
#' \emph{accepted} names (\emph{sensu} \code{TNRS}).
#' 
#' @param species_list a vector containing names of plant species 
#' @param TRAITS a vector containing the traits to be downloaded (used as a check for
#' @param rearranged a variable which passes the already downloaded LEDA dataset if this is available (NULL otherwise)
#' #' tr8_gui() created variables)
#' @return dataframe containing traits data and species names as row.names
#' @author Gionata Bocci <boccigionata@@gmail.com>
#' @seealso \code{\link{LEDA_df}}
#' @references Kleyer, M., Bekker, R.M., Knevel, I.C., Bakker, J.P, Thompson, K., Sonnenschein, M.,
#' Poschlod, P., Van Groenendael, J.M., Klimes, L., Klimesova, J., Klotz, S., Rusch, G.M., Hermy, M.,
#' Adriaens, D., Boedeltje, G., Bossuyt, B., Dannemann, A., Endels, P., Götzenberger, L., Hodgson, J.G.,
#' Jackel, A-K., Kühn, I., Kunzmann, D., Ozinga, W.A., Römermann, C., Stadler, M., Schlegelmilch, J.,
#' Steendam, H.J., Tackenberg, O., Wilmann, B., Cornelissen, J.H.C., Eriksson, O., Garnier, E., Peco, B. (2008):
#' The LEDA Traitbase: A database of life-history traits of Northwest European flora.
#' Journal of Ecology 96: 1266-1274.
#' @examples \dontrun{
#' 
#' #My_traits<-leda(species_list=c("Abies alba"))
#' }
leda<-function(species_list,TRAITS,rearranged){

    res<-new("results")
    ## there are some problems with duplicate TNRS names for the LEDA
    ## database, thus I will use the original names
    ##    results<-LEDA_df[LEDA_df$acceptedname%in%species_list,1:7]
    ##    row.names(results)<-results$acceptedname
    
    ## data("LEDA_df")
    ## the final df will be composed of the following traits
    env<-new.env(parent = parent.frame())
    data(column_list,envir = env)
    column_list<-get("column_list",envir=env)

    data(leda_lookup,envir = env)
    leda_lookup<-get("leda_lookup",envir=env)

    leda_traits<-ldply(column_list)
    leda_traits<-with(leda_traits,.id[V3=="LEDA"])

    leda_lu<-ldply(leda_lookup)
    
    if(is.null(TRAITS)){
        res@results<-NULL
    }else{
        if(length(TRAITS)==0){
            tr_list<-leda_traits
        }else{
            tr_list<-TRAITS
        }
        ##TRAITS<-list()

        #column_variables<-leda_lu$V4[leda_lu$.id%in%tr_list]
        column_variables<-leda_lu$V4[leda_lu$.id%in%tr_list]
        
        if(!is.null(rearranged)){
            spec<-leda_extract_from_local_df(rearranged,species_list,column_variables)
        }else{

            leda_subset<-leda_lu[leda_lu$.id%in%tr_list,]
            
            spec<-as.data.frame(species_list)
            row.names(spec)<-species_list

            for(trait in 1:nrow(leda_subset)){
                    extract<-leda_subset[trait,]
                    leda_temp<-leda_general(url=extract$V1 , skip_row=as.numeric(extract$V2), column=extract$V3, out_name=extract$V4,species=species_list)
                    spec<-merge(spec,leda_temp,by.x=0,by.y=0,all.x=TRUE)
                    row.names(spec)<-spec$Row.names
                    names_of_column<-names(spec)[!(names(spec)%in%c("Row.names","species_list"))]
                    spec<-as.data.frame(spec[,!(names(spec)%in%c("Row.names","species_list"))],row.names = species_list)
                    names(spec)<-names_of_column
            }

        }
        names(spec)<-mapvalues(names(spec),leda_lu$V4,leda_lu$.id,warn_missing = FALSE)
        res@results<-spec
    }
    
    ##res@results<-obj@results
    ##obj<-new("Leda",results=spec)
    ##obj@not_valid<-species_list[!species_list%in%LEDA_df$SBS.name]
    stringa<-"Kleyer, M., Bekker, R.M., Knevel, I.C., Bakker, J.P, Thompson, K., Sonnenschein, M., Poschlod, P., \nVan Groenendael, J.M., Klimes, L., Klimesov\303\241, J., Klotz, S., Rusch, G.M., Hermy, M., Adriaens, D.,\nBoedeltje, G., Bossuyt, B., Dannemann, A., Endels, P., G\303\266tzenberger, L., Hodgson, J.G., Jackel, A-K.,\nK\303\274hn, I., Kunzmann, D., Ozinga, W.A., R\303\266mermann, C., Stadler, M., Schlegelmilch, J., Steendam, H.J.,\nTackenberg, O., Wilmann, B., Cornelissen, J.H.C., Eriksson, O., Garnier, E., Peco, B., 2008.\nThe LEDA Traitbase: A database of life-history traits of Northwest European flora.\nJournal of Ecology 96: 1266-1274.\n"
    Encoding(stringa)<-"unicode"
    res@bibliography<-stringa
    remove(list=c("column_list","leda_lookup"), envir = env)    
    return(res)
}


