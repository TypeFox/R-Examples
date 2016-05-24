##' Retrieve data about AMF potential inoculation for a list of
##' plant species passed as argument
##'
##' The returned dataframe comprises one column: \code{Myco_infection}
##'
##' \code{Myco_infection} :  a numeric vector containing percentage of infection as provided
##' by Akhmetzhanova et al.
##' @title retrieve_amf 
##' @param species a vector containing names of plant species 
##' @param TRAITS a vector containing the traits to be downloaded (used as a check for
##' tr8_gui() created variables)
##' @return a data frame
##' @references \itemize{
##' \item Asem A. Akhmetzhanova, Nadejda A. Soudzilovskaia, Vladimir G. Onipchenko,
##' Will K. Cornwell, Vladimir A. Agafonov, Ivan A. Selivanov, and Johannes H. C. Cornelissen.
##' 2012. A rediscovered treasure: mycorrhizal intensity database for 3000 vascular plants
##' species across the former Soviet Union. Ecology 93:689.
##' \samp{http://esapubs.org/Archive/ecol/E093/059/default.htm}
##' }
##' @author Gionata Bocci <boccigionata@@gmail.com>
##'
##' @examples \dontrun{
##' ##My_traits<-retrieve_amf(species_list=c("Abies alba"))
##' }
retrieve_amf<-function(species,TRAITS,rest,data_myco){

    res<-new("results")
    if(is.null(TRAITS)){
        results<-NULL
    }else{

        temp_df<-as.data.frame(species,row.names = species)
        list_temp<-list()
        res_df<-data_myco[data_myco$Species%in%species,c("Species","Intensity.of.mycorrhizal.infection")]
        temp_df<-merge(temp_df,res_df,by.x=0,by.y="Species",all.x=TRUE)
        ## remove NAs
        ##res_df<-res_df[!is.na(res_df$Intensity.of.mycorrhizal.infection),]
        temp_df$species<-as.character(temp_df$species)
        temp_df<-with(temp_df,tapply(Intensity.of.mycorrhizal.infection,species,FUN = function(x){paste(unique(x),collapse = "-")}))
        temp_df<-as.data.frame(temp_df)
        names(temp_df)<-"Myco_infection"
        results<-temp_df
    }
    res@results<-results
    ## remove(list=c("myco"),pos =".GlobalEnv")
    
        
    res@bibliography<-"Akhmetzhanova, A. A., Soudzilovskaia, N. A., Onipchenko, V. G.,\n Cornwell, W. K., Agafonov, V.A., Selivanov, I. A., and Cornelissen, J. H. C., 2012.\nA rediscovered treasure: mycorrhizal intensity database for 3000 vascular plants\nspecies across the former Soviet Union. Ecology 93:689.\n"
    return(res)
}
