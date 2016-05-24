##' This function will make a simple GUI appear: this window containts
##' a tab  for each trait database: the user che choos which traits should be
##' downloaded from the \code{tr8} function.
##'
##' 
##' @title \code{tr8_config} a GUI to configure TR8 package.
##' @return  the function will return a list containind the selected traits for each source of information
##' @seealso tr8()
##' @author Gionata Bocci <boccigionata@@gmail.com>
##' @export tr8_config

tr8_config<-function(){
    ## Here "gWidgets::gbasicdialog" is used since it makes
    ## the window modal (ie. the script has to wait for the
    ## user's choiche before moving on)

    ## get the list of traits 
    ## data(column_list)
    ## convert it to a data frame
    env<-new.env()
    data("column_list",envir=env)
    column_list<-get("column_list",envir = env)
    
    temp_dframe<-ldply(column_list)
    names(temp_dframe)<-c("long_code","short_code","description","db")

    ## get traits for ecoflora
    temp_Ecoflora<-temp_dframe[temp_dframe$db=="Ecoflora",c("long_code","description")]
    ## get traits for pignatti
    temp_Pignatti<-temp_dframe[temp_dframe$db=="Pignatti",c("long_code","description")]
    ## get traits for biolflor
    temp_BiolFlor<-temp_dframe[temp_dframe$db=="BiolFlor",c("long_code","description")]
    ## get traits for leda
    temp_LEDA<-temp_dframe[temp_dframe$db=="LEDA",c("long_code","description")]
    ## get traits for AMF
    temp_AMF<-temp_dframe[temp_dframe$db=="AMF",c("long_code","description")]
    ## get traits from catminat
    temp_Catminat<-temp_dframe[temp_dframe$db=="Catminat",c("long_code","description")]
    ## get traits from BROT
    temp_BROT<-temp_dframe[temp_dframe$db=="BROT",c("long_code","description")]
    ## get traits from Electronic flora of California
    temp_efloracal<-temp_dframe[temp_dframe$db=="EFlora_Cal",c("long_code","description")]
    ## get traits from PLANTS:
    ## traits from PLANTS are split into 2 panels because
    ## they do not all fit in a single one
    temp_PLANTS<-temp_dframe[temp_dframe$db=="PLANTS",c("long_code","description")]
    split_df<-ceiling(nrow(temp_PLANTS)/2)
    temp_PLANTS_A<-temp_PLANTS[1:split_df,]
    temp_PLANTS_B<-temp_PLANTS[(split_df+1):nrow(temp_PLANTS),]
    
    
    
    ## create the main MODAL window
    window<-gbasicdialog(title="Traits selector for TR8")
    glabel("\nSelect the traits you want to download \nfrom the various databases\n",container=window)
    ## create a group which will containt sheets
    g<-ggroup(container=window)
    ## create a 'notebook' object which contains different sheets
    nb=gnotebook(container=window)
    ## create a sheet of traits for 'biolflor' 
    res_BiolFlor<-gcheckboxgroup(temp_BiolFlor$description,container=nb,label="BiolFlor")
    ## create a sheet of traits for 'leda' 
    res_LEDA<-gcheckboxgroup(temp_LEDA$description,container=nb,label="LEDA")
    res_Ecoflora<-gcheckboxgroup(temp_Ecoflora$description,container=nb,label="Ecoflora")
    res_Pignatti<-gcheckboxgroup(temp_Pignatti$description,container=nb,label="Pignatti")
    res_AMF<-gcheckboxgroup(temp_AMF$description,container=nb,label="AMF")
    res_Catminat<-gcheckboxgroup(temp_Catminat$description,container=nb,label="Catminat")
    res_BROT<-gcheckboxgroup(temp_BROT$description,container=nb,label="BROT")
    res_efloracal<-gcheckboxgroup(temp_efloracal$description,container=nb,label="Eflora Calif.")
    ##res_PLANTS<-gcheckboxgroup(temp_PLANTS$description,container=nb,label="PLANTS")
    ## two panels for PLANTS traits
    res_PLANTS_A<-gcheckboxgroup(temp_PLANTS_A$description,container=nb,label="PLANTS (1)")
    res_PLANTS_B<-gcheckboxgroup(temp_PLANTS_B$description,container=nb,label="PLANTS (2)")
    
    visible(window,TRUE)

    
    
    res_BiolFlor<-svalue(res_BiolFlor)
    res_LEDA<-svalue(res_LEDA)
    res_Ecoflora<-svalue(res_Ecoflora)
    res_Pignatti<-svalue(res_Pignatti)
    res_AMF<-svalue(res_AMF)
    res_Catminat<-svalue(res_Catminat)
    res_BROT <- svalue(res_BROT)
    res_PLANTS <- c(svalue(res_PLANTS_A),svalue(res_PLANTS_B))
    res_efloracal<-c(svalue(res_efloracal))
    
    res_BiolFlor<-fix_values(res_BiolFlor,temp_BiolFlor)
    res_LEDA<-fix_values(res_LEDA,temp_LEDA)
    res_Ecoflora<-fix_values(res_Ecoflora,temp_Ecoflora)
    res_Pignatti<-fix_values(res_Pignatti,temp_Pignatti)
    res_AMF<-fix_values(res_AMF,temp_AMF)
    res_Catminat<-fix_values(res_Catminat,temp_Catminat)
    res_BROT <- fix_values(res_BROT,temp_BROT)
    ## get the chosen traits from the 2 PLANTS panel and merge them in a single vector
    ##res_PLANTS<-c(res_PLANTS,temp_PLANTS)
    res_PLANTS<- fix_values(res_PLANTS,temp_PLANTS)
    res_efloracal<-fix_values(res_efloracal,temp_efloracal)
    
    traits_list<-list("BiolFlor"=res_BiolFlor,"LEDA"=res_LEDA,"Ecoflora"=res_Ecoflora,"Pignatti"=res_Pignatti,"AMF"=res_AMF,"Catminat"=res_Catminat,"BROT"=res_BROT,"PLANTS"=res_PLANTS,"efloracal"=res_efloracal)
    return(traits_list)

}

fix_values<-function(TEMP_VAR,DF){
    
    if(length(TEMP_VAR)==0){
        TEMP_VAR<-c()}else{
    TEMP_VAR<-with(DF,long_code[description%in%TEMP_VAR])
     }
    
    return(TEMP_VAR)
}
