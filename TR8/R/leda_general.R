leda_download<-function(url,skip_row,column,out_name){
#    base_url<-"http://www.leda-traitbase.org/LEDAportal/objects/Data_files/"
    ##base_url<-"http://www.uni-oldenburg.de/en/landeco/research/projects/LEDA/Data Files/"
    base_url<-"http://www.uni-oldenburg.de/fileadmin/user_upload/biologie/ag/landeco/download/LEDA/Data_files/"
    url<-paste(base_url,url,sep="")
    downloaded<-read.csv(url,row.names=NULL,skip=skip_row,sep=";",check.names="F")
    
    ## extract only those column of interest
    rearranged<-downloaded[,c("SBS name",column)]
    names(rearranged)<-c("Species","variable")
    
    ## if the variable is numeric, we calculate the mean
    ## if it's a factor, we paste the various levels in
    ## one string
    if(is.numeric(rearranged$variable)){
        rearranged<-with(rearranged,tapply(variable,Species,FUN=function(x){mean(x,na.rm=TRUE)}))
    } else if (is.factor(rearranged$variable)){
        rearranged<-with(rearranged,tapply(variable,Species,FUN = function(x){paste(unique(x),collapse = " + ")}))
    }
    rearranged<-as.data.frame(rearranged)
    ## if varaible should be considered as a factor it
    ## must be converted as such
    if(!is.numeric(rearranged$rearranged)){rearranged$rearranged<-as.factor(rearranged$rearranged)}
    names(rearranged)<-c(out_name)
    return(rearranged)
}





leda_general<-function(url,skip_row,species,column,out_name){
    ## download the original traits
    rearranged<-leda_download(url=url,skip_row=skip_row,column=column,out_name="variable")

    ########
    ######## Mettere qui le opzioni di utilizzo del db scaricato
    #######  Ricorda di chiamare la colonna "variable"
    
    temp_df<-as.data.frame(species,row.names = species)
    ## merge the species list with the retrieved data
    temp_df<-merge(temp_df,rearranged,by.x=0,by.y=0,all.x=TRUE)
    row.names(temp_df)<-as.character(temp_df$Row.names)
    results<-data.frame(temp_df$variable)
    row.names(results)<-row.names(temp_df)
    names(results)<-out_name
    return(results)
}

## age_of_first_flowering<-leda_general(url="age%20of%20first%20flowering.txt",skip_row = 4,column="age of first flowering",out_name="age_first_flowering",species=species)

## branching<-leda_general(url="branching.txt",skip_row = 4,column="branching",out_name="branching",species=species)

## bud_bank_seasonality_soil<-leda_general(url="buds%20seasonality.txt",skip_row = 4,column="budb seas. at soil surface",out_name="bud_seasonality_soil_surface",species=species)

## #bud_bank_seasonality_soil<-leda_general(url="buds%20seasonality.txt",skip_row = 4,column="budb. seas at layer -10-0cm",out_name="bud_seasonality_0-10cm",species=species)


## ## buds_above_10<-leda_general(url="buds%20vertical%20dist.txt",skip_row = 4,column="buds in layer >10 cm",out_name="buds_above_10",species=species)


## ## buds_1_10<-leda_general(url="buds%20vertical%20dist.txt",skip_row = 4,column="buds in layer 1-10cm",out_name="buds_1_10",species=species)


## buoyancy<-leda_general(url="buoyancy.txt",skip_row = 6,column="gen. dispersal type",out_name="dispersal_type",species=species)

                        
## canopy_height<-leda_general(url="canopy%20height.txt",skip_row = 4,column="single value [m]",out_name="mean_height[m]",species=species)


## dispersal<-leda_general(url="dispersal%20type.txt",skip_row = 6,column="dispersal type",out_name="dispersal_type",species=species)


## leaf_distribution<-leda_general(url="leaf%20distribution.txt",skip_row = 7,column="leaf distribution",out_name="leaf_distribution",species=species)


## leaf_dmc<-leda_general(url="LDMC%20und%20Geo.txt",skip_row = 4,column="single value [mg/g]",out_name="leaf_dmc[mg/g]",species=species)

## leaf_mass<-leda_general(url="leaf%20mass.txt",skip_row = 4,column="single value [mg]",out_name="mean_leaf_mass[mg]",species=species)

## leaf_size<-leda_general(url="leaf%20size.txt",skip_row = 4,column="single value [mm^2]",out_name="mean_leaf_size[mm^2]",species=species)

## dispersal_morphology<-leda_general(url="morphology%20dispersal%20unit.txt",skip_row = 6,column="diaspore type",out_name="diaspore_type",species=species)

## growth_form<-leda_general(url="plant%20growth%20form.txt",skip_row = 14,column="plant growth form",out_name="growth_form",species=species)

## life_span<-leda_general(url="plant%20life%20span.txt",skip_row = 14,column="plant lifespan",out_name="life_span",species=species)

## releasing_height<-leda_general(url="releasing%20height.txt",skip_row = 4,column="single value [m]",out_name="releasing_height[m]",species=species)

## sbank<-leda_general(url="seed%20bank.txt",skip_row = 6,column="seed bank type",out_name="seed_bank",species=species)

## ##seed_longevity<-leda_general(url="seed%20longevity.txt",skip_row = 4,column="max longevity",out_name="seed_longevity",species=species)

## seed_mass<-leda_general(url="seed%20mass.txt",skip_row = 4,column="single value [mg]",out_name="mean_seed_mass[mg]",species=species)

## shoot_growth_form<-leda_general(url="shoot%20growth%20form.txt",skip_row = 7,column="shoot growth form",out_name="shoot_growth_form",species=species)

## seed_number_per_shoot<-leda_general(url="SNP.txt",skip_row = 0,column="single value",out_name="seed_number_per_shoot",species=species)

## woodiness<-leda_general(url="ssd.txt",skip_row = 4,column="woodiness",out_name="woodiness",species=species)

## terminal_velocity<-leda_general(url="TV.txt",skip_row = 4,column="single value [m/s^2]",out_name="terminal_velocity[m/s^2]",species=species)


## ##merge dei vaeri dataframe
## vai<-Reduce(function(first,second){
##     merged<-merge(first,second,by.x=0,by.y=0)
##     row.names(merged)<-merged$Row.names
##     merged<-merged[,-1]
## },list(age_of_first_flowering,branching,bud_bank_seasonality_soil,buoyancy,canopy_height,dispersal,leaf_distribution,leaf_dmc,leaf_mass,leaf_size,dispersal_morphology,growth_form,life_span,releasing_height,sbank,seed_mass,shoot_growth_form,seed_number_per_shoot,woodiness,terminal_velocity))
