### R code from vignette source 'TR8.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: options
###################################################
options(width = 60)


###################################################
### code chunk number 2: install (eval = FALSE)
###################################################
## install.packages("TR8",dependencies = TRUE)


###################################################
### code chunk number 3: load (eval = FALSE)
###################################################
## library(TR8)


###################################################
### code chunk number 4: devtools (eval = FALSE)
###################################################
## ## install the package
## install.packages("devtools")
## ## load it
## library(devtools)
## ## activate dev_mode
## dev_mode(on=T)
## ## install TR8
## install_github("GioBo/TR8",ref="master")
## ## load it
## library(TR8)
## ## you can now work with TR8 functions
## 
## ## if you want to go back and use the CRAN version
## ## already installed, simply deactivate dev_mode
## dev_mode(on=F)


###################################################
### code chunk number 5: usage (eval = FALSE)
###################################################
## ## a vector containing a list of plant species names
## my_species<-c("Apium graveolens","Holcus mollis","Lathyrus sylvestris")
## ## a vector of traits
## to_be_downloaded<-c("reprod_B","strategy")
## ## now run tr8 and store the results in the my_traits object
## my_traits<-tr8(species_list = my_species,download_list = to_be_downloaded)


###################################################
### code chunk number 6: load_TR8
###################################################
library(TR8)


###################################################
### code chunk number 7: available_traits
###################################################
## see the firs lines of available_tr8 database
head(available_tr8)


###################################################
### code chunk number 8: first_tr8 (eval = FALSE)
###################################################
## my_species<-c("Salix alba","Populus nigra")
## my_traits<-c("h_max","le_area","li_form")
## my_Data<-tr8(species_list = my_species, download_list = my_traits)


###################################################
### code chunk number 9: usage2 (eval = FALSE)
###################################################
## ## see the downloaded data
## print(my_Data)


###################################################
### code chunk number 10: extract_traits (eval = FALSE)
###################################################
## traits_dataframe<-extract_traits(my_Data)


###################################################
### code chunk number 11: lookup (eval = FALSE)
###################################################
## lookup(my_Data)


###################################################
### code chunk number 12: lookup2 (eval = FALSE)
###################################################
## my_lookup<-lookup(my_Data)
## head(my_lookup)


###################################################
### code chunk number 13: issues (eval = FALSE)
###################################################
## my_species<-c("Salix alba","Populus nigra")
## my_traits<-c("h_max","le_area","li_form")
## my_Data<-tr8(species_list = my_species, download_list = my_traits)
## issues(my_Data)


###################################################
### code chunk number 14: dataset (eval = FALSE)
###################################################
## ## suppose veg_data is our dataframe with
## ## plant species as columns and sites as rows
## 
## ## extract species names
## specie_names<-names(veg_data)
## ## use the tr8() function
## ## and tick those traits of interest in the pop-up window
## my_traits<-tr8(species_names,gui_config=TRUE)
## ## print the results
## print(my_traits)


###################################################
### code chunk number 15: bib (eval = FALSE)
###################################################
## bib(my_traits)


###################################################
### code chunk number 16: import (eval = FALSE)
###################################################
## My_data<-read.csv("my_veg_data.csv",
##                   header=T,row.names=1,check.names=F)


###################################################
### code chunk number 17: one (eval = FALSE)
###################################################
## species_names<-names(veg_data)
## checked_names<-tnrs(species_names,source="iPlant_TNRS")
## print(checked_names)


###################################################
### code chunk number 18: tr8_ex1 (eval = FALSE)
###################################################
## my_traits<-tr8(species_names,gui_config = TRUE)
## print(my_traits)


###################################################
### code chunk number 19: issue_workflow (eval = FALSE)
###################################################
## my_traits<-tr8(species_names,gui_config = TRUE)
## issues(my_traits)


###################################################
### code chunk number 20: extract (eval = FALSE)
###################################################
## traits_df<-extract_traits(my_traits)


###################################################
### code chunk number 21: store_to_csv (eval = FALSE)
###################################################
## save(traits_df,file="traits_df.csv")


###################################################
### code chunk number 22: synonims (eval = FALSE)
###################################################
## my_species<-c("Salix alba","Inula viscosa")
## my_traits<-c("h_max","le_area","li_form")
## my_Data<-tr8(species_list = my_species, download_list = my_traits, synonyms=TRUE)


###################################################
### code chunk number 23: vignette (eval = FALSE)
###################################################
## vignette("TR8_workflow")


