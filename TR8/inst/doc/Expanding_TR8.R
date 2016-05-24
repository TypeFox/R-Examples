### R code from vignette source 'Expanding_TR8.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: aa (eval = FALSE)
###################################################
## res<-new("results")


###################################################
### code chunk number 2: zz (eval = FALSE)
###################################################
## res@results<-NULL


###################################################
### code chunk number 3: bb (eval = FALSE)
###################################################
## res@results<-datatraits


###################################################
### code chunk number 4: cc (eval = FALSE)
###################################################
## res@bibliography<-"Di Sarli, C. and Troilo, A., 2014. 
##    TRAITS: A new web traitbase for the flora of Argentina. 
##    http://www.pichuco.edu"


###################################################
### code chunk number 5: dd (eval = FALSE)
###################################################
##  return(res)


###################################################
### code chunk number 6: tt (eval = FALSE)
###################################################
## it_flowering<-get_italian_flowering(species_list,
##                           TRAITS=traits_list$Pignatti,rest=rest)
## my_exp<-funct_tr(species_list,TRAITS=traits_list$New_db,rest=rest)


###################################################
### code chunk number 7: zx (eval = FALSE)
###################################################
##         for(i in c(eco_traits,biolflor_traits,
##                    leda_traits,pignatti_traits,it_flowering,amf_traits)){
##             ## merge the dataframes only if they contain data
##        ...
##    }


###################################################
### code chunk number 8: zy (eval = FALSE)
###################################################
##         for(i in c(my_exp,eco_traits,biolflor_traits,
##                    leda_traits,pignatti_traits,it_flowering,amf_traits)){
##             ## merge the dataframes only if they contain data
##        ...
##    }


###################################################
### code chunk number 9: a (eval = FALSE)
###################################################
## column_list<-list(
##     ## already existing traits
##     ## ...
##     ## ...
##     "height"=c("height","height of a species","New_db"),
##     "dispersal_type"=c("disp_type","Typology of dispersal","New_db"),
##     "clonality"=c("clonality","Type of clonal species","New_db")
##     )


###################################################
### code chunk number 10: sa (eval = FALSE)
###################################################
## library(plyr)
## tp<-ldply(column_list)[2:4]
## names(tp)<-c("short_code","description","db")
## save(tp,file="TR8/man/available_tr8.Rd")


