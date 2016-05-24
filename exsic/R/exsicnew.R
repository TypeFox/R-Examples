#' Provides botanists with convenience functions to create exsiccatae indices
#' 
#' The tool allows creating simple specimen indices as found in
#' taxonomic treatments based on a table of specimen records. 
#' An example file of tabulated speciment data is provided. In addition,
#' four different exsiccatae styles are provided.
#' The naming of the columns in the specimen table follows largely the conventions used in the BRAHMS 
#' software package.
#' Each specimen record must at least have content in the following nine fields: 
#' id, genus, species, collcite, number, colldate, country, majorarea, minorarea.
#' If not present, the fields are added and filled with dummy values like 's.d.' for no date or 'Unknown
#' country/area'.
#' Highly recommended fields include: collector, addcoll.
#' Optional fields include: locnotes, phenology, elevation, latitude, longitude, and dups
#' The produced indices will sort countries and species alphabetically. Within a country 
#' records will be sorted alphabetically by 'majorarea' (if present) and by collector and 
#' collecting nunber. 
#' A web page in standard html format is created based on a template.
#' The template may be changed and specified in most word processing software.
#' The package provides one main function 'exsic'. 
#' See the example in this section on how to access it. 
#' @name exsic-package
#' @aliases exsic-package
#' @author Reinhard Simon, David M. Spooner
#' @example inst/examples/exsic.R
#' @docType package
NULL

#library(markdown)
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".obj1", "obj2"))

globalVariables(c("format.SBMG", "format.PK", "format.ASPT", "format.NYBG"))

 
#' @name potato
#' @title Wild potato specimen sample list 
#' @format Contains 1000 observations of 16 variables. The variable names conform largely to the
#'  BRAHMS standard.
#' \itemize{
#'  \item{"genus"} {Genus name}
#'  \item{"species"} {Species name}
#' }
#' @note Access the data using system.file("samples/exsic.csv", package = "exsic")
#' @docType data
#' @author David M. Spooner with format edits from R. Simon
#' @family datasets
#' @aliases potato
#' @export
#' @keywords dataset
NULL


#' @name format.SBMG
#' @title A formatting scheme following SBMG (Systematic Botany Monographs) conventions for an exsiccatae record.
#'
#' @format A dataframe with 9 entries and three variables
#' \itemize{
#'  \item{"field"} {A column name present in the exsiccatae table.}
#'  \item{"style"} {The text style applied to that element. One of: none, bold, italics, capitals, 
#'  underline or underscore, uppercase, and () or parentheses}
#'  \item{"sept"}  {A seperator following the element. One of '', ' ', ', ', '. ', ': ' and ';' }
#'  \item{"comments"} {}  
#' }
#' 
#' @docType data
#' @author R. Simon
#' @family formats
#' @aliases format.SBMG
#' @keywords dataset
#' @export
NULL

#' @name format.ASPT
#' @title A formatting scheme following ASPT (American Society of Plant Taxonomists) conventions for an exsiccatae record.
#'
#' @format A dataframe with 9 entries and three variables
#' \itemize{
#'  \item{"field"} {A column name present in the exsiccatae table.}
#'  \item{"style"} {The text style applied to that element. One of: none, bold, italics, capitals, 
#'  underline or underscore, uppercase, and () or parentheses}
#'  \item{"sept"}  {A seperator following the element. One of '', ' ', ', ', '. ', ': ' and ';' }
#'  \item{"comments"} {}  
#' }
#' 
#' @docType data
#' @author R. Simon
#' @family formats
#' @source \url{http://www.aspt.net/publications/sysbot/checklist_systbot.php}
#' @aliases format.ASPT
#' @keywords dataset
NULL


#' @name format.NYBG
#' @title A formatting scheme following NYBG (New York Botanical Garden) conventions for an exsiccatae record.
#' @description In comparison to the two other included specimen citation formats this includes an additional field
#' to handle phenological information at the time of observation.
#'
#' @format A dataframe with 10 entries and three variables
#' \itemize{
#'  \item{"field"} {A column name present in the exsiccatae table.}
#'  \item{"style"} {The text style applied to that element. One of: none, bold, italics, capitals, 
#'  underline or underscore, uppercase, and () or parentheses}
#'  \item{"sept"}  {A seperator following the element. One of '', ' ', ', ', '. ', ': ' and ';' }
#'  \item{"comments"} {}  
#' }
#' 
#' @docType data
#' @keywords datasets
#' @family formats
#' @author R. Simon
#' @source \url{http://www.nybg.org/botany/ofn/fn-gd3.htm}
#' @aliases format.NYBG
#' @keywords dataset
NULL


#' @name format.PK
#' @title A formatting scheme following PK (Phytokeys) conventions for an exsiccatae record.
#'
#' @format A dataframe with 9 entries and three variables
#' \itemize{
#'  \item{"field"} {A column name present in the exsiccatae table.}
#'  \item{"style"} {The text style applied to that element. One of: none, bold, italics, capitals, 
#'  underline or underscore, uppercase, and () or parentheses}
#'  \item{"sept"}  {A seperator following the element. One of '', ' ', ', ', '. ', ': ' and ';' }
#'  \item{"comments"} {}  
#' }
#' 
#' @docType data
#' @author R. Simon
#' @family formats
#' @source \url{http://www.pensoft.net/journals/phytokeys/about/Author Guidelines}
#' @aliases format.PK
#' @keywords dataset
NULL



#' @name sort.specs
#' @title A sorting and filtering configuration data.frame
#' 
#'  
#' @docType data
#' @author R. Simon
#' @aliases sort.specs
NULL

sort.specs = as.data.frame(cbind(
  species = "chacoense;acaule;brevicaule",
  country = "Peru;Bolivia;Argentina"
  ), stringsAsFactors = FALSE)


format.element <- function (elt, styles) {
  elt = str_trim(elt)
  if( !is.na(elt) && elt != "NA" && elt!="" && elt != " "){
  sty = str_split(styles,";")[[1]]
  sty = unique(sty)
  for(ii in 1:length(sty)){
    
      if(str_length(elt) > 0){
        if(sty[ii] == "uppercase") {
          elt = toupper(elt)
        }
        if(sty[ii] == "italics") {
          elt = paste("*", elt, "*", sep="")
        }
        if(sty[ii] == "bold") {
          elt = paste("**", elt, "**", sep="")
        }
        if(sty[ii] == "[]") {
          elt = paste("[", elt, "]", sep="")
        }
        
        if(sty[ii] == "()" | sty[ii] == "parentheses") {
          elt = paste("&#40;", elt, "&#41;", sep="")
        }
        if(sty[ii] == "underline" | sty[ii] == "underscore") {
          elt = paste("<span style=\"text-decoration:underline\">", elt, "</span>", sep="")
        }
        if(sty[ii] == "capitals") {
          elt = paste("<span style=\"font-variant: small-caps\">", elt, "</span>", sep="")
        }
        #print(el)
      } else {
        elt = ""
      } 
  }
  } else {
    elt = ""
  }
 # write(elt,"out.txt",app=T)
  return(elt)
}



format.exsic.rec = function( record = "", form.specs = format.SBMG) {
  if(!is.data.frame(record)) return("")
  n = ncol(record)
  for(i in 1:n) record[i] = as.character(record[i])
  fsp = form.specs
  fsp = fsp[!fsp$field %in% c("species","country","majorarea","minorarea"),]
  res = ""
    for(i in 1:nrow(fsp)){
      x = ""
      #try({
        nms = names(record)
        if(fsp$field[i] %in% nms){
        elt = record[fsp$field[[i]]]
        elt = str_trim(as.character(elt))
        x = ""
        styles = fsp$style[i]
       
        elt = format.element(elt, styles)
        
        if(elt!=""){
          y = str_trim( paste(x, elt, sep="") )
          x = paste(y, fsp$sept[i],  sep="")   
        } 
        
        }
      lst = str_sub(res,str_length(res), str_length(res))
      if(lst == " " & (x==", " | x == ";")) {
        res = str_trim(res)
      }
      lst = str_sub(res,str_length(res), str_length(res))

      if(lst == "," & x==", ") {
        res = paste(res, " ", sep="")
      } else   res = paste(res, x, sep="")
      
      
    }
  res = str_trim(res)
  sf = str_sub(res, str_length(res),str_length(res))
  #print(sf)
  if(sf==",") {
    res = paste(str_sub(res, 1,(str_length(res)-1)),fsp$sept[length(fsp$sept)],sep="")
  }
  if(sf!="," & sf!=";") {
   
    res = paste(res,";",sep="")
  }
  return(res)
}


compl.table <- function(recs){
  must.have = c("genus", "species", "country", "majorarea", "minorarea","collcite","number")
  missing.r = c("Unknown genus",      "sp.","Unknown country","Unknown major area","Unknown minor area","Anonymous","s.n.")
  n = nrow(recs)
  nr= names(recs)
  nr = nr[!is.na(nr)]
  must.add  = must.have[!must.have %in% nr]
  must.use  = missing.r[!must.have %in% nr]
  if(length(must.add) > 0){
  for(i in 1:length(must.add)){
    acol = rep(must.use[i],n)
    recs = cbind(recs, acol)
    names(recs)[ncol(recs)] = must.add[i]
  }
  for(i in 1:ncol(recs)) recs[,i] = as.character(recs[,i])
  }
  
  atbl = recs[,must.have]
  for(i in 1:ncol(atbl)){
    atbl[atbl[,i]=="",i] = missing.r[i]
  }
  
  if(any( atbl$genus == "Unknown" )){
    atbl[atbl$genus == "Unknown","species"] = "sp."
  }
  recs[,must.have] = atbl
  recs
}

sort.tbl <- function(atbl, sortfilter=NULL){
  ss = sortfilter
  s.ct = "all"
  s.sp = "all"
  if(!is.null(ss)){
    if("genus" %in% names(ss)){
      if(str_detect(ss$genus,";")) {
        s.gr = str_split(ss$genus,";")[[1]]    
      } else {
        s.gr = ss$genus
      }
    } else {s.gr = NA}
    if("species" %in% names(ss)){
      if(str_detect(ss$species, ";")){
        s.sp = str_split(ss$species,";")[[1]]  
      } else {
        s.sp = ss$species
      }
      
    } else {s.sp= NA}
    if("country" %in% names(ss)){
      if(str_detect(ss$country, ";")){
        s.ct= str_split(ss$country,";")[[1]]  
      } else {
        s.ct = ss$country
      }
    } else {s.ct = NA}
  } # is null filter
  
  
  tbl.nm = names(atbl)
  if("number" %in% tbl.nm) {
    #atbl=  arrange(atbl, number)
    atbl = atbl[order(atbl[,"number"]) ,]
  }
  if("collcite" %in% tbl.nm) {
    atbl = atbl[order(atbl[,"collcite"]) ,]
  } 
  
  adf = atbl
  
  if("minorarea" %in% tbl.nm) {
    adf = adf[-c(1:nrow(adf)),]
    s.mi = sort(unique(atbl$minorarea))
    s.mi = s.mi[s.mi!='Unknown']
    for(i in 1:length(s.mi)){
      adf = rbind(adf, atbl[atbl[,"minorarea"] == s.mi[i],])
    }
    adf = rbind(adf, atbl[atbl[,"minorarea"] == "Unknown",])
    atbl= adf
    
  }
  
  if("majorarea" %in% tbl.nm){
    adf = adf[-c(1:nrow(adf)),]
    s.ma = sort(unique(atbl$majorarea))
    s.ma = s.ma[s.ma!='Unknown']
    for(i in 1:length(s.ma)){
      adf = rbind(adf, atbl[atbl[,"majorarea"] == s.ma[i],])
    }
    adf = rbind(adf, atbl[atbl[,"majorarea"] == "Unknown",])
    atbl= adf
  }
  
  if(s.sp[1] == "all"){
    atbl = atbl[order(atbl[,"species"]) ,]
  } else {
    adf = atbl
    adf = adf[-c(1:nrow(adf)),]
    for(i in 1:length(s.sp)){
      adf = rbind(adf, atbl[atbl[,"species"] == s.sp[i],])
    }
    atbl = adf
  }
  
  if(s.ct[1] == "all"){
    atbl = atbl[order(atbl[,"country"]) ,]
  } else {
    adf = atbl
    adf = adf[-c(1:nrow(adf)),]
    #print(s.ct)
    for(i in 1:length(s.ct)){
      #print(s.ct[i])
      adf = rbind(adf, atbl[atbl[,"country"] == s.ct[i],])
    }
    atbl = adf
  }
  
  
  atbl
}

#' Prepare exsic table
#' 
#' Complements missing minimal columns and missing content. If applicable (not null)
#' applies a list of countries and species in the given order to sort the table. Countries and 
#' species not present in a given sortfilter data.frame will be filtered out from the underlying
#' table.
#' 
#' @param atbl a data.frame table
#' @param sortfilter a data frame with two columns (one for species and one for genus; maybe NULL)
#' @return a table in the format expected by the exsic function
#' @aliases prepare.table
#' @author Reinhard Simon
#' @export
#' @family main
prepare.table <- function(atbl, sortfilter = NULL){
  rec = compl.table(atbl)    # fill empty space with meaningfull words
  rec = sort.tbl(rec, sortfilter)
  if(nrow(rec) > 0) class(rec) = c("exsic","data.frame")
  rec = rec[row.names(rec)!="NA", ]
  rec
}

#' Creates an index of specimen citations
#' 
#' 
#' @aliases index.citations
#' @param atable a table with exsiccatae records
#' @param format a data frame containing formatting options
#' @return a string in markdown format
#' @author Reinhard Simon
#' @family main
#' @export
index.citations = function(atable, format= format.SBMG #, sortfilter = sort.specs
                            ) {
  stopifnot(is.format(format))
  group.ma = TRUE
  group.sty= " - "
  try({
    if(format[format$field == "group.majorarea","style"] == "no") group.ma = FALSE
    group.sty = format[format$field == "group.majorarea","sept"]
  })
  
  
  
  recs = atable
  fsp = format
  rec = recs
  if(nrow(rec) < 1){
    res = "**The filter applied did not give leave any records in the table. Please check the sortfilter.**"
  } else {
  
  pb <- txtProgressBar(0, nrow(rec), style=3)
  #update.pb(pb, 1, est)
  
  gsp = paste(rec$genus, rec$species)
  rec = cbind(rec, gsp)
  gsp = unique(paste(rec$genus, rec$species))
  species.lev = gsp
  country.lev = unique(rec$country)

  xy=0
  res = ""
    for(i in 1:length(gsp)){
      ele = paste(i,". ", gsp[i], sep="")
      sty = fsp[fsp=="species","style"]
      ele = format.element(ele, sty)
      res = paste(res," \n\n\n ",ele," ", sep="")
      
#      res = paste(res," \n\n\n **_",i, ". ",gsp[i],"_** ", sep="")
      s.rec = rec[rec$gsp == species.lev[i],]
      
      sp.country = unique(s.rec$country)

        scountry.lev =country.lev[ country.lev %in% sp.country]
        for(j in 1:length(scountry.lev)){
          
          #ele = paste(scountry.lev[j], fsp[fsp=="country","sept"], sep="")
          ele = scountry.lev[j]
          sty = fsp[fsp=="country","style"]
          sep = fsp[fsp=="country","sept"]
          ele = format.element(ele, sty)
          if(ele != ""){
            res = paste(res," \n ",ele,sep, sep="")  
          }
          
          
          #res = paste(res," \n **",scountry.lev[j],".** ", sep="")
          sc.rec = s.rec[s.rec$country == scountry.lev[j],]
            ma.lev = unique(sc.rec[,"majorarea"])  
            for(k in 1:length(ma.lev)){
              if(length(ma.lev) == 0 ){
              } else if (length(ma.lev)==1 ) {
                txt = ma.lev[k]
                if(ma.lev[k]=="Unknown" | is.na(ma.lev[k])) txt = "Unknown major area"
                 
              } else {
                txt = ma.lev[k]
                #res = paste(res," <span style=\"font-variant: small-caps\">",ma.lev[k],".</span> ", sep="")  
              }
              #txt = paste(txt,fsp[fsp=="majorarea","sept"] ,sep="")
              sep = fsp[fsp=="majorarea","sept"] 
              sty = fsp[fsp=="majorarea","style"] 
              tx1 = format.element(txt, sty)
              if(txt != ""){
                res = paste(res,tx1, sep, sep="")   
              } else {
                res = paste(res,txt, sep, sep="")     
              }
              
              
              scm.rec = sc.rec[sc.rec$majorarea == ma.lev[k],]
              mia.lev = unique(scm.rec[,"minorarea"])
              pp = length(mia.lev)
              if(pp > 0){
              
              for(p in 1:pp){
                if (!is.na(mia.lev[p]) & p>1 ) {
                  txt = mia.lev[p]
                  if(str_detect(mia.lev[p], "Unknown") ) {
                    txt <- "Unknown minor area"
                  }  
                } else if ((!str_detect(mia.lev[p], "Unknown") | is.na(mia.lev[p])) & p ==1){
                  if(!str_detect(ma.lev[k],"Unknown")) {
                    txt =  mia.lev[p]
                  }
                } else txt = "Unknown minor area" #res = paste(" ",res, "Unknown minor area", ": ", sep="")  
                
                #txt = paste(txt,fsp[fsp=="minorarea","sept"] ,sep="")
                sep = fsp[fsp=="minorarea","sept"] 
                sty = fsp[fsp=="minorarea","style"] 
                tx1 = format.element(txt, sty)
                if(tx1 != ""){
                  res = paste(res,tx1, sep, sep="") 
                } else {
                  res = paste(res,txt, sep="")   
                }
                
                
                scmi.rec = scm.rec[scm.rec$minorarea == mia.lev[p],]
                for(l in 1:nrow(scmi.rec)){
                  xy = xy+1
                  ###################################
                  # Central call to format one record
                  x = format.exsic.rec(scmi.rec[l,], fsp)  
                  update.pb(pb, xy, "")
                  res = paste(res,x," ", sep="")
                  #print(res)
                }
                res = str_trim(res)
                
                fs2 = str_sub(res, str_length(res)-1, str_length(res))
                if(fs2==".;") res = paste(str_sub(res, 1 , str_length(res)-1)," ", sep="")
                fs1 = str_sub(res, str_length(res), str_length(res))
                if(fs1==";") res = paste(str_sub(res, 1 , str_length(res)-1),". ", sep="")
              } # end for p
              
              } else {
                for(l in 1:nrow(scm.rec)){
                  x = format.exsic.rec(scm.rec[l,], fsp)  
                  res = paste(res,x," ", sep="")
                  #print(res)
                 }
              }
              if(k < length(ma.lev) && group.ma ){
                gs = " -- "
                if(group.sty != " - ") gs = group.sty
                res = paste(res,gs,sep="")
              }
            }

        
      }
    }
  }
  res = paste(res,"\n\r", sep="")
  close(pb)
  res
}


#' Creates a simple index of countries.
#' 
#' 
#' @aliases index.countries
#' @param atable a table with exsiccatae records
#' @return a string in markdown format
#' @author Reinhard Simon
#' @family main
#' @export
index.countries <-function(atable=NULL) {
  if(is.null(atable)) {
    res = "\n\n*No table given.*"
  }
  nms = names(atable)
  if("country" %in% nms) {
    cn = unique(atable$country)
    n = length(cn)
    res=paste(1:n,". ",cn,sep="")
    res = paste(res,collapse="\n")
  } else {
    res = "\n\n*Table does not contain a 'country' column*"
  }
  res
}

#' Creates a simple index of species.
#' 
#' 
#' @aliases index.species
#' @param atable a table with exsiccatae records
#' @return a string in markdown format
#' @author Reinhard Simon
#' @family main
#' @export
index.species <-function(atable=NULL) {
  if(is.null(atable)) {
    res = "\n\n*No table given.*"
  }
  nms = names(atable)
  if("species" %in% nms & "genus" %in% nms) {
    cl = "\n"
    gspecies = unique(paste(substr(atable$genus,1,1),". ",atable$species,sep=""))
    sp = paste(". *",gspecies,"*",sep="")
    n = 1:length(sp)
    x=paste(n,sp, sep="")
    res = paste(x,collapse=cl)
  } else {
    res = "\n\n*Table does not contain a 'species' column*"
  }
  res = paste(res, "\n\r", sep="")
  res
}


#' A condensed list of species (listed within a line).
#'
#' Should only be used within the template file.
#'
#' @aliases index.species.short
#' @param atable a table with exsiccatae records
#' @author Reinhard Simon
#' @return a string in markdown format
#' @family main
#' @export
index.species.short<-function(atable=NULL){
  if(is.null(atable)) {
    res = "\n\n*No table given.*"
  }
  nms = names(atable)
  if("species" %in% nms & "genus" %in% nms) {
    cl = ", "
    gspecies = unique(paste(substr(atable$genus,1,1),". ",atable$species,sep=""))
    sp = paste(". *",gspecies,"*",sep="")
    n = 1:length(sp)
    x=paste(n,sp, sep="")
    res = paste(x,collapse=cl)
    res = paste("[",res,"]\n\r", sep="")
  } else {
    res = "\n\n*Table does not contain a 'species' column*"
  }
  res
}


#' Creates the index of collectors and their specimes.
#' 
#' Should only be used within the template file.
#' 
#' @aliases index.collections
#' @param atable a table with exsiccatae records
#' @param format a format.exsic table (check with is.format)
#' @return a string in markdown format
#' @author Reinhard Simon
#' @family main
#' @export
index.collections <-function(atable = NULL, format = format.SBMG){
  res=""
  stopifnot(is.format(format))
  group.specimens.bool = TRUE
  group.join = format[format$field=="group.specimens","sept"]
  species.name = FALSE
  species.sty = "();number;none"
  species.sept=", "
  
  try({
    if(format[format$field=="group.specimens","style"]=="no") {
      group.specimens.bool = FALSE
    } 
    
  }, silent = TRUE)
  
  try({
    species.sty = format[format$field=="species.referral","style"]
    species.sept= format[format$field=="species.referral","sept"]
    if(str_detect(species.sty,"name")) species.name = TRUE
  }, silent = TRUE)
  
  if(nrow(atable)!=0){
    col.ord = sort(unique(atable$collcite))
    #dbs = arrange(mstr,mstr$cite, mstr$number)
    try({
      atable[atable$number=="s.n.","number"]=""
      atable$number = as.numeric(atable$number)
      atable = atable[order(atable[,"number"]) ,]  
      atable$number = as.character(atable$number)
      atable[is.na(atable$number),"number"]="s.n."
    }, silent= T)
    
    
    #atable = atable[order(atable[,"collcite"]) ,]
    atable$species = str_trim(atable$species)
    atable$number = str_trim(atable$number)
    atable$collcite = str_trim(atable$collcite)
    
    atable = atable[with(atable, order(collcite, species, number)) ,]
    
    
    
    spi = unique(atable$species)
    dbs = atable
    
    for(i in 1:length(col.ord)){
      res = paste(res,col.ord[i]," ",sep="")
      dbt = dbs[dbs$collcite==col.ord[i],]
      
      res1 = res
      res =""
      for(k in 1:nrow(dbt)){
        if(species.name){
          elt = format.element(dbt$species[k], species.sty)
          res = paste(res,dbt$number[k]," ",elt,sep="")
          #res = paste(res,dbt$number[k]," (*",dbt$species[k],"*)",sep="")
        } else {
          elt = format.element(which(spi==dbt$species[k]), species.sty)
          res = paste(res,dbt$number[k]," ",elt,sep="")
          
          #res = paste(res,dbt$number[k]," (",which(spi==dbt$species[k]),")",sep="")  
        }
        
        if (k<nrow(dbt)){
          res = paste(res,species.sept,sep="")
        } else {
          res = paste(res,".",sep="")
        }
      }
      #print(res)
      if(group.specimens.bool){
        res = group.specimens(res, group.join)  
      }
      res=paste(res1, res,"\n\n ")
    }
  } else {
    res = "**No records found for this combination of countries and species.**\n"
  }
  res = paste(res, "\n\r", sep="")
  res
}

#' Creates a section header
#' 
#' A simple wrapper to return a section header formatted for inclusion in the final document.
#' 
#' @aliases section.exsic
#' @param title a section title
#' @return a string in markdown format
#' @author Reinhard Simon
#' @family main
#' @export
section.exsic <- function(title="My title") {
  paste(str_trim(title),"\n","-----\n",sep="")
}

#' Writes the exsiccatae text
#' 
#' Writes the text string into a file formatted as HTML.
#' 
#' @aliases write.exsic
#' @param text a text in markdown format
#' @param file a file name
#' @author Reinhard Simon
#' @family main
#' @export
write.exsic <- function(text, file){
  markdownToHTML(text = text, output = file)
}

#' Checks table format
#' 
#' Checks if the table is in the expected format by the principal exsic function.
#' 
#' The format should be created either using reading a table via read.exsic or prepare.table.
#' 
#' @aliases is.exsic
#' @param obj an object
#' @author Reinhard Simon
#' @family helper
#' @export
is.exsic <-function(obj) {
  inherits(obj,"exsic")
}


#' Reads a table in .csv format.
#' 
#' The table should have a certain set of minimum columns following the conventions of 
#' the BRAHMS software.
#' 
#' The function tries to be forgiving and makes a compatible table even if none of the found 
#' columns is compliant with the expected names.
#' 
#' @param file a file path
#' @param sortfilter a dataframe containing instructions on how to sort and filter the table
#' @return a special dataframe conforming to the exsic expectations
#' @author Reinhard Simon
#' @aliases read.exsic
#' @export
#' @family main
read.exsic <- function (file, sortfilter=NULL) {
  exsic.file = file
  exsic.sortfilter = sortfilter
  stopifnot(file.exists(exsic.file))
  data = read.csv(exsic.file, header = TRUE, stringsAsFactors = FALSE)
  stopifnot(nrow(data) > 0 && ncol(data) > 0)
  data = prepare.table(data, sortfilter = exsic.sortfilter)
  data
}

#' Checks if the data frame conforms to the expectations to a formatting configuration
#' data frame.
#' 
#' @param fmt a special data frame
#' @return boolean TRUE if ok
#' @author Reinhard Simon
#' @seealso format.SBMG, format.NYBG
#' @aliases is.format
#' @export
#' @family helper
is.format <- function(fmt) {
  #if(!is.data.frame(fmt)) return(FALSE)
  stopifnot(is.data.frame(fmt))
  nms = names(fmt)
  #if(!all( c("field", "style", "sept") %in% nms)) return(FALSE)
  stopifnot(all( c("field", "style", "sept", "required", "comments") %in% nms))
  
  # check if all obligatory fields are present
  fields = c("species","country","majorarea","minorarea","collcite","number","group.majorarea",
             "species.referral", "group.specimens")
  stopifnot(all(fields %in% fmt$field))
  stopifnot(all(fmt$required %in% c("optional","obligatory")))
  
  x = paste(fmt$style, collapse=";")
  x = unique( str_split(x,";")[[1]] )
  #print(x)
  #if(!all(x %in% c("none","underline","underscore","bold","italics","uppercase","capitals", "()", "number","name")  ) ) return(FALSE)
  stopifnot(all(x %in% c("none","yes","no","underline","underscore","bold","italics","uppercase","capitals", "()", "[]","number","name")  ))
  #x = paste(fmt$sept, collapse="+")
  x = unique( fmt$sept )
  #if(!all( x %in% c(""," ",", ",";",".",":","-"," -- ")  ) ) return(FALSE)
  #print(x)
  #print(x %in% c(""," ",", ",";",".",":","-"," - ") )
  stopifnot(all( x %in% c(""," ",", ",";",". ",":",": ","-"," - ")  ))
  return(TRUE)
}

#' Checks if the data frame conforms to the expectations to sort and filter the exsic database.
#' 
#' @param sf a sortfilter data frame
#' @return boolean TRUE if ok
#' @author Reinhard Simon
#' @aliases is.sortfilter
#' @export
#' @family helper
is.sortfilter <-function(sf){
  if(!is.data.frame(sf)) return(FALSE)
  nms = names(sf)
  if(!all( c("country", "species") %in% nms)) return(FALSE)
  if(nrow(sf) > 1) return(FALSE)
  return(TRUE)  
}

#' Creates three botanical indices (exsiccatae or index of specimens; numerical list of species; 
#' index to numbered collections).
#' 
#' It uses a data.frame expecting a minimum set of columns; if those are not found they will be
#' added and filled with 'unknown' or similar values as will be other columns with missing data
#' but used for sorting the final indices. 
#' 
#' @param data a table in exsic format
#' @param file Path to the .csv file containing the database
#' @param html Path to the resulting .html file
#' @param sortfilter a table containing filters for country and species
#' @param format a table containing format options for elements in the exsiccatae index
#' @param headers A list of text lines for labeling the indices
#' @return boolean TRUE if all steps executed successfully
#' @author Reinhard Simon
#' @aliases exsic
#' @example inst/examples/exsic.R
#' @family main
#' @export
exsic <- function(data, 
                  file=NULL, 
                  html="exsic.html", 
                  sortfilter=NULL, 
                  format=format.SBMG,
                  headers = 
                    c("Citations of Specimens", "Numerical List of Species", 
                      "Index to Numbered Collections",
                      "*Numbers refer to species in the Numerical List.*\n\r")){
  res = ""
  chk = FALSE
  exsic.format = format
  exsic.file = file
  out.file = html
  exsic.sortfilter = sortfilter
  try({
    if (is.null(exsic.format)) {
      exsic.format = format.SBMG
    }
    stopifnot(is.format(exsic.format))
    if(!is.null(exsic.sortfilter)) stopifnot(is.sortfilter(exsic.sortfilter))
    
    if(!is.null(exsic.file)){
      data = read.exsic(exsic.file, exsic.sortfilter)
    } else {
      data = prepare.table(data, exsic.sortfilter)
    }
    stopifnot(nrow(data)>0)
    
    stopifnot(is.exsic(data))
  
    res = section.exsic(headers[1])
    out = index.species.short(data)
    res = paste(res, out, sep="")
    
    out = index.citations(data, format = exsic.format)
    res = paste(res, out, sep="")
    
    out = section.exsic(headers[2])
    res = paste(res, out, sep="")
    
    out = index.species(data)
    res = paste(res, out, sep="")
    
    out = section.exsic(headers[3])
    txt = headers[4]
    res = paste(res, out,txt, sep="")
    
    out = index.collections(data, format = exsic.format)
    res = paste(res, out, sep="")
    
    write.exsic(res, out.file)
    chk = TRUE
  })
  return(chk)
}



