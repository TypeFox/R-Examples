#' Provides botanists with convenience functions to create exsiccatae indices
#' 
#' The tool allows creating simple specimen indices as found in
#' taxonomic treatments. An example file of tabulated data is provided. 
#' 
#' Each specimen
#' record must at least have content in the following seven fields: 
#' 
#' id, GENUS, SPECIES, collector, number, COLLDATE and ORIGCTY 
#' (names in upper case correspond to fields defined in the MCPD standard). 
#' 
#' Additional usable fields include:
#' 
#' addcoll, majorarea, minorarea, COLLSITE, ELEVATION, DECLATITUDE, DECLONGITUDE and
#' dups.
#' 
#' The tabulated data are checked for presence of those fields and also against a
#' data dictionary. The check against the data dictionary is a crucial part of 
#' checking the consistency and fitness of the data. It can be modified. 
#' 
#' The data dictionary consists of a set of rules. They are
#' defined in standard R statements one per line. A data profile tabulated report is produced 
#' detailing for each rule whether it could be executed, how many total records were found
#' non-compliant and which ones. Any error found will result in only a minimal 
#' error report - no indices will be produced. The rules may
#' point to lookup files to check against: for example a list of countries, genera, 
#' species, collectors, etc.
#' 
#' The produced indices will sort countries and species alphabetically. Within a country 
#' records will be sorted alphabetically by 'majorarea' (if present) and by collector and 
#' collecting nunber. Missing 'numbers' must be indicated by 's.n.'; missing collecting
#' dates following the MCPD standard (e.g.: 20130101 or 201301-- or 2013---- or --------).
#' This will be according transformed into January, 1st 2013; January, 2013 or 's.d.'. 
#' Other optional but missing information will just be omitted.
#'  
#' The geographic coordinate fields 'DECLATITUDE' and
#' DECLONGITUDE' are assumed to be decimal and in WGS84 projection. 
#' 
#' A web page in standard html format is created based on a template
#' in markdown format. The template may be changed and specified. The
#' generated web-page can be edited in most word processing software.
#' All end-results are compiled in sub-directories named after the 
#' 'data file + '__exsic_results_exsic_' + version number + two-digit number'
#' under the current working directory.
#' 
#' The package provides just one main function 'exsic' and a template for creating
#' the indices. See the example in this section on how to access it. The remaining functions
#' are meant for use in the index template.
#' 
#' @name exsic-package
#' @aliases exsic-package
#' @author Reinhard Simon, David M. Spooner
#' @example inst/examples/exsic.R
#' @docType package
NA

library(datadict)
library(mcpd)
#library(R.oo)
library(knitr)
library(markdown)

#options(encoding = 'UTF-8')

char.enc = "UTF-8"


is.windows <- function(){
  str_detect(Sys.getenv("OS"),"Windows")
}


read_file <- function(file){
    res = read.csv(file, stringsAsFactors = FALSE, encoding = char.enc)
    #if(ncol(res)<2) throw("exsic package: file ",file, " has too few columns.")
    return(res)
}

sniff_exsic <- function(x){
  res = TRUE

    found_names = names(x)
    if(length(found_names)<7) {res = FALSE}
  return(res)
}

has_minimumVariables <- function(x){
  names_found = names(x)
  names_expct = c("id","ORIGCTY", "SPECIES", "COLLDATE",
                  "GENUS", "collector", "number")
  res = all(names_expct %in% names_found)
  return(res)
}

has.exsicConfig <- function(exsicFile){
  res = FALSE
  if(!file.exists(exsicFile)) stop() #throw("exsic package: file not found.")
  #bn = basename(exsicFile)
  cfn = str_replace(exsicFile, ".csv", "-config.csv")
  if(file.exists(cfn)) res = TRUE
  return(res)
}

make.exsicConfig <- function(exsicFile){
  cfn = str_replace(exsicFile, ".csv", "-config.csv")
  frm = system.file("samples/config.csv", package = "exsic")
  res = FALSE
  try({
    file.copy(frm, cfn)  
    res = TRUE
  })
  return(res)
}

is.exsicConfig <- function(exsicFile){
  res = FALSE
  try({
    cfn = str_replace(exsicFile, ".csv", "-config.csv")
    df = read.csv(cfn, stringsAsFactors = FALSE, encoding = char.enc)
    res1 = names(df) %in% c("Parameter", "Value")
    
    pre = c("Title","Author(s)", "Version", "Comments", "Owner", "Publisher","Type",
            "Format", "License", "Funding", "Keywords", "Species-order", "Countries")
    pro = df$Parameter
    res2 = pre %in% pro
  })
  res = all(res1, res2)
  return(res)
}


#' is this a file with a table in the expected format?
#' 
#' @aliases is.exsicFile
#' @param file a file name
#' @param config a config file name; if not provided see details
#' @param rules a file containing data variable rules
#' @return boolean
#' @author Reinhard Simon
#' @family exsic
#' @export
is.exsicFile <- function(file = NULL, 
                         config = system.file("samples/config.csv", package = "exsic"), 
                         rules = system.file("samples/rules.R", package = "exsic")) {
  
  # Progressively work through conditions to make sure file exists in required format
  # Protocol any non-compliances for re-use
#   print("is.exsicFile")
#   print(file)
#   print(rules)
  
  res = NA
  msg = c("File name not null: ok", "Valid file name: ok", "File name exists: ok",
          "File reading in text format: ok", "File in tabular format: ok")
  bas.acc = res.notnull = res.ischar = res.fileex = res.canread = FALSE
  res1 = res2 = res3 = res4 = FALSE
  dp = NULL
  
  if(is.null(file)) {
    msg[1] = "File name is null"
  } else {
    res.notnull = TRUE
  }
  
  if(!is.character(file)) {
    msg[2] = paste("exsic package: file ",file, " name is not valid.", sep="")
  } else {
    res.ischar = TRUE
  }
  
  bas.acc = all(res.notnull, res.ischar)
   
  if(bas.acc){  

    if(!file.exists(file)) {
      msg[3] = paste("exsic package: file ",file, " name does not exist.", sep="")
    } else {
      res.fileex = TRUE
    }
    
    try({
      res = read_file(file)
      res.canread = TRUE
      
      if(!res.canread)  msg[4] = "Error reading file."
      
    }, silent=TRUE)
    
  }
  
  bas.acc = all(bas.acc, res.fileex, res.canread)
  
  if(bas.acc){
    res1 = sniff_exsic(res)
    if(!res1) {   
      msg[5] = paste("exsic package:: file: '",file,"'' is not in a tabular comma separated format!",
                     sep="")  
    } 
  }
  
  bas.acc = all(bas.acc, res1)
    
  if(bas.acc){
     res2 = has_minimumVariables(res)
    
    # check against rules
    res3 = FALSE
    if(file.exists(rules)){
      try({
        at = read.csv(file, stringsAsFactors = FALSE)
        ad = read.rules(rules)
        dp = datadict.profile(at, ad)
        res3 = !has.ruleErrors(dp$checks)    
      })
    }
    
    # check if config file exists and if so if it conforms to expected content
    res4 = FALSE
    if(!has.exsicConfig(file)) make.exsicConfig(file)

    try({
      res4 = is.exsicConfig(file)
    })
  }
 
  #print("ok")
  res = all(bas.acc, res2, res3, res4)
  out = list(file= file, ok = res, messages = msg, rule.checks = dp$checks)
  #print(out)
  return(out)

}

makeErrTargetFileName <- function (from, td) {
  file = basename(from)
  to   = file.path(td, file)
}

copy.errTemplate <- function(td, from){
  res = FALSE
  try({
    to   = makeErrTargetFileName(from, td)
    file.copy(from, to)
    res = TRUE
  })
  return(res)
}

make.errReport <- function(td){
  res = FALSE
  try({
    from = system.file("templates/error.Rmd", package = "exsic")
    res1 = copy.errTemplate(td, from)
    to   = makeErrTargetFileName(from, td)
    out  = str_replace(to, ".Rmd", "_out.md")
    knit(to, out)
    html = str_replace(out, ".md", ".html")
    x = markdownToHTML(file=out, output=html)
    res = TRUE
  })
  
  return(res)
}


##############################################

#For internal use in index functions
#Returns a table with Variable, Rule, Path
get.lookup.files <- function(file = "rule.checks.csv"){
  #print("check dir of lookups")
  #print(getwd())
  lu   = read.csv(file, stringsAsFactors = FALSE, encoding = char.enc)
  vars = lu[str_detect(lu$Rule,"is.oneOf"),c("Variable","Rule")]
  fils = str_extract(vars$Rule,"([A-Za-z0-9_-]{1,}).csv")
  hasf = rep(FALSE, nrow(vars))
  # get the full file path either in the current
  # try first in the current wd .. if not in systems directory
  #print(nrow(vars))
  #print(length(fils))
  if(nrow(vars)>0){
  for(i in 1:nrow(vars)){
    try({
      hasf[i] = file.exists(fils[i])
      #print(hasf)
      
      if(!hasf[i]) {
        sf = system.file(file.path("samples",fils[i]), package="exsic")
        hasf[i] = file.exists(sf)
        if(hasf[i]) {
          fils[i] = sf
        } else {
          fils[i] = NA
        }
        
      }# if 
      
    }) # try
  }# for
  }
  Path = fils
  vars = cbind(vars, Path)
  vars$Path = as.character(vars$Path)
  vars = vars[hasf,]
  return(vars)
}

.lookup <- function(aterm, atable){
  atable[atable$VALUE == aterm,"LABELS"]
}

# For internal use
# As part of load.data routine before saving the mstr file
# return processed mstr table
apply.lookup <- function(mstr){
  glf = get.lookup.files()
  if(nrow(glf) == 0) return(mstr)
  
  n = nrow(glf)
  for(i in 1:n){
    fn = glf$Path[i]
    at = read.csv(fn, stringsAsFactors = FALSE, encoding = char.enc)
    mstr[[glf$Variable[i]]] = as.character(sapply(mstr[[glf$Variable[i]]], .lookup, at))
  }
  return(mstr)
}

get.coll.cite = function(coll, addcoll="", use.initial="x"){
  coll.last = strsplit(coll,",")[[1]][1]
  res = coll.last
  addc.cnt = addcoll
  if(str_detect(addcoll,";")) addc.cnt = str_split(addcoll,";")[[1]]
  if(length(addc.cnt)>0){
    n = length(addc.cnt)
    if(n>1){
      res=paste(res," & et al.",sep="")
    }
    if(n==1){
      if(str_detect(addc.cnt[1],",")){
        add.ln = str_split(addc.cnt[1],",")[[1]][1]
        res = paste(res," & ", add.ln, sep="")
      }
    }
  }
#   if(!is.null(use.initial)){
#     if("x"==use.initial){
#       res=str_replace(",","",coll)
#     }
#   }
  res = str_replace(res," & NA","")
  res
}


load.data <-function(data.file, rules){
  try({
    res = is.exsicFile(data.file, rules = rules)
    #print("\n#1")
    write.csv(res$rule.checks,"rule.checks.csv", fileEncoding = char.enc)
    #print(getwd())
    if(!res$ok) return(res)
#     print(data.file)
#     print(rules)
#     print(list.files())
    try({
      mstr = read.csv(data.file, stringsAsFactors = FALSE, encoding=char.enc)  
    })
    
    #print("#-2")
    #print(mstr)
    try({
      mstr = apply.lookup(mstr)  
    })
    
    #print("#-1")
    conf.file = str_replace(data.file,".csv", "-config.csv")
    conf = read.csv(conf.file, stringsAsFactors = FALSE, encoding=char.enc)
    
    #print("#0")
    #sind = conf[conf$Parameter == "Species-order",]
    #if(sind$Value == "alphabetically"){ # TODO handle more cases
    gsms = paste(mstr$GENUS, mstr$SPECIES)
    gsui = sort(unique(gsms))
    n = length(gsui)
    #print("#1")
    sind = matrix("", nrow=n, ncol=2)
    sind = as.data.frame(sind)
    sind[,1] = as.character(sind[,1])
    sind[,2] = as.character(sind[,2])
    names(sind) = c("GENUS", "SPECIES")
    gs = str_split(gsui, " ")
    for(i in 1:n){
      sind[i, 1] = gs[[i]][1]
      sind[i, 2] = gs[[i]][2]
    }
    #print("#2")
    cntr = conf[conf$Parameter == "Countries",]
    
    # apply filters
    #if(cntr$Countries!="all"){
      uctr = sort(unique(mstr$ORIGCTY))
      mstr <- mstr[mstr$ORIGCTY %in% uctr,]
    #}
     uspc = sort(unique(mstr$SPECIES))
     mstr <- mstr[mstr$SPECIES %in% uspc,]
    
    if(length(which(mstr$remove=="x"))>0) {
      mstr <- mstr[-which(mstr$remove=="x"),]
    }
    #print("#3")
    cite = 1:nrow(mstr)
    spi = rep("",nrow(mstr))
    has.addcoll = "addcoll" %in% names(mstr)
    for(i in 1:nrow(mstr)){
      if(has.addcoll){
        cite[i] = get.coll.cite(mstr[i,"collector"], mstr[i,"addcoll"])
      } else {
        cite[i] = get.coll.cite(mstr[i,"collector"])
      }
        
      sp=mstr[i,"SPECIES"]
      #spi[i] = sind[sind$SPECIES==sp,c("order")]
      spi[i] = paste(str_sub(mstr$GENUS[i],1,1),". ", mstr$SPECIES[i], sep="")
    }
    #print(str(mstr))
    #print(str(spi))
    mstr <- cbind(mstr, cite, spi)
    #print("#5")
    #mstr <- cbind(mstr, cite)
    mstr$cite <- as.character(mstr$cite)
    #print("#6")
    mstr$spi <- as.character(mstr$spi)
    
    #print("before saving")
    #save(mstr,file = "mstr.Rda")
    write.csv(mstr,"mstr.csv", fileEncoding="UTF-8")
    save(cntr,file = "cntr.Rda")
    save(sind,file = "sind.Rda")
    #print("before returning TRUE")
    return(res)
    
  },silent=T)
  #print("some error")
  return(res)
}



format.date=function(d,m=NULL,y=NULL){
  dd="s.d."
  if(d=="--------") return(dd)
  if(!is.na(d)){
    if(is.null(m) & is.null(y)){
      #then assume MCPD format definition for COLLDATE
      dmy = d
      y = str_sub(dmy, 1,4)
      m = str_sub(dmy, 5,6)
      d = str_sub(dmy, 7,8)
      if(m=="00" | m=="--") m = NA
      if(d=="00" | d=="--") d = NA
      d = as.integer(d)
      m = as.integer(m)
      y = as.integer(y)
    }
    
    mm = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
    dd=""
    if(!is.na(d)){
      dd = paste(d,sep="")
    }
    if(!is.na(m)){
      dd = paste(dd," ",mm[m],sep="")
    } else {
      dd=""
    }
    if(!is.na(y)){
      dd = paste(dd," ",y,sep="")
    }
  } else {
    dd=""
  }
  
  dd
}

format.degr=function(lat,lon){
  if(is.na(lat) | is.na(lon)){
    return("n.d.")
  } else {
    latl = ""
    if(lat<0) {
      latl="S"
    } else if(lat>0){
      latl="N"
    }
    latd = abs(round(lat,2))
    latm = round(abs(floor(latd)-latd)*60,0)
    latd = floor(latd)
    lats = paste(latd,"\u00B0",latm,"'",latl,sep="")
    
    lonl = ""
    if(lon<0) {
      lonl="W"
    } else if(lon>0){
      lonl="E"
    }
    lond = abs(round(lon,2))
    lonm = round(abs(floor(lond)-lond)*60,0)
    lond = floor(lond)
    lons = paste(lond,"\u00B0",lonm,"'",lonl,sep="")
    paste(lats,", ",lons,sep="")
  }
}

format.alt=function(alt){
  res=""
  if(!is.na(alt)) {
    res=paste(alt," m",sep="")
  }
  res
}


format.coll=function(coll, num){
  res=paste(coll," ",num,sep="")
  res
}

format.dups=function(dups){
  res=""
  if(!is.na(dups)){
    #sort dups!
    dps=paste(" (",dups,")",sep="")
    if(dps==" ()") dps=""
    res = paste(res,dps, sep="")
  }
  res
}


get.record=function(i, data1, ui=NULL){
  has.addcoll = "addcoll" %in% names(data1)
  has.decgis  = all(c("DECLONGITUDE","DECLATITUDE") %in% names(data1))
  has.elev = "ELEVATION" %in% names(data1)
  has.dups = "dups" %in% names(data1)
  has.site = "COLLSITE" %in% names(data1)
  
  refnms = c("collector","addcoll","number","COLLDATE",
    "ELEVATION","COLLSITE","DECLATITUDE","DECLONGITUDE","dups")
  nms = refnms[refnms %in% names(data1)]
  
  r = data1[data1$id==i, nms]
  if(has.addcoll) {
    cl = get.coll.cite(r[["collector"]], r[["addcoll"]])
  } else {
    cl = get.coll.cite(r[["collector"]])
  }

  dt = format.date(r[["COLLDATE"]])
  
  if(has.decgis){
    dg = format.degr(r[["DECLATITUDE"]],r[["DECLONGITUDE"]])  
  } else dg = ""

  if(has.elev){
    da = format.alt(r[["ELEVATION"]])  
  } else da=""
  
  dc = format.coll(cl,r[["number"]])
  
  if(has.dups){
    dp = format.dups(r[["dups"]])  
  } else dp = ""
  
  if(has.site){
    ds = r["COLLSITE"]
  } else ds=""

  # concatenating
  if(ds!="") {
    rec1 = ds
  } else rec1=""
  
  if(dg!="") {
    if(rec1 != ""){
      rec1 = paste(rec1, dg, sep=", ")  
    } else rec1=dg
  } 
  if(da!="") {
    if(rec1 != ""){
      rec1 = paste(rec1, da, sep=", ")  
    } else rec1=da
  } 
  if(dt!="") {
    if(rec1 != ""){
      rec1 = paste(rec1, dt, sep=", ")  
    } else rec1=dt
  } 
  #rec1=paste(ds,dg,da,dt,sep=", ")
  #rec1=gsub(" ,","",rec1)
  
  rec1=paste(rec1,", *",dc,sep="")
  if(dp != "") {
    rec1=paste(rec1,"* ",str_trim(dp),"",sep="")  
  } else {
    rec1=paste(rec1,"*",sep="")
  }
  
  rec1
}


md2html <-function(file){
  file = str_replace(file,".Rmd",".md")
  txt = readLines(file, encoding = char.enc)
  txt = str_replace(txt,'`','')
  writeLines(txt,file, useBytes=TRUE)
  fn = str_split(file,"\\.")[[1]][1]
  fh = paste(fn,".html",sep="")
  markdownToHTML(file,fh)

  if(is.windows()){
    # tweak output for cross operating system readability
    txt = readLines(fh, encoding = char.enc)
    txt = str_replace(txt,'; charset=utf-8','')
    writeLines(txt,fh, useBytes=TRUE)
  }
  
}

# ' Creates a simple index of countries.
# ' 
# ' Should only be used within the template file.
# ' 
# ' @aliases index.countries
# ' @return a string in markdown format
# ' @author Reinhard Simon
# ' @export
.index.countries <-function # Creates a simple index of countries.
### Creates a simple index of countries as text string.
(){
  ##note<< Should only be used within the template file.
  cntr = NULL
  mstr = NULL
  load("cntr.Rda")
  #load("mstr.Rda")
  mstr = read.csv("mstr.csv", stringsAsFactors = FALSE, encoding = char.enc)
#  if(cntr$ORIGCTY[1]=="all"){
    cn = sort(unique(mstr$ORIGCTY))
#   } else {
#     cn = sort(unique(cntr$ORIGCTY))  
#   }
  
  n = length(cn)
  res=paste(1:n,". ",cn,sep="")
  paste(res,collapse="\n")
  ### A string containing the species formatted 
  ### for the template file in markdown markup.
}

generic.list <-function
(cl="\n") {
  sind = NULL
  load("sind.Rda")
  sp = paste(". *",substr(sind$GENUS,1,1),". ",sind$SPECIES,"*",sep="")
#   sind = NULL
#   load("mstr.Rda")
#   sp = unique(mstr$spi)
  
  n = 1:length(sp)
  x=paste(n,sp, sep="")
  paste(x,collapse=cl)
}

# ' Creates a simple index of species.
# ' 
# ' Should only be used within the template file.
# ' 
# ' @aliases index.species
# ' @return a string in markdown format
# ' @author Reinhard Simon
# ' @export
.index.species <-function # Creates a simple index of species.
### Creates an index of species as text string.
##note<< Should only be used within the template file.
() {
  generic.list()
  ### A string containing the species formatted 
  ### for the template file in markdown markup.
}

# ' Creates the index of collectors and their specimes.
# ' 
# ' Should only be used within the template file.
# ' 
# ' @aliases index.collections
# ' @return a string in markdown format
# ' @author Reinhard Simon
# ' @export
.index.collections <-function # Creates the index of collectors and their specimes.
### Creates the index of collectors and their specimes.
##note<< Should only be used within the template file.
(){
  res=""
  mstr = NULL
  #load("mstr.Rda")
  mstr = read.csv("mstr.csv", stringsAsFactors = FALSE, encoding = char.enc)
  
  if(nrow(mstr)!=0){
    col.ord = sort(unique(mstr$cite))
    #dbs = arrange(mstr,mstr$cite, mstr$number)
    mstr = mstr[order(mstr[,"number"]) ,]
    mstr = mstr[order(mstr[,"cite"]) ,]
    dbs = mstr
    
    for(i in 1:length(col.ord)){
      res = paste(res,col.ord[i]," ",sep="")
      dbt = dbs[dbs$cite==col.ord[i],]
      
      for(k in 1:nrow(dbt)){
        res = paste(res,dbt$number[k]," (",dbt$spi[k],")",sep="")
        if (k<nrow(dbt)){
          res = paste(res,", ",sep="")
        } else {
          res = paste(res,".",sep="")
        }
      }
      res=paste(res,'\n')
      #print(Encoding(res))
      
      
      #res = enc2native(res)
    }
  } else {
    res = "**No records found for this combination of countries and species.**\n"
  }
  Encoding(res)<-char.enc
  #res = enc2native(res)
  res
  ### A string containing the collectors formatted 
  ### for the template file in markdown markup.
  
}

# ' A condensed list of species (listed within a line).
# '
# ' Should only be used within the template file.
# '
# ' @aliases index.species.short
# ' @author Reinhard Simon
# ' @return a string in markdown format
# ' @export
.index.species.short<-function # A condensed list of species (listed within a line).
### A condensed list of species (listed within a line).
##note<< Should only be used within the template file.
(){
  paste("[",generic.list(", "),"]",sep="")
  ### A string containing the species formatted 
  ### for the template file in markdown markup.
}

format.record <- function (db2, res) {
  
  for(l in 1:nrow(db2)){
    res=paste(res,get.record(db2[l,"id"],db2),sep="")
    if(l<nrow(db2)){
      final = "; "
    } else {
      final = "."
    }
    res=paste(res,final,sep="")
  }
  res
}


# ' The list of specimens
# '
# 'List of specimens by species and country detailing collector, ollector number, date and location.
# '
# ' @aliases index.specimens
# ' @author Reinhard Simon
# ' @return a string in markdown format
# ' @export
.index.specimens <-function # The list of specimens 
### List of specimens by species and country detailing collector,
### collector number, date and location.
##note<< Should only be used within the template file.
(){
  sind = NULL 
  load("sind.Rda")
  mstr = NULL
  #load("mstr.Rda")
  mstr = read.csv("mstr.csv", stringsAsFactors = FALSE, encoding = char.enc)
  #mstr = read.csv()
  if(nrow(mstr)==0){
    res = "**No records found for this combination of countries and species.**\n"
  } else {
    res=""  
    for(i in 1:nrow(sind)){
      
      sp = sind[i,]$SPECIES
      db = mstr[mstr$SPECIES==sp,]
      res=paste(res,"### ",i,". ",sind$GENUS[i]," ",sp,"\n",sep="")
      if(nrow(db)==0){
        res = paste(res,
                    "**No records for this species found in this database.**\n"
                    ,sep="")
      } else {
        #sort by collector last name & number
        srt = c("ORIGCTY","majorarea","collector","number")
        # use only those available in database
        nms = names(db)
        srt = srt[srt %in% nms]
        dbs = db[do.call(order,db[srt]),]
        
        
        #filter by countries
        fcntry = sort(unique(dbs$ORIGCTY))
        has.majorarea = "majorarea" %in% names(dbs)
        for(j in 1:length(fcntry)){
          if(nrow(dbs)==0){
            res = paste(res,"**No records for this species in ", fcntry[1]
                        ," in this database.**\n",sep="")
          } else {
            
            db1 = dbs[dbs$ORIGCTY==fcntry[j],]
            res=paste(res,"**",fcntry[j],".** ",sep="")
            
            #filter by admin1
            if(has.majorarea){
              fadm1 = sort(unique(db1$majorarea))
              for(k in 1:length(fadm1)){
                res = paste(res, toupper(fadm1[k]),": ",sep="")
                db2 = db1[db1$majorarea==fadm1[k],]
                db2 = db2[!is.na(db2$id),]
                
                res = format.record(db2, res)
                
                if(k<length(fadm1)){
                  res=paste(res," - ",sep="")
                }
                
              }
              
            } else {
              res = format.record(db1, res)
            }
            res=paste(res,"\n\n",sep="")
          }
        }
        
        
      }
      
    }
  }
  Encoding(res) <- char.enc
  res
  ### A string containing the specimen details formatted 
  ### for the template file in markdown markup.
}

#' Get hyperlink to database of passport data
#' 
#' Creates a hyperlink in markdown format.
#' 
#' @aliases exsic.database
#' @return a string in markdown format with a hyperlink to the local file
#' @author Reinhard Simon
#' @export
exsic.database <- function # Get hyperlink to database of passport data
### Creates a hyperlink in markdown format.
(){
  db.bn = NULL
  load("db.bn.Rda")
  paste("[Database file](",db.bn,")",sep="")
  ### A string containing the hyperlink to the file in markdown format.
}

exsic.metadata = function(){
  res = ""
  try({
    
    res = paste("Computer user:", Sys.getenv("USER"),"\n")
    res = paste(res,"\nTime stamp:",date(),"\n\n")
    if(file.exists("config.csv")){
      config = read.csv("config.csv", stringsAsFactors = FALSE)
      res = paste(res,     "Title:", config[config$Parameter=="Title","Value"],"\n\n")
      res = paste(res, "Author(s):", config[config$Parameter=="Author(s)","Value"],"\n\n")
      res = paste(res,   "Version:", config[config$Parameter=="Version","Value"],"\n\n")
      res = paste(res,  "Comments:", config[config$Parameter=="Comments","Value"],"\n\n")
      res = paste(res,     "Owner:", config[config$Parameter=="Owner","Value"],"\n\n")
      res = paste(res, "Publisher:", config[config$Parameter=="Publisher","Value"],"\n\n")
      res = paste(res,      "Type:", config[config$Parameter=="Type","Value"],"\n\n")
      res = paste(res,    "Format:", config[config$Parameter=="Format","Value"],"\n\n")
      res = paste(res,  "Keywords:", config[config$Parameter=="Keywords","Value"],"\n\n")
      res = paste(res,   "License:", config[config$Parameter=="License","Value"],"\n\n")
      res = paste(res,     "Funding:", config[config$Parameter=="Funding","Value"],"\n\n")
      
      res = paste(res,"\n\n")
      
    }
  })
  return(res)
}




# ' Create specimen indices.
# ' 
# ' Exsic uses markdown to create an html page with exsiccatae indices. 
# ' 
# ' Creates four indices based on passport data of a database of biological specimens.
# ' 
# ' @param data Path to the .csv file containing the database
# ' @param template Path to the template file
# ' @param rules a rules file for checking the data quality
# ' @param lookups a set of file paths used in the rules file in is.oneOf functions
# ' @param config a config file with metadata and parameters
# ' @return void
# ' @author Reinhard Simon
# ' @aliases exsic
# ' @example inst/examples/exsic.R
# ' @export
.exsic <-function (
  data=system.file("samples/exsiccatae.csv",package='exsic'), 
  template=system.file("templates/template-simple.Rmd",package='exsic'), ##<< Path to the template file
  rules  = system.file("samples/rules.R",package='exsic'),
  lookups= c(system.file("samples/countries.csv", package='exsic'), 
             system.file("samples/genus.csv", package='exsic')),
  config = system.file("samples/config.csv", package='exsic')
){
  data.file = data
  template.file = template
  stopifnot(!any(is.null(data.file), is.null(template.file), is.null(rules), is.null(lookups),
           is.null(config)))
  stopifnot(!any(is.na(data.file), is.na(template.file), is.na(rules)))
  stopifnot(!any(data.file=="", template.file=="", rules==""))
  stopifnot(file.exists(data.file))
  stopifnot(file.exists(rules))
  stopifnot(file.exists(template.file))
  
  td = tempdir()
  owd = getwd()
  org.rules = file.path(owd, rules)
  if(is.null(data.file)) {
    data.file = file.choose()
  }

  est = basename(data.file)
  
  pb <- txtProgressBar(0, 100, style=3)
  update.pb(pb, 1, est)
  
  unlink(file.path(td,"*"))
  unlink(file.path(td,"*.*"))
  list.files(td)
  # start copying helper files ...
  #org.rules = rules
  try({
    if(file.exists(rules)) file.copy(rules, file.path(td, "rules.R"), overwrite=TRUE)
    rules = "rules.R"
  })
  try({
    if(file.exists(config)) file.copy(config, file.path(td, "config.csv"), overwrite=TRUE)
    config = "config.csv"
  })
    n = length(lookups)
  #print(lookups)
  #print(n)
    if(lookups[1]!=""){
      for(i in 1:n){
        try({
          file.copy(lookups[i], file.path(td, basename(lookups[i])), overwrite=TRUE)          
        })
        
      }
    }
  
  
  
  adf = basename(data.file)  
  setwd(	td )
  list.files()
  file.copy(file.path(owd,data.file), adf, overwrite = TRUE)
  fr = template.file
  bn = basename(fr)
  tf.to = file.path(getwd(),bn)
  file.copy(fr, tf.to, overwrite=T)
  fr = data.file
  db.bn <- basename(fr)
  save(db.bn, file = "db.bn.Rda")
  db.to = file.path(getwd(),db.bn)
  if(!file.exists(db.to)) file.copy(fr, db.to, overwrite=T)
  
  #print(config)
  try({
    if(file.exists(config)){
      conf = read.csv(config,stringsAsFactors = FALSE)
      write.csv(conf,file="config.csv", row.names = FALSE)
    }
  })
  #print("x")
  kn.to = file.path(getwd(),db.bn)
  kn.to = str_replace(kn.to,".csv",".md")
  hm.to = str_replace(kn.to,".md",".html")

  
  update.pb(pb, 10, est) 
 #print(rules)
#   print(org.rules)
 #print(getwd())
 #print(db.bn)
  #print(file.path(td,rules))
  res = load.data(db.bn, rules = rules ) 
  
  update.pb(pb, 20, est)
  file_check_res = file.path(getwd(), "exsic.check.Rda")
  save(res, file = file_check_res)
  make.errReport(getwd())
  
  ignore.file = c("mstr.csv", "sind.Rda", "cntr.Rda", "db.bn.Rda","error.Rmd",
                  "error_out.md", "exsic.check.Rda", "template-simple.Rmd", basename(kn.to))  
  
  
  if(res$ok) {
    update.pb(pb, 30, est)
    mstr = read.csv("mstr.csv", stringsAsFactors = FALSE)
  
  if(!is.null(mstr)){
  if(interactive()) {update.pb(pb,40, est)}
    knit(tf.to, kn.to )
  if(interactive()) update.pb(pb, 70, est)
    md2html(kn.to)
  } else {
    unlink(hm.to)
    aline = paste("<html><body><p><b>Database file ",db.bn," could not be read.
                  </b></p></body></html>")
    write(aline, hm.to)
    
  }
  } else {
    #print(res)
    te = sum(res$rule.checks$Error.sum)
    warning(paste("The file contains", te,"rule violation(s)."))
  }

  res_dir = basename(data.file)
  res_dir = str_replace(res_dir,".csv","")
  res_dir = paste(res_dir,"__exsic_",pkg.version("exsic"),"_results_",sep="")
  # check potential dir name if exists
  # get all files
  res.set = list.files(path = owd, pattern = res_dir )
  if(length(res.set)==0){
    res_dir = paste(res_dir,"01",sep="")
  } else {
    # extract number and increase; pad with leading zero
    res.set = sort(res.set)
    res.lst = res.set[length(res.set)]
    rs = str_split(res.lst,"_")[[1]]
    rn = length(rs)
    ln = as.integer(rs[rn])
    nn = ln + 1
    nn = str_pad(as.character(nn), width = 2, pad="0")
    res_dir = paste(res_dir,nn,sep="")
  }
  
  nd = file.path(owd,res_dir)
  
  lf = list.files(include.dirs=FALSE)
  lf = lf[!(lf %in% ignore.file)]
  lf = lf[!str_detect(lf,"rs-graphics-")]
  lf = lf[!str_detect(lf,"-config")]
  
  if(!file.exists(nd)) dir.create(nd)
  for(i in 1:length(lf))  {
    nf = file.path(nd,lf[i])
    file.copy(lf[i], nf, overwrite = TRUE)
  }
  if(interactive()) update.pb(pb, 100, est)
  setwd(owd)
  ### Show results page or where it is
  hm.to = file.path(nd, adf)
  hm.to = str_replace(hm.to,".csv", ".html")
  if(interactive()) close(pb)
  if(res$ok){
    message("No data inconsistent with rules found.")
    message(paste("\n\nGo to '", hm.to,"' to review the results."), sep="")   
  } else {
    message(paste("\n\nGo to ",
                  file.path(dirname(hm.to),"error_out.html"),
                  " to review the results."), sep="") 
  }
  
  if(file.exists(hm.to) & interactive() & is.windows()) {shell.exec(hm.to)}
}






