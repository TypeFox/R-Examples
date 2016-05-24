
# library(stringr)
# library(plyr)
# library(XML)

# fips.state=read.csv("./data/FIPS_state_file.txt", sep="|", stringsAsFactors=F)
# fips.county=read.csv("./data/FIPS_county_file.txt", sep=",", stringsAsFactors=F)

## this file needs extra help, due to commas in some of the subdiv
## names, which goofs up the read.csv import
# fips.county.subdivision=read.csv("./data/FIPS_countysubdivision_file.txt", sep=",", stringsAsFactors=F)
# correct some problem with commas in names
# index=nchar(fips.county.subdivision$FUNCSTAT)>1
# fips.county.subdivision$COUSUBNAME[index]=paste(fips.county.subdivision$COUSUBNAME[index], fips.county.subdivision$FUNCSTAT[index], sep=":")
# fips.county.subdivision$FUNCSTAT[index]="N"
# fips.county.subdivision=fips.county.subdivision[!is.na(fips.county.subdivision$STATEFP),]

# fips.place=read.csv("./data/FIPS_place_file.txt", sep="|", stringsAsFactors=F)

# .acs.unit.levels
    #
    # includes all valid types of units for acs estimates
    
    .acs.unit.levels=c("count","dollars","proportion", "ratio", "other")
    
    
    # .acs.dimnames()
    #
    # ensures that a returned acs object includes proper row (geography)
    # and column (col.names) labels
    
    .acs.dimnames=function(acs.obj){
      dimnames(acs.obj@estimate)=list(acs.obj@geography[[1]],acs.obj@acs.colnames)
      dimnames(acs.obj@standard.error)=list(acs.obj@geography[[1]],acs.obj@acs.colnames)
      acs.obj}
    
    # .acs.combine.headers()
    # 
    # create metadata and row/column names for operations that merge two
    # acs units into one
    
    .acs.combine.headers=function(e1, e2, operator) {
      if (endyear(e1)==endyear(e2)) 
        ENDYEAR=endyear(e1)
      else
        ENDYEAR=NA_integer_
      if (currency.year(e1)==currency.year(e2)) 
        CURRENCY.YEAR=currency.year(e1)
      else
        CURRENCY.YEAR=NA_integer_
      if (span(e1)==span(e2))
        SPAN=span(e1)
      else
        SPAN=NA_integer_
      if (identical(geography(e1), geography(e2))) GEOGRAPHY=geography(e1)
      else {
        GEOGRAPHY=geography(e1)
        for (i in 1:length(geography(e1))){
          if (!identical(geography(e1)[[i]], geography(e2)[[i]]))
            GEOGRAPHY[[i]]=paste("(", GEOGRAPHY[[i]], operator, geography(e2)[[i]], ")", sep=" ")
        }
      }
      if (identical(acs.colnames(e1), acs.colnames(e2)))
        {ACS.COLNAMES=acs.colnames(e1)
         ACS.UNITS=acs.units(e1)}
      else
        {ACS.COLNAMES=paste("(", acs.colnames(e1), operator, acs.colnames(e2), ")", sep=" ")
         ACS.UNITS=factor("other", levels=.acs.unit.levels)
       }
      header=list(endyear=ENDYEAR, span=SPAN, currency.year=CURRENCY.YEAR , geography=GEOGRAPHY,
        acs.colnames=ACS.COLNAMES, acs.units=ACS.UNITS)
      header
    }
    
    
    # .acs.identify.units()
    # 
    # used to set units in acs object; initially assumes all units are
    # "counts", then changes some to "dollars" is the word "dollars"
    # matches in colnames.
    
    .acs.identify.units=function(acs.colnames)
    {
      acs.units=rep("count", length(acs.colnames))
      dollar.index=grep(pattern="dollars",x=acs.colnames, fixed=T)
      acs.units[dollar.index]="dollars"
      factor(acs.units, levels=.acs.unit.levels)
    }
    
    
    # .acs.make.constant.object()
    # 
    # use this to create an acs object with some constant value in the
    # estimate and 0 for all the standard errors -- helpful, for example,
    # to add a certain number to every value.
    
    .acs.make.constant.object=function(value, template){
      acs.obj=template
    # if given a vector, replaces by row, not by column
      acs.obj@estimate[,]=value
      acs.obj@standard.error[,]=0
      acs.obj@acs.colnames[1:length(acs.colnames(acs.obj))]=as.character(value)
      acs.obj=.acs.dimnames(acs.obj)
      acs.obj
        }

# an internal function to generate urls variables look like
# c("B21003_001E","B21003_001M") geo.call should have "for" and "in"
# "for" is pairlist, like pairlist(county=05) or c("block+group"="*")
# "in" is pairlist including none or more of: state, county, tract,
# each either number or "*"
api.url.maker=function(endyear, span, key, variables, dataset, geo.call) {
  variables=paste0(variables, collapse=",")
  if(span==0) {span=""}
  api.for=paste0("for=",paste(names(api.for(geo.call)), api.for(geo.call), sep=":"))
  api.in=paste0(paste(names(api.in(geo.call)), api.in(geo.call), sep=":"), collapse="+")
  if (!identical(api.in, "")) api.in=paste0("&in=",api.in)
  api.url=paste0("http://api.census.gov/data/", endyear, "/" , dataset, span, "?key=", key, "&get=", variables, ",NAME&", api.for, api.in)
  api.url
}

# a function to install a users key for use in this and future
# sessions
# writes to package extdata directory

api.key.install=function(key, file="key.rda") {
    dir.create(paste(system.file(package="acs"), "extdata", sep="/"), showWarnings=F)
    save(key, file=paste(system.file("extdata", package="acs"), file, sep="/"))
}

# a function to migrate a previously installed api.key after package update
# writes to package extdata directory

api.key.migrate=function() {
  key.path=system.file("../key-holder.rda", package="acs")
  if(file.exists(key.path)) {
      dir.create(paste(system.file(package="acs"), "extdata", sep="/"), showWarnings=F)
      file.copy(from=key.path, to=paste(system.file("extdata", package="acs"), "key.rda", sep="/"), overwrite=T)
  } else {warning("No archived key found;\n  try api.key.install with new key.")}
}
 

# a function to download XML variable lookup files,
# to speed up future acs.lookup calls.

acs.tables.install=function() {
    filelist=readLines("http://web.mit.edu/eglenn/www/acs/acs-variables/filelist.txt")
    dir.create(paste(system.file(package="acs"), "extdata", sep="/"), showWarnings=F)
    for (i in filelist) {
        download.file(url=paste("http://web.mit.edu/eglenn/www/acs/acs-variables/", i, sep=""), destfile=paste(system.file("extdata/", package="acs"), i, sep=""))}
}

setClass(Class="acs", representation =
  representation(endyear="integer", span="integer",
  geography="data.frame", acs.colnames="character", modified="logical"
  , acs.units="factor", currency.year="integer", estimate="matrix",
  standard.error="matrix"), prototype(endyear=NA_integer_,
  span=NA_integer_, acs.units=factor(levels=.acs.unit.levels),
  currency.year=NA_integer_, modified=F))
  
is.acs=function(object){
  if (class(object)=="acs") {TRUE}
  else{FALSE}}

globalVariables(c("fips.state","fips.school","fips.county.subdivision", "fips.american.indian.area","fips.county", "fips.place"))


  setClass(Class="geo", representation =
    representation(api.for="list", api.in="list", name="character", sumlev="numeric"), prototype())
  
  is.geo=function(object){
    if (class(object)=="geo") {TRUE}
    else{FALSE}}
  
  if (!isGeneric("api.for")) {
    setGeneric("api.for", def=function(object){standardGeneric("api.for")})}else{}
  setMethod("api.for", "geo", function(object) object@api.for)
  
  if (!isGeneric("api.in")) {
    setGeneric("api.in", def=function(object){standardGeneric("api.in")})}else{}
  setMethod("api.in", "geo", function(object) object@api.in)
  
  if (!isGeneric("name")) {
    setGeneric("name", def=function(object){standardGeneric("name")})}else{}
  setMethod("name", "geo", function(object) object@name)
  
  if (!isGeneric("sumlev")) {
    setGeneric("sumlev", def=function(object){standardGeneric("sumlev")})}else{}
  setMethod("sumlev", "geo", function(object) object@sumlev)
  
  setClass(Class="geo.set", representation =
           representation(geo.list="list", combine="logical", combine.term="character"),
           prototype(combine=T, combine.term="aggregate"))
  is.geo.set=function(object){
    if (class(object)=="geo.set") {TRUE}
    else{FALSE}}
  if (!isGeneric("combine")) {
    setGeneric("combine", def=function(object){standardGeneric("combine")})}else{}
  setMethod("combine", "geo.set", function(object) object@combine)
  if (!isGeneric("combine.term")) {
    setGeneric("combine.term", def=function(object){standardGeneric("combine.term")})}else{}
  setMethod("combine.term", "geo.set", function(object) object@combine.term)
  if (!isGeneric("geo.list")) {
    setGeneric("geo.list", def=function(object){standardGeneric("geo.list")})}else{}
  
  setMethod("geo.list", "geo.set", function(object) {
    if(length(object@geo.list)==1) object@geo.list[[1]]
    else object@geo.list})
  
  setMethod("show", signature(object = "geo"), function(object) {
    cat("\"geo\" object: ")
    print(object@name)
  })
    
  # setMethod("show", signature(object = "geo.set"), function(object) {
  #   cat("An object of class \"geo.list\"\n\n")
  #   cat("Slot \"geo.list\":\n")
  #   print(lapply(X=object@geo.list, FUN=name))
  #   cat("Slot \"combine\":\n")
  #   print(object@combine)
  #   cat("Slot \"combine.term\":\n")
  #   print(object@combine.term)
  # })
    
    
  # method to combine geos and geo.sets
  setMethod("+", signature(e1 = "geo", e2 = "geo"), function(e1, e2) {
    geo.set.obj=new(Class="geo.set",
      combine=FALSE,
      geo.list=c(e1, e2))
    geo.set.obj
  })
  setMethod("+", signature(e1 = "geo.set", e2 = "geo"), function(e1, e2) {
    geo.set.obj=new(Class="geo.set",
      combine=combine(e1),
      geo.list=c(geo.list(e1), e2))
    geo.set.obj
  })
  setMethod("+", signature(e1 = "geo", e2 = "geo.set"), function(e1, e2) {
    geo.set.obj=new(Class="geo.set",
      combine=combine(e2),
      geo.list=c(e1, geo.list(e2)))
    geo.set.obj
  })
  
  # Note on "+"
  # if geo.list(e1), geo.list(e2) are geo.sets, then take the geolists of them, recusively; 
  # should always yield flattend set.
  
  setMethod("+", signature(e1 = "geo.set", e2 = "geo.set"), function(e1, e2) {
    if(is.geo(geo.list(e2))) combine=combine(e1) else 
    if(is.geo(geo.list(e1))) combine=combine(e2) else
    combine=F
    geo.set.obj=flatten.geo.set(c(e1,e2))
    combine(geo.set.obj)=combine
    combine.term(geo.set.obj)=paste(combine.term(e1), " + ", combine.term(e2))
    geo.set.obj
  })
  
  flatten.geo.set=function(x) {
    if (!is.geo.set(x)) return(NA)
    if (length(x)==1){
      if (is.geo(geo.list(x))) {
        return(x)}
      if (is.geo.set(geo.list(x))) {
        return(flatten.geo.set(geo.list(x[1])))}
    }
    if (length(x)>1){
      a=flatten.geo.set(x[1])
      b=flatten.geo.set(x[2:length(x)])
      new(Class="geo.set", combine=combine(x),
    combine.term=combine.term(x), geo.list=c(geo.list(a), geo.list(b)))
    }
  }
        
  setMethod("c", signature(x = "geo.set"), function(x, y, ..., combine=F, combine.term="aggregate", recursive = FALSE) {
    if (recursive) {
      if (missing(y)) geo.set.obj=x
      else geo.set.obj=x+c(y, ..., recursive=T)
      geo.set.obj@combine=combine
      geo.set.obj@combine.term=combine.term
    } else {
      if (missing(y)) {geo.set.obj=x} else {
        if (length(y)==1) {
          geo.set.obj=c((x+geo.list(y)),...)
          geo.set.obj@combine=combine
          geo.set.obj@combine.term=combine.term        
        }
        else {
          geo.set.obj=new(Class="geo.set", combine=combine, combine.term=combine.term, geo.list=list(x,y, ...))}}
    }
    geo.set.obj
  })
         
  # NOTE: changed this to prevent extra nesting when only one
  # geo.set is subsetted
  
  setMethod(f="[", signature="geo.set", definition=function(x,i,j,...,drop=FALSE){
    if (missing(i)) i=j
    if (missing(j)) j=i
    if (length(i)==1 && is.geo.set(x@geo.list[[i]])) return(x@geo.list[[i]])
    else
      new(Class="geo.set",
          geo.list=x@geo.list[i], combine=combine(x),
          combine.term=paste(combine.term(x), "(partial)", sep=" "))})
  
  setMethod(f="[[", signature="geo.set", definition=function(x,i,j,...,drop=FALSE){
    if (missing(i)) i=j
    if (missing(j)) j=i
    x@geo.list[[i]]})
  
  # need to work on to allow to change combine values -- seem to not like when you replace more than one
  setReplaceMethod(f="[", signature="geo.set",
    definition=function(x,i,j,value){
    if (missing(i)) i=j
    if (missing(j)) j=i
    if (length(i)==1){
      x@geo.list[i]=value
      validObject(x)
      return(x)}
    else {
      for (a in i) {
        x@geo.list[i]=value@geo.list[a]
      }
      validObject(x)
      return(x)
    }
  })
  
  setReplaceMethod(f="[[", signature="geo.set",
    definition=function(x,i,j,value){
    if (missing(i)) i=j
    if (missing(j)) j=i
    x@geo.list[i]=value
    validObject(x)
    return(x)})
    
  
  ##  for some reason, can't declare new generic for "length"
  # if (!isGeneric("length")) {
  #  setGeneric("length", def=function(x){standardGeneric("length")})}else{}

  setMethod("length", "geo.set", function(x) length(x@geo.list))
  
  if (!isGeneric("combine<-")) {
            setGeneric("combine<-", def=function(object, value){standardGeneric("combine<-")})}else{}
  
  setReplaceMethod(f="combine", signature="geo.set",
                           definition=function(object,value){
                             object@combine=value
                             validObject(object)
                             return(object)
                           })
  
  if (!isGeneric("combine.term<-")) {
            setGeneric("combine.term<-", def=function(object, value){standardGeneric("combine.term<-")})}else{}
  
  setReplaceMethod(f="combine.term", signature="geo.set",
                           definition=function(object,value){
                             object@combine.term=value
                             validObject(object)
                             return(object)
                           })
  
  geo.lookup=function(state, county, county.subdivision, place, american.indian.area, school.district, school.district.elementary, school.district.secondary, school.district.unified) {
    # first deal with american indian areas -- only one with no state
    if (!missing(american.indian.area)) {
      if (nargs()>1) {warning("american indian area selected; no other options allowed; \n  returning NA")
                      return(NA)}
      if(is.character(american.indian.area)){
        fips.american.indian.area=fips.american.indian.area[grepl(paste(american.indian.area, collapse="|"), fips.american.indian.area$American.Indian.Area.Name), ]}
      else {
        fips.american.indian.area=fips.american.indian.area[fips.american.indian.area$American.Indian.Area.Code %in% american.indian.area]}
      american.indian.area=fips.american.indian.area[,1]
      american.indian.area.name=fips.american.indian.area[,2]    
      results=data.frame(american.indian.area, american.indian.area.name, stringsAsFactors=F)
      return(results)
    }
    # all remaining need state
    state.name=NA
    if (missing(state)){warning("state required for geo.lookup with these options; \n  returning NA")
                        return(NA)}
    for (i in 1:length(state)){
      if (is.character(state[i])){
        if (nchar(state[i])==2) {
          state[i]=fips.state[fips.state$STUSAB==state[i],1]}
        else {
          state[i]=fips.state[grep(paste("^", state[i], sep=""), fips.state$STATE_NAME),1]}
      }
    }
    state=state[state %in% fips.state$STATE] # remove non-matches
    state.name=fips.state[fips.state$STATE %in% state ,3]
    if (length(state)==0) {
      warning("No valid state names match search string;\n  returning NA")
      return(NA)}
    if(length(state) > 1){
      state=as.integer(state)
      a=geo.lookup(state=state[1], county=county, county.subdivision=county.subdivision, place=place, school.district=school.district, school.district.elementary=school.district.elementary, school.district.secondary=school.district.secondary, school.district.unified=school.district.unified)
      b=geo.lookup(state=state[2:length(state)], county=county, county.subdivision=county.subdivision, place=place, school.district=school.district, school.district.elementary=school.district.elementary, school.district.secondary=school.district.secondary, school.district.unified=school.district.unified)
      return(rbind.fill(a,b))}
    results=data.frame(state=state, state.name=state.name, stringsAsFactors=F)
    # check counties
    fips.county.sub=fips.county.subdivision[fips.county.subdivision$STATEFP==state,]
    if (!missing(county)) {
      fips.county=fips.county[fips.county$State.ANSI==state,]
      if(is.character(county)){
        fips.county=fips.county[grepl(paste(county, collapse="|"), fips.county$County.Name), ]}
      else {
        fips.county=fips.county[fips.county$County.ANSI %in% county, ]
      }
      county=fips.county[,3]
      # need to fix for when county is numeric vector, here and below
  #    else {county=county[county %in% fips.county$County.ANSI]} 
      county.name=fips.county[, 4]
      if(length(county)>0){
        results=rbind.fill(results, data.frame(state, state.name, county, county.name, stringsAsFactors=F))}}
    # check subdivisions, when no county given; state still required
    if (missing(county) && !missing(county.subdivision)){
      if(is.character(county.subdivision)){
        fips.county.sub=fips.county.sub[grepl(paste(county.subdivision,
      collapse="|"), fips.county.sub$COUSUBNAME) , ]}
      else {
        fips.county.sub=fips.county.sub[fips.county.sub$COUSUBFP %in% county.subdivision , ]}
      county.subdivision=fips.county.sub[,5]
      subdivision.name=fips.county.sub[ , 6]
      this.county=fips.county.sub[ , 3]
      this.county.name=fips.county.sub[ , 4]
      if(length(county.subdivision)>0){
        results=rbind.fill(results, data.frame(state, state.name, county=this.county, county.name=this.county.name, county.subdivision, county.subdivision.name=subdivision.name, stringsAsFactors=F))}}
    # check subdivisions, when county is given
    if (!missing(county) && !missing(county.subdivision)){
      if(is.character(county.subdivision)){
        fips.county.sub=fips.county.sub[grepl(paste(county.subdivision, collapse="|"), fips.county.sub$COUSUBNAME) & fips.county.sub$COUNTYFP %in% county, ]}
      else {
        fips.county.sub=fips.county.sub[fips.county.sub$COUSUBFP %in% county.subdivision & fips.county.sub$COUNTYFP %in% county ,]}
      county.subdivision=fips.county.sub[,5]
      subdivision.name=fips.county.sub[ , 6]
      this.county=fips.county.sub[ , 3]
      this.county.name=fips.county.sub[ , 4]
      if(length(county.subdivision)>0){
        results=rbind.fill(results, data.frame(state, state.name, county=this.county, county.name=this.county.name, county.subdivision, county.subdivision.name=subdivision.name, stringsAsFactors=F))}}
    # check place
    if (!missing(place)){
        fips.place=fips.place[fips.place$STATEFP==state,]    
        if(is.character(place)) {
          fips.place=fips.place[grepl(paste(place, collapse="|"), fips.place$PLACENAME), ]}
        else fips.place=fips.place[fips.place$PLACEFP %in% place,]
        place=fips.place[,3]
        place.name=fips.place[, 4]
        this.county.name=fips.place[, 7]
        if(length(place)>0){
          results=rbind.fill(results, data.frame(state, state.name, county.name=this.county.name, place, place.name, stringsAsFactors=F))}}
    # check schools
    ## elementary
    if (!missing(school.district.elementary)){
        fips.school.elementary=fips.school[fips.school$STATEFP==state & fips.school$TYPE=="Elementary",]    
        if(is.character(school.district.elementary)) {
          fips.school.elementary=fips.school.elementary[grepl(paste(school.district.elementary, collapse="|"), fips.school.elementary$SDNAME), ]}
        else fips.school.elementary=fips.school.elementary[fips.school.elementary$LEA %in% school.district.elementary,]
        school.district.elementary=fips.school.elementary[,3] # fips code
        school.district.elementary.name=fips.school.elementary[, 4] # name
        school.district.elementary.type=fips.school.elementary[, 5] # type (elem, secondary, unified)
        if(length(school.district.elementary)>0){
          results=rbind.fill(results, data.frame(state, state.name, school.district.elementary, school.district.elementary.name, school.district.elementary.type, stringsAsFactors=F))}}
    ## secondary
    if (!missing(school.district.secondary)){
        fips.school.secondary=fips.school[fips.school$STATEFP==state  & fips.school$TYPE=="Secondary",]    
        if(is.character(school.district.secondary)) {
          fips.school.secondary=fips.school.secondary[grepl(paste(school.district.secondary, collapse="|"), fips.school.secondary$SDNAME), ]}
        else fips.school.secondary=fips.school.secondary[fips.school.secondary$LEA %in% school.district.secondary,]
        school.district.secondary=fips.school.secondary[,3] # fips code
        school.district.secondary.name=fips.school.secondary[, 4] # name
        school.district.secondary.type=fips.school.secondary[, 5] # type (elem, secondary, unified)
        if(length(school.district.secondary)>0){
          results=rbind.fill(results, data.frame(state, state.name, school.district.secondary, school.district.secondary.name, school.district.secondary.type, stringsAsFactors=F))}}
    ## unified
    if (!missing(school.district.unified)){
        fips.school.unified=fips.school[fips.school$STATEFP==state  & fips.school$TYPE=="Unified",]    
        if(is.character(school.district.unified)) {
          fips.school.unified=fips.school.unified[grepl(paste(school.district.unified, collapse="|"), fips.school.unified$SDNAME), ]}
        else fips.school.unified=fips.school.unified[fips.school.unified$LEA %in% school.district.unified,]
        school.district.unified=fips.school.unified[,3] # fips code
        school.district.unified.name=fips.school.unified[, 4] # name
        school.district.unified.type=fips.school.unified[, 5] # type (elem, secondary, unified)
        if(length(school.district.unified)>0){
          results=rbind.fill(results, data.frame(state, state.name, school.district.unified, school.district.unified.name, school.district.unified.type, stringsAsFactors=F))}}
    ## any type
    if (!missing(school.district)){
        fips.school.any=fips.school[fips.school$STATEFP==state,]    
        if(is.character(school.district)) {
          fips.school.any=fips.school.any[grepl(paste(school.district, collapse="|"), fips.school.any$SDNAME), ]}
        else fips.school.any=fips.school.any[fips.school.any$LEA %in% school.district,]
        school.district=fips.school.any[,3] # fips code
        school.district.name=fips.school.any[, 4] # name
        school.district.type=fips.school.any[, 5] # type (elem, secondary, unified)
        if(length(school.district)>0){
          results=rbind.fill(results, data.frame(state, state.name, school.district, school.district.name, school.district.type, stringsAsFactors=F))}}
    results
  }
  
  
  geo.make=function(us, region, division, state, county, county.subdivision, place, tract,
    block.group, msa, csa, necta, urban.area, congressional.district, state.legislative.district.upper, state.legislative.district.lower, puma, zip.code, american.indian.area, school.district.elementary, school.district.secondary, school.district.unified, combine=F, combine.term="aggregate", check=FALSE, key="auto") {
    .geo.unit.make=function(us, region, division, state, county, county.subdivision, place, tract,
    block.group, msa, csa, necta, urban.area, congressional.district, state.legislative.district.upper, state.legislative.district.lower, puma, zip.code, american.indian.area, school.district.elementary, school.district.secondary, school.district.unified) {
      geo.obj=NA
      nargs=nargs()
      ## geos with only one argument: sumlev 010, 020, 030, 350, 400, 860
      # sumlev 010 -- all US
      if (!missing(us) && (us==1 || us=="*" || us==TRUE)) {
        if (nargs>1) {warning("entire U.S. selected; no other geographic arguments allowed;\n  returning NA")
                      return(NA)}
        if (us==TRUE) us=1
        geo.obj=new(Class="geo", api.for=list("us"=us), api.in=list(), name="US", sumlev=10)
        return(geo.obj)
      }
      # sumlev 020 -- region
      if (!missing(region)) {
        if (nargs>1) {warning("region selected; ; no other geographic arguments allowed;\n  returning NA")
                      return(NA)}
        geo.obj=new(Class="geo", api.for=list("region"=region), api.in=list(), name=paste("Region ", region, sep=""), sumlev=20)
        return(geo.obj)
      }
      # sumlev 030 -- division
      if (!missing(division)) {
        if (nargs>1) {warning("division selected; no other geographic arguments allowed;\n  returning NA")
                      return(NA)}
        geo.obj=new(Class="geo", api.for=list("division"=division), api.in=list(), name=paste("Division ", division, sep=""), sumlev=30)
        return(geo.obj)
      }
      # sumlev 250 -- american indian area/alaska native area/hawaiian home land
      ## check for american indian area names and convert to codes
      if (!missing(american.indian.area) && is.character(american.indian.area) && american.indian.area !="*"){
        american.indian.area=fips.american.indian.area[grepl(paste("^", american.indian.area, sep=""), fips.american.indian.area$American.Indian.Area.Name), 1]
        if (length(american.indian.area)>1){
          warning("More than one American Indian area/Alaska Native area/Hawaiian Home Land name matches search string;\n  returning NA;\n  Perhaps try 'geo.lookup()...?")
          return(NA)}
        if (length(american.indian.area)==0){
          warning("No American Indian area/Alaska Native area/Hawaiian Home Land name matches search string;\n  returning NA;\n  Perhaps try 'geo.lookup()...?")
          return(NA)}
      }
      ## deal with wildcards
      if (hasArg(american.indian.area) && american.indian.area!="*")
        american.indian.area.name=fips.american.indian.area[fips.american.indian.area$American.Indian.Area.Code==american.indian.area, 2]
      else american.indian.area.name="All American Indian/Alaska Native/Hawaiian Home Land areas"  # later, not used unless hasArg(county)
      ## make geo for sumlev 250
      if (!missing(american.indian.area)) {
        if (nargs>1) {warning("american.indian.area selected; \n  no other geographic arguments allowed;\n  returning NA")
                      return(NA)}
        geo.obj=new(Class="geo", api.for=list("american+indian+area/alaska+native+area/hawaiian+home+land"=american.indian.area), api.in=list(), name=american.indian.area.name, sumlev=250)
        return(geo.obj)
      }
      # sumlev 310 -- metropolitan statistical area/micropolitan statistical area **NO STATE**
      if (!missing(msa) && missing(state)) {
        if (nargs>1) {warning("msa selected; ignoring all other geographic arguments")}
        geo.obj=new(Class="geo", api.for=list("metropolitan+statistical+area/micropolitan+statistical+area"=msa), api.in=list(), name=paste("MSA ", msa, sep=""), sumlev=310)
        return(geo.obj)
      }
      # sumlev 330 -- combined statistical area **NO STATE**
      if (!missing(csa) && missing(state)) {
        if (nargs>1) {warning("csa selected; ignoring all other geographic arguments")}
        geo.obj=new(Class="geo", api.for=list("combined+statistical+area"=csa), api.in=list(), name=paste("CSA ", csa, sep=""), sumlev=330)
        return(geo.obj)
      }
      # sumlev 350 -- new england city and town area
      if (!missing(necta)) {
        if (nargs>1) {warning("necta selected;  no other geographic arguments allowed;\n  returning NA")
                      return(NA)}
        geo.obj=new(Class="geo", api.for=list("new+england+city+and+town+area"=necta), api.in=list(), name=paste("NECTA ", necta, sep=""), sumlev=350)
        return(geo.obj)
      }
      # sumlev 400 -- urban area
      if (!missing(urban.area)) {
        if (nargs>1) {warning("urban.area selected;  no other geographic arguments allowed;\n  returning NA")
                      return(NA)}
        geo.obj=new(Class="geo", api.for=list("urban+area"=urban.area), api.in=list(), name=paste("Urban Area ", urban.area, sep=""), sumlev=400)
        return(geo.obj)
      }
      # sumlev 860 -- zip.code
      if (!missing(zip.code)) {
        if (nargs>1) {warning("zip.code selected, no other geographic arguments allowed;\n  returning NA")
                    return(NA)}
        geo.obj=new(Class="geo", api.for=list("zip+code+tabulation+area"=zip.code), api.in=list(), name=paste("Zip Code Tabulation Area ", zip.code, sep=""), sumlev=860)
        return(geo.obj)
      }
      # all other geos need a valid state 
      if(missing(state)) {
        warning("state required")
        return(NA)}
      # check for state abbreviations/text names and convert to codes
      if (is.character(state) && state !="*"){
        if (nchar(state)==2) state=fips.state[fips.state$STUSAB==state,1]
        else state=fips.state[grep(paste("^", state, sep=""), fips.state$STATE_NAME),1]
        if (length(state)==0) {
          warning("No valid state names match search string;\n  returning NA")
          return(NA)}
        if (length(state)>1) {
          warning("More than one state name matches search string;\n  returning NA;\n  Perhaps try 'geo.lookup()...?")
          return(NA)}
      }
      # check for county names and convert to codes
      if (!missing(county) && is.character(county) && county !="*"){
        county=fips.county[grepl(paste("^", county, sep=""), fips.county$County.Name) & (fips.county$State.ANSI==state), 3]
        if (length(county)>1){
          warning("More than one county name matches search string in state ", state,";\n  returning NA;\n  Perhaps try 'geo.lookup()...?")
          return(NA)}
        if (length(county)==0){
          warning("No county name matches search string in state ", state,";\n  returning NA;\n  Perhaps try 'geo.lookup()...?")
          return(NA)}
      }
      # check for county subdivision names and convert to codes
     if (!missing(county.subdivision) && is.character(county.subdivision) && county.subdivision !="*"){
        county.subdivision=fips.county.subdivision[grepl(paste("^", county.subdivision, sep=""), fips.county.subdivision$COUSUBNAME) & (fips.county.subdivision$STATEFP==state) , 5]
        if (length(county.subdivision)>1){
          warning("More than one county subdivision name matches search string in state ", state,";\n  returning NA;\n  Perhaps try 'geo.lookup()...?")
          return(NA)}
        if (length(county.subdivision)==0){
          warning("No county subdivision name matches search string in state ", state,";\n  returning NA;\n  Perhaps try 'geo.lookup()...?")
          return(NA)}
      }
      # check for place names and convert to codes
      if (!missing(place) && is.character(place) && place !="*"){
        place=fips.place[grepl(paste("^", place, sep=""), fips.place$PLACENAME) & (fips.place$STATEFP==state), 3]
        if (length(place)>1){
          if(isTRUE(all.equal(max(place),min(place)))
             && isTRUE(all.equal(max(place), min(place)))) # i.e., all match same place #
            {
              place=place[1]
            } else {
              warning("More than one place name matches search string in state ", state,";\n  returning NA;\n  Perhaps try 'geo.lookup()...?")
              return(NA)}
        }
        if (length(place)==0){
          warning("No place name matches search string in state ", state,";\n  returning NA;\n  Perhaps try 'geo.lookup()...?")
          return(NA)}
      }
      # check for school district names and convert to codes
      ## elementary
      if (!missing(school.district.elementary) && is.character(school.district.elementary) && school.district.elementary !="*"){
        school.district.elementary=fips.school[grepl(paste("^", school.district.elementary, sep=""), fips.school$SDNAME) & (fips.school$STATEFP==state) & (fips.school$TYPE=="Elementary"), 3]
        if (length(school.district.elementary)>1){
          warning("More than one elementary school district name matches search string in state ", state,";\n  returning NA;\n  Perhaps try 'geo.lookup()...?")
          return(NA)}
        if (length(school.district.elementary)==0){
          warning("No elementary school district name matches search string in state ", state,";\n  returning NA;\n  Perhaps try 'geo.lookup()...?")
          return(NA)}
      }
      ## secondary
      if (!missing(school.district.secondary) && is.character(school.district.secondary) && school.district.secondary !="*"){
        school.district.secondary=fips.school[grepl(paste("^", school.district.secondary, sep=""), fips.school$SDNAME) & (fips.school$STATEFP==state) & (fips.school$TYPE=="Secondary"), 3]
        if (length(school.district.secondary)>1){
          warning("More than one secondary school district name matches search string in state ", state,";\n  returning NA;\n  Perhaps try 'geo.lookup()...?")
          return(NA)}
        if (length(school.district.secondary)==0){
          warning("No secondary school district name matches search string in state ", state,";\n  returning NA;\n  Perhaps try 'geo.lookup()...?")
          return(NA)}
      }
      ## unified
      if (!missing(school.district.unified) && is.character(school.district.unified) && school.district.unified !="*"){
        school.district.unified=fips.school[grepl(paste("^", school.district.unified, sep=""), fips.school$SDNAME) & (fips.school$STATEFP==state) & (fips.school$TYPE=="Unified"), 3]
        if (length(school.district.unified)>1){
          warning("More than one unified school district name matches search string in state ", state,";\n  returning NA;\n  Perhaps try 'geo.lookup()...?")
          return(NA)}
        if (length(school.district.unified)==0){
          warning("No unified school district name matches search string in state ", state,";\n  returning NA;\n  Perhaps try 'geo.lookup()...?")
          return(NA)}
      }
      # deal with names for wildcards
      if (state != "*") state.name=fips.state[fips.state$STATE==state,3] else state.name="All states"
      if (hasArg(county) && county !="*")
        county.name=fips.county[fips.county$County.ANSI==county & fips.county$State.ANSI==state, 4]
      else county.name="All counties"  # later, not used unless hasArg(county) 
      if (hasArg(place)  && place !="*") {
        place.name=fips.place[fips.place$PLACEFP==place &
          fips.place$STATEFP==state, 4]
        if(length(place.name)>1) place.name=place.name[1]}
      else place.name="All places"  # later, not used unless hasArg(place)
      if (hasArg(county.subdivision) && county.subdivision !="*")
        subdivision.name=fips.county.subdivision[fips.county.subdivision$COUSUBFP==county.subdivision & fips.county.subdivision$STATEFP==state, 6]
      else subdivision.name="All subdivisions" # later, not used unless hasArg(county.subdivision)
      if (hasArg(school.district.elementary) && school.district.elementary!="*")
        school.district.elementary.name=fips.school[fips.school$LEA==school.district.elementary & fips.school$STATEFP==state, 4]
      else school.district.elementary.name="All elementary school districts"  # later, not used unless hasArg(county)
      if (hasArg(school.district.secondary) && school.district.secondary!="*")
        school.district.secondary.name=fips.school[fips.school$LEA==school.district.secondary & fips.school$STATEFP==state, 4]
      else school.district.secondary.name="All secondary school districts"  # later, not used unless hasArg(county)
      if (hasArg(school.district.unified) && school.district.unified!="*")
        school.district.unified.name=fips.school[fips.school$LEA==school.district.unified & fips.school$STATEFP==state, 4]
      else school.district.unified.name="All unified school districts"  # later, not used unless hasArg(county)
      # done with looking up names and codes
      # sumlev 40 (state)
      if (nargs==1 && hasArg(state)) 
        geo.obj=new(Class="geo", api.for=list(state=state), api.in=list(), name=state.name, sumlev=40)
      # sumlev 50 (state-county)
      if (nargs==2 && hasArg(county))
        geo.obj=new(Class="geo", api.for=list(county=county), api.in=list(state=state), name=paste(county.name, state.name, sep=", "), sumlev=50)
      # sumlev 160 (state-place) -- only sumlev with place
      if (hasArg(place)) {
        if(nargs>2) warning("Using sumlev 160 (state-place)\n  Other levels not supported by census api at this time")
        geo.obj=new(Class="geo", api.for=list(place=place), api.in=list(state=state), name=paste(place.name, state.name, sep=", "), sumlev=160) }
      # sumlev 60 (state-county-county.subdivision)
      if (nargs==3 && hasArg(county) && hasArg(county.subdivision))
        geo.obj=new(Class="geo", api.for=list("county+subdivision"=county.subdivision), api.in=list(state=state, county=county), name=paste(subdivision.name, county.name, state.name, sep=", "), sumlev=60)
      # sumlev 140 (state-county-tract)
      if (nargs==3 && hasArg(county) && hasArg(tract))
        geo.obj=new(Class="geo", api.for=list(tract=tract), api.in=list(state=state, county=county), name=paste("Tract ", tract, ", ", county.name, ", ", state.name, sep=""), sumlev=140)
      # sumlev 150 (state-county-tract-block.group)
      if (nargs==4 && hasArg(county) && hasArg(tract) && hasArg(block.group))
        geo.obj=new(Class="geo", api.for=list("block+group"=block.group), api.in=list(state=state, county=county, tract=tract), name=paste("Tract ", tract, ", Blockgroup ", block.group, ", ", county.name, ", ", state.name, sep=""), sumlev=150)
      # sumlev 320 (state-msa)
      if (nargs==2 && hasArg(msa))
       geo.obj=new(Class="geo", api.for=list("metropolitan+statistical+area/micropolitan+statistical+area"=msa), api.in=list(state=state), name=paste("MSA ", msa, ", ", state.name, sep=""), sumlev=320)
      # sumlev 340 (state-csa)
      if (nargs==2 && hasArg(csa))
       geo.obj=new(Class="geo", api.for=list("combined+statistical+area"=csa), api.in=list(state=state), name=paste("CSA ", csa, ", ", state.name, sep=""), sumlev=340)
      # sumlev 500 (state-congressional.district)
      if (nargs==2 && hasArg(congressional.district))
       geo.obj=new(Class="geo", api.for=list("congressional+district"=congressional.district), api.in=list(state=state), name=paste("Congressional District ", congressional.district, ", ", state.name, sep=""), sumlev=500)
      # sumlev 510 (state-congressional.district-county)
      if (nargs==3 && hasArg(county) && hasArg(congressional.district))
       geo.obj=new(Class="geo", api.for=list("county"=county), api.in=list(state=state, "congressional+district"=congressional.district), name=paste(county.name, ", Congressional District ", congressional.district, ", ", state.name, sep=""), sumlev=510)
      # sumlev 610 (state-state.leg.upper)
      if (nargs==2 && hasArg(state.legislative.district.upper)) {
        # turn to string exactly three characters long
        if (state.legislative.district.upper !="*") state.legislative.district.upper=str_sub(paste("000",state.legislative.district.upper,sep=""), -3,-1)
        geo.obj=new(Class="geo", api.for=list("state+legislative+district+(upper+chamber)"=state.legislative.district.upper), api.in=list(state=state), name=paste("State Legislative District (upper chamber) ", state.legislative.district.upper, ", ", state.name, sep=""), sumlev=610)}
      # sumlev 620 (state-state.leg.lower)
      if (nargs==2 && hasArg(state.legislative.district.lower)) {
        # turn to string exactly three characters long
        if (state.legislative.district.lower !="*")  state.legislative.district.lower=str_sub(paste("000",state.legislative.district.lower,sep=""), -3,-1)
        geo.obj=new(Class="geo", api.for=list("state+legislative+district+(lower+chamber)"=state.legislative.district.lower), api.in=list(state=state), name=paste("State Legislative District (lower chamber) ", state.legislative.district.lower, ", ", state.name, sep=""), sumlev=620)}
      # sumlev 795 (state-puma)
      if (nargs==2 && hasArg(puma))
       geo.obj=new(Class="geo", api.for=list("public+use+microdata+area"=puma), api.in=list(state=state), name=paste("Public Use Microdata Area ", puma, ", ", state.name, sep=""), sumlev=795)
      # sumlev 950 (state-school.elementary)
      if (nargs==2 && hasArg(school.district.elementary))
       geo.obj=new(Class="geo", api.for=list("school+district+(elementary)"=school.district.elementary), api.in=list(state=state), name=paste(school.district.elementary.name, ", ", state.name, sep=""), sumlev=950)
      # sumlev 960 (state-school.secondary)
      if (nargs==2 && hasArg(school.district.secondary))
       geo.obj=new(Class="geo", api.for=list("school+district+(secondary)"=school.district.secondary), api.in=list(state=state), name=paste(school.district.secondary.name, ", ", state.name, sep=""), sumlev=960)
      # sumlev 970 (state-school.unified)
      if (nargs==2 && hasArg(school.district.unified))
       geo.obj=new(Class="geo", api.for=list("school+district+(unified)"=school.district.unified), api.in=list(state=state), name=paste(school.district.unified.name, ", ", state.name, sep=""), sumlev=970)
      # fail if still no match
      if (!is.geo(geo.obj)) {
        warning("No valid geography from these arguments; returning NA\n  Perhaps add or drop levels from this request?")
      }
      geo.obj
    }
    ## after here is actual function to recur through more complex
    ## requests
    # recycle!
    if (key=="auto" && check==T) load(system.file("extdata/key.rda", package="acs"))
    arglist=as.list(environment())
    missing.args=unlist(lapply(arglist, is.symbol))
    arglist=arglist[!missing.args]
    arglist=arglist[names(arglist)!="combine"]
    arglist=arglist[names(arglist)!="combine.term"]
    arglist=arglist[names(arglist)!="key"]
    arglist=arglist[names(arglist)!="check"]
    max.length=max(sapply(arglist, length))
    arglist=lapply(arglist, rep, length=max.length)
    if(max.length==1) {
      geo.obj=do.call(.geo.unit.make, arglist)
      geo.set.obj=new(Class="geo.set", geo.list=list(geo.obj), combine=combine, combine.term=combine.term)
    } else {
      geo.a=do.call(.geo.unit.make, lapply(arglist, head, 1))
      geo.b=do.call(geo.make, lapply(arglist, tail, -1))
      geo.set.obj=geo.a+geo.b
    }
    if(is.geo.set(geo.set.obj)) {
      geo.set.obj@combine=combine
      geo.set.obj@combine.term=combine.term}
    if (check==T){
      for (i in 1:length(geo.set.obj)){
        cat(paste("Testing geography item ", i,": ", name(geo.list(geo.set.obj[i])), " .... ", sep=""))
        obj.test=acs.fetch(endyear=2010, geography=geo.set.obj[i], key=key, variable="B01001_001",  lookup=F)
        cat("OK.\n")
      }
    }
    geo.set.obj
  }

setClass(Class="acs.lookup", representation =
      representation(endyear="numeric", span="numeric", args="list", results="data.frame"))
    if (!isGeneric("endyear")) {
      setGeneric("endyear", def=function(object){standardGeneric("endyear")})}else{}
    setMethod("endyear", "acs.lookup", function(object) object@endyear)
    if (!isGeneric("span")) {
      setGeneric("span", def=function(object){standardGeneric("span")})}else{}
    setMethod("span", "acs.lookup", function(object) object@span)
    if (!isGeneric("results")) {
      setGeneric("results", def=function(object){standardGeneric("results")})}else{}
    setMethod("results", "acs.lookup", function(object) object@results)
    # could add other slots, plus "show" method
    
    setMethod("show", "acs.lookup", function(object) {
      cat("An object of class \"acs.lookup\"\n")
      cat("endyear=", endyear(object)," ; span=", span(object), "\n\n")
      cat("results:\n")
      print(results(object))
      cat("\n")
      })
    
    
    is.acs.lookup=function(object){
      if (class(object)=="acs.lookup") {TRUE}
      else {FALSE}}
    
    setMethod("+", signature(e1 = "acs.lookup", e2 = "acs.lookup"), function(e1, e2) {
      e3=rbind(e1@results, e2@results)
      new(Class="acs.lookup", endyear=e1@endyear, args=list(e1@args, e2@args), span=e1@span, results=e3)})
    
    setMethod("c", signature(x = "acs.lookup" ), function(x, y, ...,
         recursive = FALSE) {
      if(missing(y)) x
      else x + c(y, ...)})
    
    setMethod(f="[", signature="acs.lookup", definition=function(x,i,j,...,drop=FALSE){
      if (missing(i)) i=j
      if (missing(j)) j=i
      new(Class="acs.lookup", endyear=x@endyear, args=x@args, span=x@span, results=x@results[i,])})
    
    acs.lookup=function(endyear, span=5, dataset="acs", keyword, table.name, table.number, case.sensitive=T) {
      arglist=as.list(environment())
      if (!missing(table.number)){
        if (!missing(table.name)) warning("Cannot specify both table.name and table.number; using table.number")
   # in future?: consider changing next line to table.name="", and let table.number drive the train
        if(endyear!=1990){table.name=paste(table.number, ".", sep="")} else {table.name=table.number}
}
      if (missing(table.name) && missing(keyword)) {
        warning("No search terms provided; returning NA")
        return(NA)}
      if (!missing(keyword) && sum(unlist(lapply(X=keyword, FUN=grepl, "Margin Of Error For", ignore.case=T)))>0) {
        warning("'keyword' marching string 'Margin Of Error For' not permitted\n  Returning NA")
        return(NA)}
      if (!case.sensitive) {
        if (!missing(table.name)) table.name=tolower(table.name)
        if (!missing(keyword)) keyword=tolower(keyword)}
      else {insensitive=F}
  # find correct XML variables
  # new way / updated for v2.0
  # doc.string is xml file name when saved locally or on eglenn archive
  # doc.url is path to census file for XML variables
      if(dataset=="acs") {
        doc.string=paste(dataset,"_", span,"yr_", endyear,"_var.xml.gz", sep="")
        doc.url=paste("http://api.census.gov/data/", endyear,"/acs",span,"/variables.xml", sep="")
      }
      if(dataset=="sf1" | dataset=="sf3"){
        doc.string=paste(dataset,"_", endyear,".xml.gz", sep="")
        doc.url=paste("http://api.census.gov/data/", endyear, "/", dataset, "/variables.xml", sep="")
        span=0
      }
      # first look for XML table internally
      if(file.exists(system.file(paste("extdata/", doc.string, sep=""), package="acs")))
      {
          doc=xmlInternalTreeParse(system.file(paste("extdata/", doc.string, sep=""), package="acs"))
      }
      # next check online at census site
      else if(url.exists(doc.url))
      {
          doc=xmlInternalTreeParse(doc.url)
      }
      # finally, check personal eglenn archive
      else if(url.exists(paste("http://web.mit.edu/eglenn/www/acs/acs-variables/", doc.string, sep="")))
      {
          # since only here is issues, give some advice
          warning(paste("XML variable lookup tables for this request\n  seem to be missing from '", doc.url, "';\n  temporarily downloading and using archived copies instead;\n  since this is *much* slower, recommend running\n  acs.tables.install()"), sep="")  
          doc.download=tempfile()
          download.file(url=paste("http://web.mit.edu/eglenn/www/acs/acs-variables/", doc.string, sep=""), destfile=doc.download)
          doc=xmlInternalTreeParse(doc.download)
          unlink(doc.download)
      }
      # if found nowhere, issue warning and return NA
      else {
         warning("As of the date of this version of the acs package\n  no variable lookup tables were available\n  for this dataset/endyear/span combination;\n  perhaps try a different combination...?\n  Returning NA;")
        return(NA)
      }
      
      if (!missing(keyword)){
          if (!case.sensitive) {str.a="contains(translate(@label, 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'),'"}
        else {str.a="contains(@label, '"}
          str.b=paste(str.a,keyword,"')", sep="")
          str.c=paste(str.b, collapse=" and ")
          keyword=str.c
      } else {keyword=""}
# add in stanza using table number
      if (!missing(table.name)) {
          if (!case.sensitive) {str.a="contains(translate(@concept, 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'),'"}
        else {str.a="contains(@concept, '"}
          str.b=paste(str.a,table.name,"')", sep="")
          str.c=paste(str.b, collapse=" and ")
          table.name=str.c
          if(endyear==1990){table.name=gsub(table.name, pattern="@concept", replacement="@xml:id")}
      } else {table.name=""}

      if(identical(table.name, "")){
             STRING=paste0("//ns:var[", keyword, "]")
      }
      else if(identical(keyword,"")){
             STRING=paste0("//ns:var[", table.name, "]")
      }
      else {
             STRING=paste0("//ns:var[", paste(table.name, keyword, sep=" and "), "]")
      }
# get variable codegs ("id")
      names=suppressWarnings(xpathSApply(doc, STRING, namespaces="ns", xmlGetAttr, "xml:id"))
      names=gsub("!!!!"," ",names)
      names=gsub("!!"," ",names)
      my.index=order(names)  # added for 2012, since data not sorted
      names=names[my.index]  # added for 2012, since data not sorted
      names=gsub("E$", "", names) # remove "E" from variable name
      if(dataset=="acs" && length(names)>0) {names=names[seq(1,length(names),2)]} # only want every other
    
# get table names 
      table.names=suppressWarnings(xpathSApply(doc, STRING, namespaces="ns", xmlGetAttr, "concept"))
      table.names=gsub("!!!!"," ",table.names)
      table.names=gsub("!!"," ",table.names)
      table.numbers=regmatches(table.names, m=regexpr(table.names, pattern="^.*\\.")) # find table numbers
      if(endyear==1990){
        table.numbers=substr(names, 1, 5)
        table.numbers=gsub(x=table.numbers, pattern="0$", replacement="")
        table.numbers=gsub(x=table.numbers, pattern="$", replacement=".")
      }
      table.names=gsub(x=table.names, pattern="^.*\\.  ", replacement="") # remove table numbers from names 
      table.names=gsub(x=table.names, pattern="* \\[.*\\]$", replacement="") # remove bracketed variable counts from SF1/SF3 tables 
      table.names=table.names[my.index] # added for 2012, since data not sorted
      if(dataset=="acs" && length(table.names)>0) {table.names=table.names[seq(1,length(table.names),2)]} # only want every other
      
# get table numbers
      table.numbers=substr(table.numbers, 1, unlist(lapply(table.numbers, nchar))-1) # remove trailing period
      if(dataset=="acs"){  # sf1/sf3 table.numbers have already been reordered!
        table.numbers=table.numbers[my.index]
      } # added for 2012, since data not sorted
      if(dataset=="acs" && length(table.numbers)>0) {table.numbers=table.numbers[seq(1,length(table.numbers),2)]} # only want every other
    
# get variable names
      values=suppressWarnings(xpathSApply(doc, STRING, namespaces="ns", xmlGetAttr, "label"))
      values=gsub("!!!!"," ",values)
      values=gsub("!!"," ",values)
      values=gsub(x=values, pattern="*\\[.*\\]", replacement=":") # remove bracketed variable counts from SF1/SF3 tables 
      values=values[my.index] # added for 2012, since data not sorted
      if(dataset=="acs"  && length(values)>0) {values=values[seq(1,length(values),2)]} # only want every other
    
      if (length(names)==0){
        warning("Sorry, no tables/keyword meets your search.\n  Suggestions:\n    try with 'case.sensitive=F',\n    remove search terms,\n    change 'keyword' to 'table.name' in search (or vice-versa)")
        return(NA)}
      free(doc)
      rm(doc)
      gc()
      new(Class="acs.lookup", endyear=endyear, span=span, args=arglist, 
            results=data.frame(variable.code=names, 
            table.number=table.numbers, table.name=table.names, variable.name=values, 
            stringsAsFactors=F))
    }

read.acs=function(filename, endyear="auto", span="auto",
  col.names="auto", acs.units="auto", geocols="auto", skip="auto")
  {
                                        # Attempt to automatically determine endyear from filename if
                                        # necessary.
    
    if (endyear=="auto")
      { # try to guess end year from filename
        endyear=2000+as.integer(str_extract(str_extract(filename, "ACS_[0-9][0-9]_"), "[0-9][0-9]"))
        if(is.na(endyear)) endyear=as.integer(str_extract(filename, "200[0-9]"))
        if (is.na(endyear) | endyear<2000 | endyear>2012)
          {
            warning("Can't determine endyear from filename;\nplease set manually before proceeding.\nSetting endyear to default of 2010.\nOperations with this acs object may not be reliable...")
            endyear=2010}
        else {
          warning("Guessing endyear of " , endyear, " based on filename...", call.=F)
        }
      }
                                        # Attempt to automatically determine span from filename if
                                        # necessary.
  
    if (span=="auto")
      {
        span=as.integer(str_extract(str_extract(filename, "_[0-9]YR"), "[0-9]"))
        if (is.na(span)) span=as.integer(substr(str_extract(filename, "[0-9]yr"),1,1))
        if (is.na(span) | span>5) {
          warning("Can't determine span from filename;\nplease set manually before proceeding.\nSetting span to default of 1.\nOperations with this acs object may not be reliable...")
          span=1
        }
        else
          {
            warning("Guessing span of ", span, " based on filename...", call.=F)
          }
      }
    span=as.integer(span)
    endyear=as.integer(endyear) 
    if (span > 5 | span < 1){
      warning("Span out of range; returning NA")
      return(NA)}
                                        # set geocols
    if (identical(geocols, "auto"))
      {geocols=3:1
       warning("Using first three columns as geographic headers.", call.=F)
     }
    
    if (identical(str_sub(filename, start=-4), ".zip")) {
      zip.filename=filename
      contents=unzip(zip.filename, list=T)
      acs.filename=grep(pattern="with_ann.csv", x=contents[[1]], value=T)
      if(length(acs.filename)==0)
        acs.filename=grep(pattern="[0-9].csv", x=contents[[1]], value=T)
      make.con=function(){unz(description=zip.filename, filename=acs.filename)
                          }
    } else {make.con=function(){file(description=filename)
                                }
          }
                                        # figure out how many rows to skip -- only call if skip="auto"
    figure.skip=function(){
      con=make.con()
      open(con)
      i=0
      while(str_sub(scan(con, what="", nlines=1, quiet=T), end=1)[1]==","){
        i=i+1
      }
      close(con)
      i+1  # need to add one more to skip extra column header; most ACS files have two rows for headers
    }
    if (identical(skip, "auto")){
      skip=figure.skip()
    }
    
                                        # figure out column names based on headers
    get.colnames=function(geocols, skip){
      con=make.con()
      open(con)
      headers=read.csv(con, nrows=skip+1, header=F)  # add one to include "header" row
      headers=headers[,-geocols]
      headers=apply(headers, FUN="paste", MARGIN=2, collapse=".")
      headers=headers[seq(1,length(headers),2)]
      close(con)
      headers
    }
    if (identical(col.names, "auto")) {
      col.names=get.colnames(geocols, skip)
    }
    # helper function to read and clean data        
    get.acs=function(geocols, skip) {
      con=make.con()
      open(con)
      in.data=read.csv(con, skip=skip, na.strings=c("-", "**",
      "***", "(X)", "N"), stringsAsFactors=F)
      # trick to get geocols to be first columns
      datacols=1:length(in.data)
      datacols=datacols[-geocols]
      in.data=in.data[,c(geocols,datacols)]
      in.data[in.data=="*****"]=0
      for (i in (length(geocols)+1):length(in.data)) {
        in.data[[i]]=gsub(",", "", in.data[[i]])
        in.data[[i]]=as.numeric(in.data[[i]])  
      }
      colnames(in.data)[-geocols]=col.names
      close(con)
      in.data
    }
      
    in.data=get.acs(geocols=geocols, skip=skip)
    acs.colnames=unname(col.names)
#  acs.colnames=sub(colnames(in.data[seq((length(geocols)+1),length(in.data),2)]), pattern="..Estimate.", replacement="")
                                          # create new acs object
    if (acs.units=="auto") {
                                          # try to guess variable types from filename
      acs.units=.acs.identify.units(acs.colnames)
    }
    acs.units=factor(acs.units, levels=.acs.unit.levels)
    acs.obj=new(Class="acs",
      endyear=endyear,
      span=span, 
      geography=as.data.frame(in.data[,1:length(geocols)]),
      acs.colnames=acs.colnames, 
      acs.units=acs.units,
      currency.year=endyear,
      standard.error=as.matrix(in.data[,seq((length(geocols)+2),length(in.data), 2)]),
      modified=F,
      estimate=as.matrix(in.data[,seq((length(geocols)+1),length(in.data), 2)]
        )
      )
    
  # convert 90% MOE into standard error, correct for 2005 flaw
    if (endyear(acs.obj)<=2005)
      acs.obj@standard.error=acs.obj@standard.error/1.65
    else
      acs.obj@standard.error=acs.obj@standard.error/1.645
    acs.obj=.acs.dimnames(acs.obj)
    acs.obj
  }

acs.fetch=function(endyear, span=5, geography, table.name,
      table.number, variable, keyword, dataset="acs", key, col.names="auto", ...)
      {
        var.max=40  # most var for acs call; keep even; 
        # check some basic stuff about arguments
        if (missing(key)){
          if (file_test("-f", system.file("extdata/key.rda", package="acs"))) {
            load(system.file("extdata/key.rda", package="acs"))
          } else {
            warning("'key' required to access Census API site for download;\n  See http://www.census.gov/developers/ to request a key\n  and/or use 'key=' (or run 'api.key.install()') to avoid this error.")
            return(NA) }
        }
         if (missing(endyear)) {
             warning("No endyear provided\n  As of version 2.0, endyear must be explicit.\n  Returning NA.")
             return(NA)
        }
        endyear=as.integer(endyear)
        span=as.integer(span)
         if (span < 0 | span > 5)
           {
             warning("span must be 1, 3, or 5 (or 0, for decennial census)")
           }
        if(dataset=="sf1" | dataset=="sf3") {span=as.integer(0)}
        if (missing(keyword) && missing(table.name) && missing(table.number) && missing(variable)) {
          warning("No search terms provided; returning NA")
          return(NA)}
        if (!missing(variable)){
          if (!missing(keyword) || !missing(table.name) || !missing(table.number)) {
            warning("Cannot specify both 'variable' and 'keyword/table.name/table.number'.\n  Using 'variable' and ignoring others.")}
        }
        # when variable is NOT provided
        if (missing(variable)) {
          arglist=as.list(environment())
          missing.args=unlist(lapply(arglist, is.symbol))
          arglist=arglist[!missing.args]
#          arglist$dataset="acs"
          arglist=arglist[names(arglist)!="variable"]
          arglist=arglist[names(arglist)!="geography"]
          arglist=arglist[names(arglist)!="key"]
          arglist=arglist[names(arglist)!="col.names"]
          arglist=arglist[names(arglist)!="var.max"]
  #        return(arglist)
          variable=do.call(acs.lookup, arglist)
          if(!isS4(variable) && is.na(variable)){return(NA)}}
        # when variable is provided
        if (is.acs.lookup(variable)) {
          variables.xml=results(variable)
          variables=variables.xml$variable.code
# next "if" added to allow for sf1/sf3 fetching
          if(dataset=="acs") {variables=paste(rep(variables, each=2), c("E","M"), sep="")}
        }
        if (!is.acs.lookup(variable)) {
          if (variable[1]=="recur") {
            variables=variable[2:length(variable)]}
          else {
            if(dataset=="acs") {variables=paste(rep(variable, each=2), c("E","M"), sep="")}
            if(dataset=="sf1" | dataset=="sf3"){variables=variable}
            }
          if (variable[1]=="") {
            variables.xml=acs.lookup(keyword=keyword, table.name=table.name, endyear=endyear, dataset=dataset, ...)
            if (!identical(NA, variables.xml)){
              variables.xml=results(variables.xml)
              variables=variables.xml$variable.code
            }
            else {warning("No results found;\n  perhaps try acs.lookup()...?")
                  return(NA)}
          }
        }
        if (length(variables)==1 && is.na(variable)) return(NA)
                                            # pretty option to pull descriptive names from XML lookup results; only take every other, since these are the headers for estimates, not MOEs
        if (identical(col.names, "pretty")) {
          if (!is.acs.lookup(variable)) {
            warning("\"pretty\" col.names not available when variable codes provided.\n  Using standard variable code names for columns.")
            col.names="auto"}
          else {
            col.names=paste(variables.xml$table.name, variables.xml$variable.name, sep=": ")
          }}
    #    if (identical(acs.units, "auto")) acs.units=rep("auto",(length(variables)/2))
        if (identical(col.names, "auto") & dataset=="acs") col.names=rep("auto",(length(variables)/2))
        if (identical(col.names, "auto") & (dataset=="sf1" | dataset=="sf3")) col.names=rep("auto",length(variables))
        # deal with too many variables -- API currently limits to < 50
        # note two versions: acs datasets have doubled variables
        if (length(variables) > var.max & dataset=="acs") {
          acs.obj=cbind(acs.fetch(endyear=endyear, span=span, geography=geography, variable=c("recur", variables[1:(var.max-2)]), key=key, dataset=dataset, col.names=col.names[1:((var.max-2)/2)]), acs.fetch(endyear=endyear, span=span, geography=geography, variable=c("recur", variables[(var.max-1):length(variables)]), col.names=col.names[(1+((var.max-2)/2)):length(col.names)], key=key, dataset=dataset))
          return(acs.obj)
        }
        # note two versions:sf1/sf3 datasets have single variables, since no E and M types
         if (length(variables) > var.max & (dataset=="sf1" | dataset=="sf3")) {
          acs.obj=cbind(acs.fetch(endyear=endyear, span=span, geography=geography, variable=c("recur", variables[1:(var.max-2)]), key=key, dataset=dataset, col.names=col.names[1:(var.max-2)]), acs.fetch(endyear=endyear, span=span, geography=geography, variable=c("recur", variables[(var.max-1):length(variables)]), col.names=col.names[(var.max-1):length(col.names)], key=key, dataset=dataset))
          return(acs.obj)
        }

        # next deal with complex geographies
        if (is.geo.set(geography) && length(geography)>1){
          acs.obj=rbind(acs.fetch(endyear=endyear, span=span, geography=geography[1], variable=c("recur", variables), key=key, dataset=dataset, col.names=col.names), acs.fetch(endyear=endyear, span=span, geography=geography[2:length(geography)], variable=c("recur", variables), dataset=dataset, key=key, col.names=col.names))
          if(combine(geography)){
            acs.obj=apply(acs.obj, FUN=sum, MARGIN=1,
          agg.term=combine.term(geography), ...)
            acs.obj@acs.units=.acs.identify.units(acs.colnames(acs.obj))
          }
          return(acs.obj)
        }
        if (is.geo.set(geography) && length(geography)==1) {
          acs.obj=acs.fetch(endyear=endyear, span=span, geography=geography[[1]], variable=c("recur", variables), key=key, dataset=dataset, col.names=col.names)
          if (combine(geography)) {
            acs.obj=apply(acs.obj, FUN=sum, MARGIN=1,
          agg.term=combine.term(geography), ...)
            acs.obj@acs.units=.acs.identify.units(acs.colnames(acs.obj))
          }
          return(acs.obj)
        }
    # after this point, geography should be always be a single geo
        api.url=api.url.maker(endyear=endyear, span=span, key=key, variables=variables, dataset=dataset, geo.call=geography)
        geo.length=length(api.in(geography))+2
   # adding check to stop bad url / maybe do this later
        url.test=url.exists(api.url, .header=T)
        if(url.test["statusMessage"]!="OK") {
            warning(call.=F, paste("No data found at:\n  ", api.url, sep=""))
        }
        in.data=suppressWarnings(read.csv(api.url, na.strings=c("-", "**", "***", "(X)", "N", "null"), stringsAsFactors=F))
        in.data=in.data[,-length(in.data)] # remove junk NA columns
                                            # set geocols
        geocols=(length(in.data)-geo.length+1):length(in.data)
        ## col.names=names(in.data)[seq(1,(length(in.data)-geo.length), 2)]
        ## col.names[1]=gsub("X..","",col.names[1])
        if (identical(col.names[1], "auto")) {  # check this!
          if(dataset=="acs"){
          col.names=names(in.data)[seq(1,(length(in.data)-geo.length), 2)]
          } else if(dataset=="sf1" | dataset=="sf3") {
          col.names=names(in.data)[1:(length(in.data)-geo.length)]
          }
          col.names[1]=gsub("X..","",col.names[1])
          col.names=gsub(pattern="E$", x=col.names, replacement="")
          
        }
        datacols=1:(length(in.data)-geo.length)
        in.data[in.data=="*****"]=0
        in.data[[1]]=gsub("[","",in.data[[1]], fixed=T )
        in.data[[length(in.data)]]=gsub("]","",in.data[[length(in.data)]], fixed=T )
                                            # clean brackets
        for (i in 1:length(datacols)) {
          in.data[[i]]=gsub(",", "", in.data[[i]])
          in.data[[i]]=as.numeric(in.data[[i]])  
        }
                                            # create new acs object
    #    if (identical(acs.units[1], "auto")) {
                                              # try to guess variable types from filename
    #    acs.units=.acs.identify.units(col.names)
    #    }
    #    acs.units=factor(acs.units, levels=.acs.unit.levels)
        GEOGRAPHY=as.data.frame(in.data[,geocols])
        names(GEOGRAPHY)=gsub(".", "", names(GEOGRAPHY), fixed=T) # remove strange trailing period
        if(dataset=="acs"){
        acs.obj=new(Class="acs",
          endyear=endyear,
          span=span, 
          geography=GEOGRAPHY,
          acs.colnames=col.names,
          acs.units=.acs.identify.units(col.names),
          currency.year=endyear,
          standard.error=as.matrix(in.data[,seq(2,(length(in.data)-geo.length), 2)]),
          modified=F,
          estimate=as.matrix(in.data[,seq(1,(length(in.data)-geo.length), 2)])
          )}
        if(dataset=="sf1" | dataset=="sf3"){
        acs.obj=new(Class="acs",
          endyear=endyear,
          span=span, 
          geography=GEOGRAPHY,
          acs.colnames=col.names,
          acs.units=.acs.identify.units(col.names),
          currency.year=endyear,
          standard.error=as.matrix(0*(in.data[1:(length(in.data)-geo.length)])),
          modified=F,
          estimate=as.matrix(in.data[1:(length(in.data)-geo.length)])
          )
         }     
      # convert 90% MOE into standard error, correct for 2005 flaw
        if (endyear(acs.obj)<=2005)
          acs.obj@standard.error=acs.obj@standard.error/1.65
        else
          acs.obj@standard.error=acs.obj@standard.error/1.645
        acs.obj=.acs.dimnames(acs.obj)
    #    if (obj.combine) {acs.obj=apply(acs.obj, FUN=sum, MARGIN=1)}
        acs.obj
      }

if (!isGeneric("endyear")) {
  setGeneric("endyear", def=function(object){standardGeneric("endyear")})}else{}
setMethod("endyear", "acs", function(object) object@endyear)

if (!isGeneric("span")) {
  setGeneric("span", def=function(object){standardGeneric("span")})}else{}
setMethod("span", "acs", function(object) object@span)

if (!isGeneric("geography")) {
  setGeneric("geography", def=function(object){standardGeneric("geography")})}else{}
setMethod("geography", "acs", function(object) object@geography)

if (!isGeneric("acs.colnames")) {
  setGeneric("acs.colnames", def=function(object){standardGeneric("acs.colnames")})}else{}
setMethod("acs.colnames", "acs", function(object) object@acs.colnames)

if (!isGeneric("currency.year")) {
  setGeneric("currency.year", def=function(object){standardGeneric("currency.year")})}else{}
setMethod("currency.year", "acs", function(object) object@currency.year)

if (!isGeneric("modified")) {
  setGeneric("modified", def=function(object){standardGeneric("modified")})}else{}
setMethod("modified", "acs", function(object) object@modified)

if (!isGeneric("acs.units")) {
  setGeneric("acs.units", def=function(object){standardGeneric("acs.units")})}else{}
setMethod("acs.units", "acs", function(object) object@acs.units)

if (!isGeneric("estimate")) {
  setGeneric("estimate", def=function(object, which, conf.lev, ...){standardGeneric("estimate")})}else{}
setMethod("estimate", "acs", function(object) object@estimate)

if (!isGeneric("standard.error")) {
  setGeneric("standard.error", def=function(object){standardGeneric("standard.error")})}else{}
setMethod("standard.error", "acs", function(object) object@standard.error)

setMethod(f="[", signature="acs", definition=function(x,i,j,...,drop=FALSE){
  if (missing(i)) i=1:dim(x@estimate)[1]
  if (missing(j)) j=1:dim(x@estimate)[2]
  new(Class="acs",
      endyear=endyear(x),
      span=span(x),
      geography=geography(x)[i,],
      acs.colnames=acs.colnames(x)[j],
      modified=modified(x),
      acs.units=acs.units(x)[j],
      currency.year=currency.year(x),
      estimate=estimate(x)[i,j, drop=F],
      standard.error=standard.error(x)[i,j, drop=F])
})

setReplaceMethod(f="[", signature="acs",
    definition=function(x,i,j,value){
      if (missing(i)) i=1:dim(x)[1]
      if (missing(j)) j=1:dim(x)[2]
# is value acs object? ## still need to check for metadata being the same
      if (is.acs(value) && all(dim(value)==c(length(i),length(j)))){
        if (endyear(x) != endyear(value)) {
          warning("original and replacement do not have same endyear;\nkeeping original value", call.=F)}
        if (span(x) != span(value)) {
          warning("original and replacement do not have same span;\nkeeping original value", call.=F)}
        if (currency.year(x) != currency.year(value)) {
          warning("original and replacement do not have same currency.year;\nkeeping original value", call.=F)}
        x@estimate[i,j]=value@estimate
        x@standard.error[i,j]=value@standard.error
# check for mismatch geo when not all cols changed
        if (!all(geography(x[i,j])==geography(value))){  # if not identical geogs
          if (dim(x)[2]<=length(j)){      # if changing all cols or more
            x@geography[i,]=geography(value)           # change all geo
            warning("geographies do not match but all columns changed;\nusing new geographies", call.=F)
          }else{
            warning("geographies do not match but some columns retained;\nkeeping original geography values", call.=F)
        }
        }
        if (!all(acs.colnames(x[i,j])==acs.colnames(value))){   # if not identical colnames
          if (dim(x)[1]<=length(i)){             # if not changing all rows or more
            x@acs.colnames[j]=acs.colnames(value)
            warning("acs.colnames do not match but all rows changes;\nusing new acs.colnames", call.=F)
          }else{
            warning("acs.colnames do not match but some rows retained;\nkeeping original acs.colnames", call.=F)
          }
        }
                                            # is value two item list?
      } else if (is.list(value) && length(value)==2) {
          if (is.null(value$estimate)) x@estimate[i,j]=value[[1]]
          else x@estimate[i,j]=value$estimate
          if (is.null(value$standard.error)){
            if (is.null(value$error)) x@standard.error[i,j]=value[[2]]
            else x@standard.error[i,j]=value$error
          } else x@standard.error[i,j]=value$standard.error
          ## standard.error --> 0
                                      # is value a single number?
          ## } else if (is.numeric(value) && (length(value)==1)) {
          ##   x@estimate[i,j]=value
          ##   x@standard.error[i,j]=0
          ##                                 # is value a vector of numbers?
          ## } else if (is.numeric(value) && (length(value)==length(x[i,j]))){
          ##   x@estimate[i,j]=value
          ##   x@standard.error[i,j]=0
          ## next stanza does the work of both of the previous:
        } else if (is.numeric(value)) {
          x@estimate[i,j]=value
          x@standard.error[i,j]=0
        } else {stop("incompatible objects or dimensions;\nunable to parse for replacement", call.=F)}
      x@modified=T
      x=.acs.dimnames(x)   # in case geography or acs.colnames changed
      validObject(x)
      return(x)
    })

cbind.acs=function(e1, e2) {
  if (e1@endyear != e2@endyear | e1@span != e2@span) {
    warning("** acs objects x and y must have same endyear and span;\nreturning NA **")
    return(NA)}
  if (identical(geography(e1), geography(e2))) GEOGRAPHY=geography(e1)
  else {warning( "geographies do not appear to match; using first geography")
        GEOGRAPHY=geography(e1)}
  NEW.ESTIMATE=cbind(estimate(e1), estimate(e2))
  NEW.ERROR=cbind(standard.error(e1), standard.error(e2))
  acs.obj=new(Class="acs", endyear=endyear(e1), span=span(e1), modified=T, geography=GEOGRAPHY, acs.units=factor(c(acs.units(e1),acs.units(e2)), levels=.acs.unit.levels), currency.year=currency.year(e1), acs.colnames=c(acs.colnames(e1), acs.colnames(e2)), estimate=NEW.ESTIMATE, standard.error=NEW.ERROR)
  acs.obj=.acs.dimnames(acs.obj)
  acs.obj
}

rbind.acs=function(e1, e2) {
  if (e1@endyear != e2@endyear | e1@span != e2@span) {
    warning("** acs objects x and y must have same endyear and span;\nreturning NA **")
    return(NA)}
  if (identical(acs.colnames(e1), acs.colnames(e2))) ACS.COLNAMES=acs.colnames(e1)
  else {warning( "columns do not appear to match; using first colnames")
        ACS.COLNAMES=acs.colnames(e1)}
  GEOGRAPHY=rbind.fill(geography(e1), geography(e2))
  NEW.ESTIMATE=rbind(estimate(e1), estimate(e2))
  NEW.ERROR=rbind(standard.error(e1), standard.error(e2))
  acs.obj=new(Class="acs", endyear=endyear(e1), span=span(e1), modified=T, geography=GEOGRAPHY, acs.units=acs.units(e1), currency.year=currency.year(e1), acs.colnames=acs.colnames(e1), estimate=NEW.ESTIMATE, standard.error=NEW.ERROR)
  acs.obj=.acs.dimnames(acs.obj)
  acs.obj
}

setMethod("+", signature(e1 = "acs", e2 = "acs"), function(e1, e2) {
    header=.acs.combine.headers(e1,e2,"+")
    NEW.ESTIMATE=estimate(e1)+estimate(e2)
    NEW.ERROR=sqrt(standard.error(e1)^2+standard.error(e2)^2)
    acs.obj=new(Class="acs",
      endyear=header$endyear,
      span=header$span,
      modified=T,
      geography=header$geography,
      acs.units=header$acs.units,
      currency.year=header$currency.year,
      acs.colnames=header$acs.colnames,
      estimate=NEW.ESTIMATE,
      standard.error=NEW.ERROR)
    acs.obj=.acs.dimnames(acs.obj)
    acs.obj
  }
           )
  
  setMethod("-", signature(e1 = "acs", e2 = "acs"), function(e1, e2) {
    header=.acs.combine.headers(e1, e2, "-")
    NEW.ESTIMATE=estimate(e1)-estimate(e2)
    NEW.ERROR=sqrt(standard.error(e1)^2+standard.error(e2)^2)
    acs.obj=new(Class="acs",
      endyear=header$endyear,
      span=header$span,
      modified=T,
      geography=header$geography,
      acs.units=header$acs.units,
      currency.year=header$currency.year,
      acs.colnames=header$acs.colnames,
      estimate=NEW.ESTIMATE,
      standard.error=NEW.ERROR)
    acs.obj=.acs.dimnames(acs.obj)
    acs.obj
  }
           )

 ## old; works, but not great.  
  # .acs.divider=function(num, den, proportion, verbose=F) {
  #   if (proportion==T) header=.acs.combine.headers(num, den, "/")
  #   else header=.acs.combine.headers(num, den, ":")
  #   p=estimate(num)/estimate(den)
  #   if (all(estimate(den)!=0) && proportion==T & all((p^2 * standard.error(den)^2)>0)){
  #     header$acs.units=factor("proportion", levels=.acs.unit.levels)
  #     if (verbose) {warning("** using formula for PROPORTIONS, which assumes that numerator is a SUBSET of denominator **")}
  #     NEW.ERROR=sqrt(standard.error(num)^2 - (p^2 * standard.error(den)^2))/estimate(den)}
  #   else {
  # #   a recommended correction when term under sqrt is negative
  #     if (proportion==T){
  #       if(verbose){warning("** due to the nature of some of the errors, using the more conservative formula for RATIOS, which assumes that numerator is not a subset of denominator **")}}
  #     header$acs.units=factor("ratio", levels=.acs.unit.levels)
  #     NEW.ERROR=sqrt(standard.error(num)^2 + (p^2 * standard.error(den)^2))/estimate(den)
  #   }
  #   acs.obj=new(Class="acs",
  #     endyear=header$endyear,
  #     span=header$span,
  #     modified=T,
  #     geography=header$geography,
  #     acs.units=header$acs.units,
  #     currency.year=header$currency.year,
  #     acs.colnames=header$acs.colnames,
  #     estimate=p,
  #     standard.error=NEW.ERROR)
  #   acs.obj=.acs.dimnames(acs.obj)
  #   acs.obj}

# new, to deal with zeroes, and more precise ratio-style correction

  .acs.divider=function(num, den, proportion, verbose=F, output="result") {
    if (proportion==T) header=.acs.combine.headers(num, den, "/")
    else header=.acs.combine.headers(num, den, ":")
    p=estimate(num)/estimate(den)
# start with proportion-style
    if (proportion==T){
        header$acs.units=rep(factor("proportion", levels=.acs.unit.levels), length(header$acs.units))
        if (verbose) {warning("** using formula for PROPORTIONS, which assumes that numerator is a SUBSET of denominator **")}  
        NEW.ERROR=suppressWarnings(sqrt(standard.error(num)^2 - (p^2 * standard.error(den)^2))/estimate(den))
        # change all that are should be ratio-stye
        # index for ratio corrections
        proportion.numerators=suppressWarnings(sqrt(standard.error(num)^2 - (p^2 * standard.error(den)^2)))
        ratio.correct.index=is.na(proportion.numerators)
      # if any fail the test
      if(any(ratio.correct.index)){
          # use ratio-style correction
          ratio.errors=sqrt(standard.error(num)^2 + (p^2 * standard.error(den)^2))/estimate(den)
          NEW.ERROR[ratio.correct.index]=ratio.errors[ratio.correct.index]
          if(verbose) {
              warning(paste("**Note: due to the nature of the errors in some cells,\n  they were divided using the more conservative formula for RATIOS\n  which assumes that numerator is not a subset of denominator**:\n  in total, ", sum(ratio.correct.index), " ratio-style divisions substituted;\n  see ?divide.acs and/or use output=\"result\" for more.", sep=""))}
          if(output=="both" | output=="div.method"){
              ratio.report=p  # just to get a matrix with the same dims
              colnames(ratio.report)=header$acs.colnames
              rownames(ratio.report)=header$geography[[1]]
              ratio.report[ratio.correct.index]="ratio"
              ratio.report[!ratio.correct.index]="proportion"
          }}}
    else {
  #   ratio style
        header$acs.units=rep(factor("ratio", levels=.acs.unit.levels), length(header$acs.units))
        NEW.ERROR=sqrt(standard.error(num)^2 + (p^2 * standard.error(den)^2))/estimate(den)
        ratio.report=p  # just to get a matrix with the same dims
        colnames(ratio.report)=header$acs.colnames
        rownames(ratio.report)=header$geography[[1]]
        ratio.report[,]="ratio"
    }
    acs.obj=new(Class="acs",
      endyear=header$endyear,
      span=header$span,
      modified=T,
      geography=header$geography,
      acs.units=header$acs.units,
      currency.year=header$currency.year,
      acs.colnames=header$acs.colnames,
      estimate=p,
      standard.error=NEW.ERROR)
    acs.obj=.acs.dimnames(acs.obj)
    if(output=="both") {list("result"=acs.obj, "div.method"=ratio.report)}
    else if(output=="div.method"){ratio.report}
    else acs.obj
  }
      

      
  setMethod("/", signature(e1 = "acs", e2 = "acs"), function(e1, e2) {
  # by default, use more conservative "ratio-style" dividing
    warning("** using the more conservative formula for ratio-type
    dividing, which does not assume that numerator is a subset of
    denominator; for more precise results when seeking a proportion
    and not a ratio, use divide.acs(..., method=\"proportion\") **")
    .acs.divider(num=e1, den=e2, proportion=F, verbose=F)
  })
  
  divide.acs=function(numerator, denominator, method="ratio", verbose=T, output="result"){
    if(identical(method, "ratio")) {
      .acs.divider(num=numerator, den=denominator, proportion=F, verbose=verbose, output=output)
    } else if(identical(method, "proportion")) {
      .acs.divider(num=numerator, den=denominator, proportion=T, verbose=verbose, output=output)
    } else {warning("Error: must set method to \"ratio\" or \"proportion\"")
            NA}
  }
  
  setMethod("*", signature(e1 = "acs", e2 = "acs"), function(e1, e2) {
    header=.acs.combine.headers(e1, e2, "*")
    NEW.ESTIMATE=estimate(e1)*estimate(e2)
    NEW.ERROR=sqrt((estimate(e1)^2*standard.error(e2)^2)+(estimate(e2)^2*standard.error(e1)^2))
    acs.obj=new(Class="acs",
      endyear=header$endyear,
      span=header$span,
      modified=T,
      geography=header$geography,
      acs.units=header$acs.units,
      currency.year=header$currency.year,
      acs.colnames=header$acs.colnames,
      estimate=NEW.ESTIMATE,
      standard.error=NEW.ERROR)
    acs.obj=.acs.dimnames(acs.obj)
    acs.obj}
  )

setMethod("+", signature(e1 = "acs", e2 = "numeric"), function(e1, e2) {
  e2=.acs.make.constant.object(value=e2, template=e1)
  e1+e2
  }
         )

# and reverse classes...

setMethod("+", signature(e1 = "numeric", e2 = "acs"), function(e1, e2) {
  e1=.acs.make.constant.object(value=e1, template=e2)
  e1+e2
  }
         )

setMethod("-", signature(e1 = "acs", e2 = "numeric"), function(e1, e2) {
  e2=.acs.make.constant.object(value=e2, template=e1)
  e1-e2
  }
         )

# ditto...

setMethod("-", signature(e1 = "numeric", e2 = "acs"), function(e1, e2) {
  e1=.acs.make.constant.object(value=e1, template=e2)
  e1-e2
  }
         )

setMethod("/", signature(e1 = "acs", e2 = "numeric"), function(e1, e2) {
  e2=.acs.make.constant.object(value=e2, template=e1)
  e1/e2
  }
         )

# ditto...

setMethod("/", signature(e1 = "numeric", e2 = "acs"), function(e1, e2) {
  e1=.acs.make.constant.object(value=e1, template=e2)
  e1/e2
  }
         )

setMethod("*", signature(e1 = "acs", e2 = "numeric"), function(e1, e2) {
  e2=.acs.make.constant.object(value=e2, template=e1)
  e1*e2
  }
         )

# ditto...

setMethod("*", signature(e1 = "numeric", e2 = "acs"), function(e1, e2) {
  e1=.acs.make.constant.object(value=e1, template=e2)
  e1*e2
  }
         )

setMethod("show", signature(object = "acs"), function(object) {
  if(is.na(span(object)) | is.na(endyear(object))) years="NO MEANINGFUL YEAR"
  else
    if(span(object)==0 | span(object)==1) years=endyear(object)
    else years=paste(endyear(object)-span(object)+1,"--",endyear(object))
  if(span(object)==0) {dataset="Decennial Census (SF1/SF3)" } else {dataset="ACS"}
  cat(dataset, "DATA: \n", years, ";\n")
  cat("  Estimates w/90% confidence intervals;\n  for different intervals, see confint()\n")
  est=estimate(object)
  err=standard.error(object)
  output=matrix(paste(est, "+/-", 1.645*err), nrow=nrow(est))
  dimnames(output)=dimnames(est)
  print(output, quote=FALSE)
})

setMethod("summary", signature(object = "acs"), function(object) {
  if(span(object)==1) years=endyear(object)
  else years=paste(endyear(object)-span(object)+1,"--",endyear(object))
  cat("ACS DATA: \n", years,"\n")
  cat("----------\nESTIMATES:\n")
  print(summary(estimate(object)))
  cat("----------\n90% MARGINS OF ERROR:\n")
  print(summary(1.645*standard.error(object)))
})

dim.acs=function(x){
  dim(estimate(x))}

length.acs=function(x){
  length(estimate(x))}

confint.acs=function(object, parm="all", level=.95, alternative="two.sided",...) {
  if (parm[1]=="all") parm=1:dim(object)[2]
  z.upper=switch(alternative, two.sided=qnorm((1+level)/2), greater=Inf, less=qnorm(level))
  z.lower=switch(alternative, two.sided=qnorm((1+level)/2), greater=qnorm(level), less=Inf)
  labels=switch(alternative, two.sided=c((1-level)/2 , 1 - (1-level)/2), less=c(0,level), greater=c(1-level,1))
  labels=paste(100*labels, "%", sep=" ")
  RESULT=list()
  for (i in parm) {
    conf.int.lower=estimate(object[,i])-standard.error(object[,i])*z.lower
    conf.int.upper=estimate(object[,i])+standard.error(object[,i])*z.upper
    RESULT[[acs.colnames(object)[i]]]=data.frame(conf.int.lower, 
            conf.int.upper,
            row.names=geography(object)[[1]])
    names(RESULT[[acs.colnames(object)[i]]])=labels
  }
  RESULT
}

setMethod("sum", signature(x = "acs"), function(x,
  agg.term=c("aggregate", "aggregate"), one.zero=FALSE, ..., na.rm=FALSE) {
  if(length(agg.term)<2){agg.term[2]=agg.term[1]}
  est=estimate(x)
  err=standard.error(x)
  if (one.zero==T && any(est==0)) {
    max.zero.error=max(err[est==0])
    err[est==0]=c(max.zero.error,rep(0,sum(est==0)-1))
  }
  if (dim(est)[1]==1){ geography=geography(x)} # single row
  else {
    geography=geography(x[1,1])
    for (i in 1:length(geography)){
      geography[,i]=agg.term[1]}
  }
  if (dim(est)[2]==1){acs.units=acs.units(x) # single column
                      acs.colnames=acs.colnames(x)
                    }
  else {acs.units=factor(levels=.acs.unit.levels)
                      acs.colnames=agg.term[2]}
  ESTIMATE=as.matrix(sum(est))
  ERROR=as.matrix(sqrt(sum(err^2)))
  acs.obj=new(Class="acs",
    endyear=endyear(x),
    span=span(x),
    modified=T,
    geography=geography,
    acs.units=acs.units,
    currency.year=currency.year(x),
    acs.colnames=acs.colnames,
    estimate=ESTIMATE,
    standard.error=ERROR)
  acs.obj=.acs.dimnames(acs.obj)
  acs.obj
})

.apply.acs=function(X, MARGIN, FUN, ...) 
   {
     FUN=match.fun(FUN)
     if (identical(MARGIN, 1)){       # apply row-wise
        acs.obj=FUN(X[,1], ...)
        if(dim(X)[2]>1){
       for (i in 2:dim(X)[2]){
         acs.obj=cbind(acs.obj, FUN(X[,i], ...))
      }}
      }
      if (identical(MARGIN, 2)){     # apply col-wise
        acs.obj=FUN(X[1,], ...)
        if (dim(X)[1]>1){
          for (i in 2:dim(X)[1]){
          acs.obj=rbind(acs.obj, FUN(X[i,], ...))
        }}
      }
 # I think the next part works, except it fails because the stuff above
 # doesn't like single rows or single columns...
 #    if (identical (MARGIN, c(1,2))){
 #      acs.obj=apply(apply(X, MARGIN=1, FUN=FUN), MARGIN=2,
 #      FUN=FUN)}

# or maybe this...?
     if (all(MARGIN==c(1,2))){
       acs.obj=FUN(apply(X, MARGIN=2, FUN=FUN, ...), ...)}
     acs.obj}
     

if (!isGeneric("apply")) {
  setGeneric("apply", def=function(X, MARGIN, FUN, ...){standardGeneric("apply")})}else{}
setMethod("apply", signature="acs", def=function(X, MARGIN, FUN, ...){.apply.acs(X, MARGIN, FUN, ...)})

if (!isGeneric("acs.colnames<-")) {
          setGeneric("acs.colnames<-", def=function(x, value){standardGeneric("acs.colnames<-")})}else{}
        
        #  setGeneric("acs.colnames<-", function(x, value) standardGeneric("acs.colnames<-"))
          
        setReplaceMethod(f="acs.colnames", signature="acs",
                         definition=function(x,value){
                           x@acs.colnames=value
                           x=.acs.dimnames(x)
                           x@modified=T
                           validObject(x)
                           return(x)
                         })
        
        if (!isGeneric("geography<-")) {
          setGeneric("geography<-", def=function(object, value){standardGeneric("geography<-")})}else{}
        
        setReplaceMethod(f="geography", signature="acs",
                         definition=function(object,value){
                           object@geography=value
                           object=.acs.dimnames(object)
                           object@modified=T
                           validObject(object)
                           return(object)
                         })
        
        
        if (!isGeneric("endyear<-")) {
          setGeneric("endyear<-", function(object, value) standardGeneric("endyear<-"))}else{}
        
        setReplaceMethod(f="endyear", signature="acs",
                         definition=function(object,value){
                           warning(paste(
                                         "Changing value of endyear from ",
                                         endyear(object),
                                         " to ",
                                         value,
                                         ".\nThis is an unusual thing to do, unless the original value was incorrect.\nAlso changing value of currency.year to",
                                         value,
                                         ", without converting currency values.\nPlease see ?endyear and ?currency.year for more information",
                                         sep="")
                                   , call.=F) 
                           object@endyear=as.integer(value)
                           object@currency.year=as.integer(value)
                           object@modified=T
                           validObject(object)
                           return(object)
                         })
        
if (!isGeneric("span<-")) {
  setGeneric("span<-", function(x, value) standardGeneric("span<-"))}else{}
        
setReplaceMethod(f="span", signature="acs",
                 definition=function(x,value){
                   warning(paste("Changing value of span from ",
                                 span(x), " to ", value, ".\nThis is an unusual
                         thing to do, unless the original value was
                         incorrect.\nSee ?acs for more information",
                                 sep=""), call.=F) 
                   x@span=as.integer(value)
                   x@modified=T
                   validObject(x)
                   return(x)
                 })
       
if (!isGeneric("acs.units<-")) {
  setGeneric("acs.units<-", function(x, value) standardGeneric("acs.units<-"))}else{}
        
setReplaceMethod(f="acs.units", signature="acs",
                 definition=function(x,value){
                   x@acs.units=factor(value, levels=.acs.unit.levels)
                   x@modified=T
                   validObject(x)
                   return(x)
                 })

if (!isGeneric("currency.year<-")) {
  setGeneric("currency.year<-", def=function(object, value){standardGeneric("currency.year<-")})}else{}

setReplaceMethod(f="currency.year", signature="acs",
                 definition=function(object, value){
                   currency.convert(object, rate="auto", newyear=value)
                 })
      
currency.convert=function(object, rate="auto", newyear=NA_integer_, verbose=F){
  if (rate=="auto"){
    .env=environment()
    data("cpi", envir=.env)
    new.rate="cpi[as.character(newyear)]"
    new.rate=eval(parse(text=new.rate))
    curr.rate="cpi[as.character(currency.year(object))]"
    curr.rate=eval(parse(text=curr.rate))
    rate=new.rate/curr.rate}
  dollar.cols=which(acs.units(object)=="dollars")
  if (verbose) {
    if (!missing(newyear)){
      output=c(paste("CPI (base 1982-84) for ", currency.year(object),
  " = ", curr.rate, sep=""),
      "\n",
      paste("CPI (base 1982-84) for ", newyear, " = ", new.rate, sep=""),
      "\n",
      paste("$1.00 in ", currency.year(object), " dollars = $", round(rate,2), " in ", newyear, " dollars.", sep=""),
        "\n")} else {output=c(paste("$1.00 in ",
                       currency.year(object),
                       " dollars = $",
                       round(rate,2),
                       " in converted dollars.",
                       sep=""),
                       "\n")}
    output=c(output, "Converting the following columns:", "\n",
  paste(acs.colnames(object)[dollar.cols], "\n", sep=""))
    warning(output, call.=F)
  }
  for (i in dollar.cols){
    object@estimate[,i]=object@estimate[,i]*rate
    object@standard.error[,i]=object@standard.error[,i]*rate
  }
  object@currency.year=as.integer(newyear)
  object@modified=T
  validObject(object)
  return(object)
}
      
      # helper function for replacing geography or acs.colnames
      
      prompt.acs=function(object, filename=NA, name=NA,
        what="acs.colnames", geocols="all", ...){
        print("To end session, enter a blank line.")
        if (what=="geography"){
          if (geocols=="all") geocols=1:dim(geography(object))[2]
          value=geography(object)
          for (j in geocols){
            for (i in 1:dim(geography(object))[1]){
              line.input=readline(prompt=paste("Change ", value[i,j]," to: \n", sep=""))
              if (line.input=="") {break} else {value[i,j]=line.input}
            }
          }
        }
        else if (what=="acs.colnames"){
          value=acs.colnames(object)
          for (i in 1:length(acs.colnames(object))){
            line.input=readline(prompt=paste("Change ", value[i]," to: \n", sep=""))
            if (line.input=="") {break} else {value[i]=line.input}
          }
        } else if (what=="acs.units"){
          value=acs.units(object)
          input=rep("", length(value))
          print("Type [c]ount, [d]ollars, [p]roportion, [r]atio, or [o]ther.")
          for (i in 1:length(value)){
            line.input=readline(prompt=paste(acs.colnames(object)[i], " is currently in these units: ", value[i],".  Change to what units?: (c,d,p,r,o)\n", sep=""))
            if (line.input=="") {break} else {input[i]=line.input}            
          }
          for (i in .acs.unit.levels){
            value[input==substr(start=1, stop=1,i)]=i}
        } else{
          value=NA
          warning(paste("prompt can only prompt for \"geography\", \"acs.units\", or \"acs.colnames\", not \""
                        , what, "\"", sep=""))}
        value
      }

setMethod("plot",
    signature(x = "acs"),
    function (x, conf.level=.95, err.col="red", err.lwd=1,
    err.pch="-", err.cex=2, err.lty=2, x.res=300, labels="auto",
    by="geography", true.min=T, ...) 
    {
      # internal plot function to plot individual x-y plots with conf
      # intervals, assume that either i or j or length 1
      plot.xy.acs=function(x, object, conf.int.upper, conf.int.lower, estimates, labels, xlab, ylab, ...){
        if(missing(xlab)) xlab=""
        if(missing(ylab)) ylab=""       
        plot(rep(x,2), c(conf.int.upper, conf.int.lower), type="n", xaxt="n", xlab=xlab, ylab=ylab,...)
        axis(side=1, labels=labels, at=x, ...)
        lines(x=matrix(c(x,x,rep(NA, length(x))), ncol=length(x), byrow=T), matrix(c(conf.int.lower, conf.int.upper, rep(NA, length(x))), ncol=length(x), byrow=T), col=err.col, lwd=err.lwd, lty=err.lty)
        points(x, conf.int.upper, pch=err.pch, cex=err.cex, col=err.col)
        points(x, conf.int.lower, pch=err.pch, cex=err.cex, col=err.col)
        points(x, estimates,...)
      }
      acs.density.plot=function(x, type="l", xlim,
        xlab=acs.colnames(x), ylab="Density Distribution",
        conf.level, col="black", err.col, err.lwd, err.lty,
        x.res, ...){
        est=estimate(x)
        err=standard.error(x)
        if (missing(xlim)) xlim=c(est-(4*err), est+(4*err))
        x.vals=seq(from=xlim[1], to=xlim[2], by=(xlim[2]-xlim[1])/x.res)
        plot(x.vals, dnorm(x.vals, mean=est, sd=err), type=type,
        xlab=xlab, ylab=ylab, col=col,...)
      if (conf.level) {abline(v=qnorm(mean=est, sd=err, p=c(((1-conf.level)/2), (1-((1-conf.level)/2)))), col=err.col, lwd=err.lwd, lty=err.lty)}
        }
      i=dim(x)[1]
      j=dim(x)[2]
      if (length(x)==1) {
        acs.density.plot(x, conf.level=conf.level, err.col=err.col,
        err.lwd=err.lwd, err.lty=err.lty, x.res=x.res, ...)
      } else if (i == 1 | j == 1){
        con=confint(x, level=conf.level)
        conf.int.upper=NA
        conf.int.lower=NA
        estimates=NA
        if (i == 1) {                               # one row
          if (identical(labels, "auto")) labels=acs.colnames(x)
          for (k in 1:j){
            conf.int.upper[k]=as.numeric(con[[k]][2])
            if (true.min==T) {
              conf.int.lower[k]=as.numeric(con[[k]][1])
            } else {if (true.min==F) {true.min=0}
                    conf.int.lower[k]=max(true.min, as.numeric(con[[k]][1]))}
            estimates[k]=estimate(x)[1,k]
            }}
          else        {
            if (identical(labels, "auto")) labels=geography(x)[[1]]
          for (k in 1:i){                       # one column
            conf.int.upper[k]=as.numeric(con[[1]][k,2])
            if (true.min==T) {
              conf.int.lower[k]=con[[1]][k,1]
            } else {if (true.min==F) {true.min=0}
                    conf.int.lower[k]=max(true.min, con[[1]][k,1])}
            estimates[k]=estimate(x)[k,1]
          }}
        plot.xy.acs(x=1:max(i,j), object=x, conf.int.upper=conf.int.upper, conf.int.lower=conf.int.lower, estimates=estimates, labels=labels,...)
      } else {
        if (by=="geography"){
          par(mfrow=c(i, 1))
        for (k in 1:i){
          plot(x[k,], sub=geography(x)[k,1], conf.level=conf.level, err.col=err.col, err.lwd=err.lwd, err.pch=err.pch, err.cex=err.cex, err.lty=err.lty, labels=labels,...)
        }        
        } else if (by=="acs.colnames") {
         par(mfrow=c(1,j))
        for (k in 1:j){
          plot(x[,k], sub=acs.colnames(x)[k], conf.level=conf.level,
      err.col=err.col, err.lwd=err.lwd, err.pch=err.pch,
      err.cex=err.cex, err.lty=err.lty, labels=labels,...)
       }        
       }
      }
    }
          )
