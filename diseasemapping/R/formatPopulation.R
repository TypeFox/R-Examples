setGeneric('formatPopulation', 
    function(popdata, 
        aggregate.by=NULL, 
        breaks=NULL, ...) {
      standardGeneric("formatPopulation")
    }
)

setMethod("formatPopulation", 
    signature("SpatialPolygonsDataFrame"),
    function(popdata, ...) {
      
      callGeneric(popdata@data, ...)
      
    }                         
)

# data is a raster.  grid is ignored
setMethod("formatPopulation", 
    signature("data.frame"),
    function(popdata, aggregate.by=NULL, breaks=NULL,...) {

#popdata <- popdata@data
ageBreaks = getBreaks(names(popdata), breaks)

####reshape the popdata:
poplong = reshape(popdata,  varying=ageBreaks$oldNames, direction="long",
	v.names="POPULATION", timevar="GROUP", times = ageBreaks$newNames)
# create age and sex variables
agecol = grep("^age$", names(poplong), value=TRUE, ignore.case=TRUE)
sexcol = grep("^sex$", names(poplong), value=TRUE, ignore.case=TRUE)


if("GROUP" %in% names(poplong)) {
    if(!length(sexcol)){
 
  poplong$sex = factor(toupper(substr(poplong$GROUP, 1, 1)))}
     if(!length(agecol)){
   ageNumeric = as.numeric(substr(poplong$GROUP, 3, 4))
   poplong$age = cut(ageNumeric, ageBreaks$breaks, right=FALSE)
  }else {
  warning("no age and sex variables found or no group variable found in popdata")
  }
}


row.names(poplong)<-NULL

if(!is.null(aggregate.by)) {

  popa <- aggregate(poplong$POPULATION, poplong[, aggregate.by, drop=FALSE], 
    sum, na.rm=TRUE)

  # change x column name to 'population'
  names(popa)[names(popa)=="x"] = "POPULATION"
  names(popa)[names(popa)=="poplong[, aggregate.by]"] = aggregate.by
  poplong<-popa

}

if(length(sexcol)) poplong[,sexcol] <- toupper(poplong[,sexcol])
attributes(poplong)$breaks = ageBreaks

poplong

}
)

setMethod("formatPopulation", 
    signature("list"),
    function(popdata, aggregate.by=NULL, breaks = NULL, 
        years=as.integer(names(popdata)), year.range=NULL,  time="YEAR", 
        personYears=TRUE,...) {
      
      time<-toupper(time)
      
      #If aggregate, then see if YEAR is there or not, if so, remove it
      if(!is.null(aggregate.by)){ 
#   agg<-aggregate.by<-toupper(aggregate.by)
        agg<-aggregate.by
        byYear<- time %in% aggregate.by
        if(byYear){aggregate.by<-aggregate.by[-which(aggregate.by==time)]}
      }
      
      
      
      listpop<-lapply(popdata,formatPopulation,aggregate.by, breaks= breaks)
      
      breaks = attributes(listpop[[1]])$breaks
      
      
      listdataframe<-lapply(listpop,as.data.frame)
      #if did not aggregate, then the data frames will have differnt columns
      pop<-NULL
      for (i in 1:length(listdataframe)){
        temp<-listdataframe[[i]]
        temp[,time]<-years[i]
        pop<-rbind(pop,temp)
      }
      
      attributes(pop)$breaks = breaks
      
      if(personYears){
        if (is.null(year.range)) {
          year.range = range(pop[,time])
        }
        times <- c(year.range[1], sort(years), year.range[2])
        times <- as.numeric(times)
        inter <- diff(times)/2
        nseq <- 1:length(inter) - 1
        mseq <- 2:length(inter)
        interval <- inter[mseq] + inter[nseq]
        names(interval) <- names(table(pop[,time]))
        pop$yearsForCensus = interval[as.character(pop[,time])]
        pop$POPULATION = pop$POPULATION * pop$yearsForCensus
        pop[,time] = factor(pop[,time], levels = unique(pop[,time]))
        pop[,time] = factor(pop[,time])
        
        pop <- pop[!is.na(pop$POPULATION),  ]
        pop <- pop[pop$POPULATION > 0,  ]
      }
      
      
      
      pop
    }
)