aggPopLong = function(popdata, poplong, sex, area.scale) {

# called by getSMR df numeric and getSMR df glm  
  
if('sex' %in% names(poplong)){
  sex = paste("^", paste(sex, collapse='$|^'), "$", sep='')
  poplong= poplong[grep(sex, poplong$sex, ignore.case=TRUE),]
}	

poplong <- aggregate(poplong$expected, 
    list(poplong[['idForAgg']]), 
    sum, na.rm=TRUE)
rownames(poplong) = as.character(poplong[,1])
poplong=poplong[poplong[,2] > 0,]


# merge results back in to the population data
# the merge function changes the order, so can't use it.
popdata$expected = NA
rownames(popdata) = as.character(popdata[,'idForAgg'])

popdata[rownames(poplong), "expected"] = poplong[,2]

surfaceArea = grep("^surfacearea$", names(popdata), ignore.case=TRUE, value=TRUE)
if (length(surfaceArea) ) {
  popdata$expected_surfaceArea <- area.scale *popdata$expected/popdata[[surfaceArea[1] ]]
  popdata$logExpected_surfaceArea = log(popdata$expected_surfaceArea)
}

popdata$logExpected = log(popdata$expected)

popdata = popdata[,
    grep('^idForAgg$', colnames(popdata), invert=TRUE)]


popdata

}

setGeneric('getSMR', 
    function(popdata, model, casedata, 
        regionCode,
        regionCodeCases, 
        area.scale=1, 
        sex=c('m','f'), ...
) {
      standardGeneric("getSMR")
    }
)




setMethod("getSMR", 
    signature("data.frame", "numeric", "missing", "ANY"),
    function(popdata, model, 
        casedata, 
        regionCode,
        area.scale=1, 
        sex=c('m','f'), ...){
      # model is a vector of rates
      
      popdata$idForAgg = 1:nrow(popdata)
 
      # check breaks for groups, make sure they line up
      rateBreaks = getBreaks(names(model))
      popBreaks = getBreaks(names(popdata))
      
      newBreaks = getBreaks(intersect(rateBreaks$newNames, popBreaks$newNames))
      
      
      newModel = data.frame(
          age=rateBreaks$age, 
          sex=rateBreaks$sex, 
          rate=model)
      newModel = formatCases(newModel, newBreaks)
      newModel = aggregate(newModel$rate, 
          newModel[,c('age','sex')],mean,na.rm=T)
      rownames(newModel) = 
          paste(newModel$sex, 
              newModel$age, 
              sep='.')
      
      poplong = formatPopulation(popdata, breaks=newBreaks$breaks)
      
      poplong$expected = poplong$POPULATION *  
          newModel[paste(poplong$sex, poplong$age, sep='.')
              ,'x']

      aggPopLong(popdata, poplong, sex, area.scale) 
    }
)

setMethod("getSMR", 
    signature("data.frame", "glm", "missing", "ANY"),
    function(popdata, model, casedata, 
        regionCode,
        area.scale=1, 
        sex=c('m','f'), ...){
      
      popdata$id = 1:nrow(popdata)
      poplong <- formatPopulation(popdata, breaks=attributes(model)$breaks$breaks)
      
      p<-grep("^population$", names(poplong), value=TRUE, ignore.case=TRUE)  
      
      poplong[is.na(poplong[,p]),p] <- 0
      
      popBreaks = attributes(poplong)$breaks
      
      # get rid of zero populations,because they dont lead to rates of exactly zero
      ## check the class of the POPULATION column and make sure it is numeric: 
      # converted to character first in case it's a factor
      # (convert factor names to integers, not factor levels)
      
      if(class(poplong[,p]) != "numeric"){
        poplong[,p] <- as.numeric(as.character(poplong[,p]))
      }
      
      poplong=poplong[poplong[,
              grep("^population$", names(poplong), value=TRUE, ignore.case=TRUE)]>0, ]     
      #changes poplong names to be consistent with model
      agevar<-grep("^age$",names(attributes((terms(model)))$dataClasses),value=TRUE,ignore.case=TRUE)
      sexvar<-grep("^sex$",names(attributes((terms(model)))$dataClasses),value=TRUE,ignore.case=TRUE)
      yearvar<-grep("^year$",names(attributes((terms(model)))$dataClasses),value=TRUE,ignore.case=TRUE)
      
      agevar1<-grep("^age$",names(poplong),value=TRUE,ignore.case=TRUE)
      sexvar1<-grep("^sex$",names(poplong),value=TRUE,ignore.case=TRUE)
      yearvar1<-grep("^year$",names(poplong),value=TRUE,ignore.case=TRUE)
      
      if(length(agevar) & length(agevar1)){
        names(poplong[[agevar1]])=agevar
      }
      if(length(sexvar) & length(sexvar1)){
        names(poplong[[sexvar1]])=sexvar
      }
      if(length(yearvar) & length(yearvar1)){
        names(poplong[[yearvar1]])=yearvar
      }
      
      # get rid of a gender
      if (length(model$sexSubset) == 1) {
        poplong = poplong[poplong[[sexvar1]] == model$sexSubset, ]
      }
      
      
      offsetvar<- grep("^offset",names(attributes((terms(model)))$dataClasses)  ,value=TRUE, ignore.case=TRUE)
      offsetvar<- substr(offsetvar,8,nchar(offsetvar)-1)
      
      
      poplong[[offsetvar]] = log(poplong$POPULATION)
#    poplong[[sexvar]]= factor(poplong[[sexvar]])
#    poplong[[agevar]] = factor(poplong[[agevar]])
      #names(poplong) <- tolower(names(poplong))
      for (Dlevel in names(model$xlevels)) {
        alllevels = levels(poplong[[Dlevel]])
        if (!all(alllevels %in% model$xlevels[[Dlevel]])) {
          tokeep = poplong[[Dlevel]] %in% model$xlevels[[Dlevel]]
          poplong = poplong[tokeep, ]
        }
      }
      interactNA <- names(model$coefficients)[is.na(model$coefficients)]
      if (length(interactNA) > 0) {
        interact <- grep(":", interactNA, value = TRUE)
        poplong$param <- paste(paste("age", poplong$age, sep = ""),
            paste("sex", poplong$sex, sep = ""), sep = ":")
        poplong = poplong[!poplong$param %in% interact, ]
        poplong$param <- NULL
      }
      if(length(yearvar1)) {
        agg<-c(yearvar1,agevar, sexvar,offsetvar)
      }else{
        agg<-c(agevar, sexvar, offsetvar)
      }
      
      # multiply population by popScale, to make it in person years
      if(any(names(attributes(popdata))=="popScale")) {
        poplong[,offsetvar]=     poplong[,offsetvar] + log(attributes(popdata)$popScale)
      }
      
      
      poplong$expected <- predict(model, poplong,#[,agg],
          type = "response")
      
      aggPopLong(popdata, poplong, sex, area.scale) 

    }
)

# set regionCode and regionCodeCases with case data
setMethod("getSMR", 
    signature("data.frame", "ANY", "data.frame", 'missing', 'missing'),
    function(popdata, model, casedata, 
        regionCode,
        regionCodeCases, area.scale=1, 
        sex=c('m','f'), ...
    ){
      
      regionCodeCases=regionCode = intersect(names(popdata), 
          names(casedata)
      )[1]
      
      
      if(!length(regionCodeCases))
        warning('cant find region code in case data')
      
      callGeneric()
      
    }
)

# set regionCode and regionCodeCases with case data
setMethod("getSMR", 
    signature("data.frame", "ANY", "data.frame", 'character', 'missing'),
    function(popdata, model, casedata, 
        regionCode,
        regionCodeCases, area.scale=1, 
        sex=c('m','f'), ...
    ){
      
      regionCodeCases=regionCode
      
      
      if(!length(grep(regionCodeCases, names(casedata))))
        warning('cant find region code in case data')
      
      callGeneric(popdata, model, casedata, regionCode,
          regionCodeCases, area.scale, sex)
      
    }
)


setMethod("getSMR", 
    signature("data.frame", "ANY", "data.frame", "character", "character"),
    function(popdata, model, casedata, 
        regionCode,
        regionCodeCases, area.scale=1, 
        sex=c('m','f'), ...
){
# have case data, compute SMR

    popdata = callGeneric(
        popdata, model,
        area.scale=area.scale,
        sex=sex
        )
      
        if(!length(regionCode)) {
          regionCode='id'
        }

        casecol = grep("^cases$|^count$|^y$", names(casedata), value=TRUE, ignore.case=TRUE)
    if(length(casecol)>1) {
      casecol=casecol[1]
      warning("more than one column which could be interpreted as case numbers, using ", casecol)
    }
    
    if(!length(casecol)) {
      #there is no case col
      casecol = "cases"
      casedata[,casecol] = 1
    }
    
    
    casedata = casedata[
        as.character(casedata[, regionCodeCases]) %in% 
            as.character(popdata[[regionCode]]), ]
    
    casedata <- aggregate(casedata[,casecol], 
        list(casedata[,regionCodeCases]), sum)
    names(casedata) = c(regionCodeCases, "observed")
    
    popdata[['observed']] = casedata[
        pmatch(
        as.character(popdata[[regionCode]]),
        as.character(casedata[[regionCodeCases]])
        ),
        'observed']
    
    # change 0's in expected to NA, so SMR is NA
    theexpected = popdata$expected
    theexpected[theexpected==0] = NA
    
    popdata$SMR <- popdata$observed/theexpected

    popdata
    }
)


setMethod("getSMR", 
    signature("SpatialPolygonsDataFrame", 'ANY','ANY', "ANY"),
    function(
        popdata, model,  casedata, 
        regionCode,
        regionCodeCases, area.scale=1, 
        sex=c('m','f'), ...
    ) {
      
      area = !length(grep("longlat", popdata@proj4string@projargs))
      
      if (area & !("surfaceArea" %in% names(popdata) ) ) {
        popdata$surfaceArea = sapply(slot(popdata, "polygons"), slot, "area")
      }	

      popdata$idForGetSmr = 1:length(popdata)
      theSP = popdata
      popdata = data.frame(popdata)
      
      forSp <- callGeneric()

      theSP@data = forSp[pmatch(
              theSP$idForGetSmr,
              forSp$idForGetSmr
              ),]
      theSP = theSP[,grep("^idForGetSmr$", names(theSP), invert=TRUE)]    
      theSP
    }
)

setMethod("getSMR", 
    signature("data.frame", 'list', 'missing', 'ANY'),
    function(popdata, model, casedata, 
        regionCode,
        regionCodeCases, area.scale=1, 
        sex=c('m','f'), ...){
      
      keepCols = c(
          'expected','expected_surfaceArea',
          'logExpected_surfaceArea', 'logExpected'
          )
      
          
      if(!length(names(model))){
        names(model) = as.character(1:length(model))
      }
      
      result =  as.data.frame(matrix(NA, nrow(popdata), 0))
      
      
      modelFull = model
       for(Drate in names(modelFull)){
        model=modelFull[[Drate]]
        toBind = callGeneric()
        toBind = toBind[,intersect(colnames(toBind), keepCols)]
        names(toBind) = paste(
            names(toBind), Drate, sep="_"
        )
        result = cbind(result,toBind)
      }
 
      popdata = cbind(popdata, result)
      popdata
    }
)

setMethod("getSMR", 
    signature("list", 'list', 'missing', 'ANY'),
    function(popdata, model, casedata, 
        regionCode,
        regionCodeCases, area.scale=1, 
        sex=c('m','f'), ...){

    Nmax = max(c(length(popdata),length(model)))

    if(!length(names(model))){
      names(model) = as.character(round(seq(from=1,to=Nmax, len=length(model))))
    }
    if(!length(names(popdata))){
      names(popdata) = as.character(round(seq(from=1,to=Nmax, len=length(popdata))))
    }
    
    Scensus = as.numeric(names(popdata))
    Srates = as.numeric(names(model))
    if(any(is.na(Scensus)) | any(is.na(Srates)))
      warning("can't turn names of popdata or model into numbers")
    
    yearBreaks = c(
        -Inf,
        Scensus[-1] - diff(Scensus)/2,
        Inf
        )
    
    rateCut = Scensus[
        as.integer(cut(Srates, breaks=yearBreaks))
      ]

      popDataFull = popdata
      modelFull = model
      
      
    result = list()
    

    for(Dcensus in Scensus){
       popdata = popDataFull[[as.character(Dcensus)]]
      model = modelFull[
          rateCut==Dcensus
          ]
       result[[
          as.character(Dcensus)
      ]] = callGeneric()
    }
    result
    }
)

setMethod("getSMR", 
    signature("list", 'ANY', 'data.frame'),
function(popdata, model, casedata, regionCode,
    regionCodeCases, area.scale=1 , 
    years = NULL, personYears=TRUE,year.range = NULL,...){
#  lennon's stuff
  #isSP = (class(popdata[[1]]) == "SpatialPolygonsDataFrame")
  
  if (is.null(years)) {
    years = as.integer(names(popdata))
  }        
  
  yearVar = grep("year", names(attributes((terms(model)))$dataClasses) ,
      value=TRUE,ignore.case=TRUE)  # find year var in the model
  #If year is not in the model as factor
  if(length(yearVar)==0) yearVar="YEAR"
  
  caseYearVar = grep("year",names(casedata),value=TRUE,ignore.case=TRUE)
  
  
  if(personYears){
    # if year.range is missing, use year range of the cases
    if(is.null(year.range)) {
      if(!length(caseYearVar))
        warning("year.range unspecified and no year column in case data")
      year.range = range(as.numeric(as.character(casedata[,caseYearVar])), na.rm=TRUE)
    }
    
    times <- c(year.range[1], sort(years), year.range[2])
    times <- as.numeric(times)
    inter <- diff(times)/2
    nseq <- 1:length(inter) - 1
    mseq <- 2:length(inter)
    interval <- inter[mseq] + inter[nseq]
    names(interval) <- names(popdata)
  }else{
    interval<-rep(1,length(popdata))
  }
  
  names(interval) <- names(popdata)
  
  
  caseSexVar =grep("^sex$",names(casedata),value=TRUE,ignore.case=TRUE)
  if (length(model$sexSubset) == 1) warning("only one sex is being used:",model$sexSubset)
  
  
  
  result=list()
  
  area = !length(grep("longlat", popdata[[1]]@proj4string@projargs))
  
  #for(Dyear in seq(1, length(popdata))) {
  for(Dyear in names(popdata)) {
    
    #Compute area
    if(area & class(popdata[[Dyear]]) == "SpatialPolygonsDataFrame") {
      areas <- unlist(computeArea(popdata[[Dyear]]))
      popdata[[Dyear]]$sqk <- areas * area.scale
    }
    
    
    # scaling factor to convert to person years
    attributes(popdata[[Dyear]]@data)$popScale = interval[names(interval)==Dyear]
    
    # add year column
    #popdata[[Dyear]][[yearVar]] = factor(rep(as.character(names(popdata[Dyear])), length(popdata[[Dyear]][1]),
    #levels=levels(model$data[,yearVar])))
    
    popdata[[Dyear]][[yearVar]] = as.integer(Dyear)
    
    
    caseThisYear = casedata[casedata[[caseYearVar]]==Dyear,]
    
    if (length(model$sexSubset) == 1) {      
      caseThisYear = caseThisYear[caseThisYear[[caseSexVar]] == model$sexSubset, ]
    }
    
    if(dim(caseThisYear)[1]==0) caseThisYear<-NULL   
    
    cat("computing",Dyear,"\n")
    
    result[[Dyear]] = getSMR(popdata[[Dyear]], model, caseThisYear,regionCode =regionCode,
        regionCodeCases = regionCodeCases, years = years, year.range = year.range,
        area = area, area.scale = area.scale)
    
  }
  names(result) = years
  
  
  return(result)       
  
}
)

