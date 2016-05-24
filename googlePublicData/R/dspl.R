.getFilesNames <- function(path, sep=";") {
################################################################################
# DESCRIPTION:
# As a function of a path and a file extension, gets all the file names there.
################################################################################  
  # Listing files
  files <- list.files(path=path, full.names=T)
  
  # Match pattern
  if (sep %in% c("xls","xlsx")) {
    pat <- "[.xls][a-z]{1}$"
  }
  else {
    if (sep %in% c(";", ",")) sep <- "csv"
    else if (sep == "tab") sep <- "tab"
    pat <- paste('[.',sep,']$',sep='')
  }
  
  # Matching excluding windows "backup" files
  valid <- grepl(pattern=pat, files) & !grepl("^[~$].*", files)
  files <- files[valid]
  
  nfiles <- length(files)
  if (nfiles==0) {
    stop(nfiles, ' files found.\nCannot\' continue') 
  }
  else {
    message(nfiles, ' files found...')
  }
  return(files)
}

genMoreInfo <- function(path, encoding=getOption("encoding"), sep=";", 
                        output=NA, action="merge", dec=".") {
################################################################################
# DESCRIPTION:
# Reads .csv and .xls(x) files, outputs a descriptive dataframe of the data and
# builds a config file.
################################################################################  
  
  options(stringsAsFactors=F)
  
  # Checks if the path exists
  .checkPath(path, "input")
  
  # Generates the filelist acording to an specific sep
  files <- .getFilesNames(path, sep)
  
  # Reads and analices the files
  x <- seekTables(files=files, encoding=encoding, sep=sep, dec=dec)
    
  # Extracts the unique list of variables
  x <- unique(
    # subset(x,subset=type != 'date' & is.dim.tab==F, 
           # select=c(id, label))
    x[
      with(x, type!='date' & !is.dim.tab),
      colnames(x) %in% c('id','label'),
      FALSE]
    )
  
  x <- data.frame(cbind(x, description=NA, topic=NA, url=NA, totalName=NA, 
                        pluralName=NA), stringsAsFactors=F)
  
  # Prints the result
  if (is.na(output)) {
    return(x) 
  } 
  
  # Case of mergin
  if (action == 'merge') {
    target.exists <- file.exists(output)
    
    if (target.exists) {
      # Merge of the tables
      x0 <- read.table(file=output, na.strings='NA', sep='\t')
      x <- subset(x, !(x$id %in% x0$id))
      x <- rbind(x0,x)
    }
    else {
      warning('The file ',output,'doesn\'t exists. It will be created')
    }
    ERR<-try(write.table(x, file=output, quote=F, na="NA", sep='\t'))
    if (class(ERR)!='try-error') {
      message("DSPL Configfile written correctly at\n",normalizePath(output))
    }
    else {
      stop("An error has occurred during the file writing at:\n",normalizePath(output))
    }
  }
  
  # Case of replacing
  else if (action %in% c('replace','merge')) {
    ERR<-try(write.table(x, file=output, quote=F, na="NA", sep='\t'))
    if (class(ERR)!='try-error') {
      message("DSPL Configfile written correctly at:\n",normalizePath(output))
    }
    else {
      stop("An error has occurred during the file writing at:\n",normalizePath(output))
    }
  }
}

seekTables <- function(files, encoding, sep, output = NA, replace = T, dec) {
################################################################################
# DESCRIPTION:
# Reads .csv and .xls(x) files, exports them as csv and outputs a descriptive da
# taframe. Also determinates which field is dim or metric.
################################################################################  
  # Checking if xls option is Ok
  #if (sep == 'xls') load.xlsx <- require(XLConnect) else load.xlsx <- TRUE
  
  #if (!load.xlsx) stop("In order to read MS Excel files ",
  #                     "you need to get the package \'XLConnect\' first.")  
  
  # Timeframe metrics
  metrics <- matrix(c(
    'dia','day','semana','week','trimestre','quarter', 'mes','month','ano', 
    'year', 'year','year','month','month'), ncol = 2, byrow=T)

  FUN <- function(x,y,z) {           
           # In the case of csv, tab 
           if (!(sep %in% c("xls", "xlsx"))) {
             
             data <- read.table(
               x, sep=sep, strip.white=T, encoding=y, fill=T,
               dec=z, header=T
             )
             
           } else {
              data <- readWorksheetFromFile(x, header=T, sheet=1)
           }
           
           cols <- colnames(data)
           
           cols <- gsub(".", " ", cols, fixed = T)
           cols <- gsub("^[[:space:]]*|[[:space:]]*$", "", cols)
           cols <- gsub("[^[:graph:]][[:space:]]*"," ", cols)
           
           colnames(data) <- cols
           
           fnames <- gsub("\\.[[:alpha:]]*$","",x)
           fnames <- gsub(".*(/|\\\\)", "", fnames)

           fnames <- rep(fnames, length(cols))
           
           # Builds descriptive matrix
           
           var <- data.frame(
             id=.cleantext(cols),
             label=cols,
             type=.fixType(unlist(lapply(data, typeof))),
             slice=fnames
           )
           
           # Creates a new column of metric vs dimm
           var <- cbind(var, concept.type='metric')
           
           var[var[,1] %in% metrics[,1], 5] <- 'dimension' # If time
           var[var[,3] %in% c('string'),5] <- 'dimension' # If string
           var[var[,1] %in% c('longitud','latitud','colour'),5] <- NA # If string
           var[var[,1] %in% metrics[,1], 3] <- 'date' # If date
            
           # Identifies which dataset is a dimension dataset
           var <- cbind(var, is.dim.tab = F)
           if (all(var[,3] != 'date')) {var['is.dim.tab'] <- T}
                
           # Replaces the date-time colnames for those proper to PDE
           for (i in 1:NROW(metrics)) {
             cols <- gsub(metrics[i,1], metrics[i,2],.cleantext(cols), fixed = T)
            }
           
           var['id'] <- cols
           
           # In the case of output, it creates a new folder
           if (!is.na(output)) {
             colnames(data) <- cols
             
             # Sorts the data acording to dimensional concepts (not date)
             ord <- var[var[,5]=='dimension' & var[,3] != 'date',1]
             if (length(ord)!=0) data <- data[do.call(order,data[ord]),]
             
             # Writes the data into csv files
             con <- file(paste(output,'/',var[1,4],'.csv',sep=''), encoding="UTF-8")
             write.table(
               x=data, file=con, na='', sep=',',quote=F,row.names=F,dec='.'
               )

             message(
               gsub(".*(/|\\\\)","",x)," analized correctly, ordered by ", ord,
               " and exported as csv"
               )
           }
           else {
             message(x,"analized correctly")
           }
           return(var)
         }
  
  # Puts it all into a single matrrx
  vars <- do.call('rbind', lapply(files, FUN, y=encoding, z=dec))
  
  # Identifies where are the correspondant tables for each dimension  
  vars <- cbind(vars, dim.tab.ref = NA)
  for (i in 1:NROW(vars)) {
    if (vars$concept.type[i] == 'dimension' & vars$type[i] != 'date' & 
      vars$is.dim.tab[i] != 'TRUE') {
      delta <- try(vars[vars[,1]==vars[i,1] & vars[,6] == 'TRUE',4])
      
      if (length(delta) == 0) {
        stop('The dimension ', vars$label[i],' needs a code (key) tab')
      } else {
        vars[i,7] <- try(delta)
      }
    }
  }
  
  # Setting which concepts extends geo:location
  
  geo <-vars[with(vars, is.dim.tab == T), colnames(vars) %in% c("id", "slice"), FALSE]
  geo <- by(
    geo, geo$slice, 
    function(x) all(c('longitude', 'latitude') %in% unlist(x['id']))
  )
  
  geo <- names(which(geo))
  vars <- cbind(vars, extends.geo=F)
  vars$extends.geo[which(vars$dim.tab.ref %in% geo)] <- T
  
  return(vars)
}

.getMoreInfo <- function(source,target, encoding="unknown", id='id', label='label') {
################################################################################
# Reads from a tab file generated by genMoreInfo as a complement info to concepts
################################################################################  
  
  # Partial fix (should work, need to see
  # http://stackoverflow.com/questions/23475309/in-r-is-it-possible-to-suppress-note-no-visible-binding-for-global-variable)
  # type <- is.dim.tab <- id <- label <- freq <- NULL
  
  if (!is.na(source)) {
    
    options(stringsAsFactors=F)
    
    # If a file, reads the moreinfo file
    if (class(source) != "data.frame") {
      source <- read.table(file=source, header=T, sep='\t', na.strings='NA',
                             encoding=encoding)
    }
    
    # Cleans up the content
    source$topicid <- .cleantext(source$topic)
    
    # Checks if there are at least moreinfo concepts as path concepts, where path
    # concepts are those captured from the files at pde(path=...)
    ntarget <- NROW(target)
    nsource <- NROW(source)
    if (ntarget < nsource) {
      warning('The number of moreinfo concepts (',nsource,') vs concepts found at path (',
              ntarget,') differ. Concepts at moreinfo won\'t appear at the metadata.')
    }
    else if (ntarget > nsource) {
      stop('Can\'t continue, there are more concepts at the path than at moreinfo')
    }
    
    # Merges source and target
    source <- source[which(source[[id]] %in% target$id),] # subset(source, id %in% target$id)
    target$label <- source[[label]] # source$label
    source <- source[,colnames(source)!=label,FALSE] #subset(source, select=-label)
    target <- merge(target, source, by=c('id'))
  }
  return(target)
}

.fixType <- function(x) {
################################################################################
# Transforms the datatypes in function of allowed datatypes for DSPL
################################################################################
  replace <- matrix(c('logical', 'integer', 'double', 'complex', 'character',
                      'boolean', 'integer', 'float', 'float', 'string'), ncol =2)
  for (i in 1:5) x <- as.character(gsub(replace[i,1], replace[i,2],x, fixed=T))
  return(x)
}

.cleantext <- function(x) {
################################################################################
# Adapts labels to IDs
################################################################################  
  
  sym <- matrix(c("$", "money", "??", "grad", "#", "n", "%", "pcent",
                  " ", "_", ".", "", ",", "_", ";", "_", ":", "_",
                  "(","_",")","_"), 
                ncol = 2, byrow = T)
  
  x <- tolower(x)
  
  x <- iconv(x, to="ASCII//TRANSLIT")
  
  for (i in 1:NROW(sym)) {
    x <- gsub(sym[i,1], sym[i,2], x, fixed = T)
  }
  
  # Extracts 
  x <- 
    
  return(x)
}

.addInfo <- function(nodename,values,parent,lang) {
################################################################################
# Function to add information and provider nodes y multiple languages
################################################################################  
  newXMLNode(nodename, parent=parent, sapply(values, function(x) {
    newXMLNode('value',attrs=c('xml:lang'=lang[which(values==x)]), x,
               suppressNamespaceWarning=T)}))
}

.addTopics <- function(nodename, values, parent, lang) {
################################################################################
# Function to create and populate the topics
################################################################################  
  values <- unique(values)
  values <- subset(values, !is.na(values[,1]))
  colnames(values) <- c('label','id')
  
  fun <- function(x, ...) {apply(x, MARGIN = 1,...)}
  
  fun(values, FUN=
    function(x) {
      tempnode0 <- newXMLNode(nodename, parent=parent, attrs=c(x['id']))
      tempnode1 <- newXMLNode('info', parent=tempnode0)
      tempnode2 <- newXMLNode('name', parent=tempnode1)
      newXMLNode('value', parent=tempnode2, attrs=c('xml:lang'=lang[1]),
                 x['label'], suppressNamespaceWarning=T)
    }
  )
}

.addConcepts <- function(val,parent,lang) {
################################################################################
# Function to create and populate the concepts
################################################################################  
  colnames(val)[3] <- 'ref'
  colnames(val)[7] <- 'geo'
  
#   if (NCOL(val) > 1) {
#     fun <- function(x, ...) {apply(x, MARGIN = 1,...)}
#   } else {
#     fun <- function(x, FUN) {FUN(x)}
#   }
#   
#   fun(val, .addConcepts_sub)
  if (ncol(val)>1) apply(val,.addConcepts_sub,MARGIN = 1,parent=parent,lang=lang)
  else .addConcepts(val,parent=parent,lang=lang)
}

.addConcepts_sub <- function(x,parent,lang) {
  function(x) {
    x['is.dim.tab'] <- gsub(' ','',x['is.dim.tab'])
    x['geo'] <- gsub(' ','',x['geo'])
    
    if (x["ref"] == 'string') {ATT <- c(x['id'], extends='entity:entity')}
    if (x['ref'] == 'string' && x['geo']=='TRUE') {ATT <- c(x['id'], extends='geo:location')}
    if (x['ref'] != 'string') {ATT <- c(x['id'])}
    
    # in the case of not being a dimensional concept
    if (x["concept.type"]!="dimension") {
      tempnode0 <- newXMLNode('concept', attrs=ATT, parent=parent)
      
      # Starts Adding info
      tempnode1 <- newXMLNode('info', parent=tempnode0)
      
      tempnode2 <- newXMLNode('name', parent=tempnode1)
      
      # Adds a description
      if (!is.na(x['description'])) {
        description <- as.character(x['description'])
        tempnode3 <- newXMLNode('description', parent=tempnode1)
        newXMLNode('value', parent=tempnode3, attrs=c('xml:lang'=lang[1]), description)
      }
      
      # URL node
      if (!is.na(x['url'])) {
        url <- as.character(x['url'])
        tempnode4 <- newXMLNode('url', parent=tempnode1)
        newXMLNode('value', parent=tempnode4, attrs=c('xml:lang'=lang[1]), url)
      }
      
      # Here should start the multilanguage loop
      newXMLNode('value', parent=tempnode2, attrs=c('xml:lang'=lang[1]),
                 suppressNamespaceWarning=T, x['label'])
      
      # Adds a topic category
      if (!is.na(x['topicid'])) {
        topicref <- as.character(x['topicid'])
        newXMLNode('topic', parent=tempnode0, attrs=c(ref=topicref))
      }
      
      # Adds the data type specification
      newXMLNode('type', attrs=c(x['ref']), parent=tempnode0)
      
      # in the case of being a dimensional concept
    } else {
      tempnode0 <- newXMLNode('concept', attrs=ATT, parent=parent)
      
      # Starts adding info
      tempnode1 <- newXMLNode('info', parent=tempnode0)
      
      # Name Node
      tempnode2 <- newXMLNode('name', parent=tempnode1)
      newXMLNode('value', parent=tempnode2, attrs=c('xml:lang'=lang[1]),
                 suppressNamespaceWarning=T, x['label'])
      
      # Description node
      if (!is.na(x['description'])) {
        description <- as.character(x['description'])
        tempnode3 <- newXMLNode('description', parent=tempnode1)
        newXMLNode('value', parent=tempnode3, attrs=c('xml:lang'=lang[1]), description)
      }
      
      # URL node
      if (!is.na(x['url'])) {
        url <- as.character(x['url'])
        tempnode4 <- newXMLNode('url', parent=tempnode1)
        newXMLNode('value', parent=tempnode4, attrs=c('xml:lang'=lang[1]), url)
      }
      
      # Plural name node
      if (!is.na(x['totalName'])) {
        pluralName <- as.character(x['pluralName'])
        tempnode6 <- newXMLNode('pluralName', parent=tempnode1)
        newXMLNode('value', parent=tempnode6, attrs=c('xml:lang'=lang[1]), pluralName)
      }        
      
      # Total name node
      if (!is.na(x['totalName'])) {
        totalName <- as.character(x['totalName'])
        tempnode5 <- newXMLNode('totalName', parent=tempnode1)
        newXMLNode('value', parent=tempnode5, attrs=c('xml:lang'=lang[1]),totalName)
      }
      
      newXMLNode('type', parent=tempnode0, attrs=c(x['ref']))
      newXMLNode('table', parent=tempnode0, attrs=
                   c(ref=paste(x['dim.tab.ref'],'_table',sep='')))
      
    }
  }
}

.addSlices <- function(tableid, sliceatt, parent) {
################################################################################
# Function to create and populate the slices
################################################################################
  colnames(sliceatt)[1] <- 'concept'
  by(data=sliceatt, INDICES=tableid,FUN=
    function(x) {
      newXMLNode(name='slice', attrs=c(id=paste(x$slice[1],'_slice',sep='')),
                 parent=parent, apply(x, MARGIN = 1,FUN=
                   function(z){
                     #z <- as.character(z)
                     # In the case of dates-time
                     if (z['type'] == 'date') {
                       newXMLNode(name=z['concept.type'], 
                                  attrs=c(concept=paste('time:',z['concept'],sep='')))
                       # Otherwise
                     } else {
                       newXMLNode(name=z['concept.type'], attrs=c(z['concept']))
                     }}), newXMLNode('table', attrs=c(ref=paste(x$slice[1],'_table',sep=''))))
    }
    )
}

.addTables <- function(tableid, tableatt, parent, format) {
################################################################################
# Function to create and populate the tables
################################################################################
  by(data=tableatt, INDICES=tableid,FUN=
    function(x) {
      newXMLNode(name='table', attrs=c(id=paste(x$slice[1],'_table',sep='')),parent=
        parent, apply(X=x, 
                      MARGIN = 1, FUN=
                        function(z){
                          if (z['type'] == 'date') {
                            newXMLNode(name='column', attrs=c(z['id'], z['type'], format=format))  
                          } else {
                            newXMLNode(name='column', attrs=c(z['id'], z['type']))
                          }}), newXMLNode(name='data', newXMLNode('file', attrs=c(format=
                            'csv', encoding='utf8'),paste(x$slice[1],'.csv',sep='')))
                 )
    }
     )
}

dspl <- function(
################################################################################
# DESCRIPTION:
# Based on an specific folder directory, the function seeks for files that match
# the specified extension (csv, tab, xls, xlsx), reads the column names, guesses
# the datatype, analyces data structure and outputs a dspl metadata file includi
# ng the corresponding csv files.
#
# VARIABLES:
# - path: Full path to the folder where the tables are saved.
# - Output: FUll path to the output folder (pde will create a subfolder call r_dspl).
# - replace: In the case of defining output, replaces the files.
# - targetNamespace:
# - timeFormat: The corresponding time format of the collection.
# - lang: A list of the languages supported by the dataset.
# - name: The name of the dataset.
# - description: The dataset description.
# - url: The corresponding URL for the dataset.
# - providerNAme
# - providerURL
# - sep: The extension of the tables in the 'path' folder.
# - encoding: The char encoding of the input tables.
################################################################################  
  path,
  output = NA,
  replace = F,
  targetNamespace = "",
  timeFormat = "yyyy",
  lang = c("es", "en"),
  name = NA,
  description = NA,
  url = NA,
  providerName = NA,
  providerURL = NA,
  sep = ";",
  dec=".",
  encoding = getOption("encoding"),
  moreinfo = NA
  ) {
  options(stringsAsFactors=F)
  
  # Initial checks
  description <- ifelse(!is.na(description),description,'No description')
  name <- ifelse(!is.na(name),name,'No name')
  providerName <- ifelse(!is.na(providerName),providerName,'No provider')
  
  # Checking if output path is Ok
  if (!is.na(output)) {
    temp.path <- tempdir()
    .checkPath(path, "input")
  }
  else temp.path <- NA
  if (!all(is.na(moreinfo))) {
    if (class(moreinfo) != "data.frame") .checkPath(moreinfo, "input")
  }
  
  # Checking timeFormat
  timeOk <- checkTimeFormat(timeFormat)
  if (!timeOk) stop("Undefined \'joda-time\' definition ", timeFormat,
                    "\nFor more information checkout\n",
                    "http://joda-time.sourceforge.net/api-release/org/joda/time/format/DateTimeFormat.html")
  # Gets the filenames
  files <- .getFilesNames(path, sep)
    
  # Variables Lists and datatypes
  vars <- seekTables(files=files, encoding=encoding, sep=sep, 
                     output=temp.path, replace=replace, dec=dec)
  
  dims <- vars[
    with(vars, concept.type=='dimension'),
    colnames(vars) %in% c('id', 'slice', 'concept.type'), FALSE
    ]
  #subset(vars, concept.type=='dimension', select=c(id, slice, concept.type)) 
  
  # Identifying if there is any duplicated slice
  with(dims, .checkSlices(dims=id, by=slice))
  
  vars <- .checkDuplConcepts(concepts=vars)
  
  # Creates a unique concept list
  varConcepts <- unique(
    # subset(vars,subset=type != 'date' & is.dim.tab==F, select=-slice) 
    vars[
      with(vars, type!='date' & is.dim.tab==FALSE),
      colnames(vars)!='slice',FALSE
      ]
    )
  
  # Checks if there is a moreinfo file
  varConcepts <- .getMoreInfo(source=moreinfo, target=varConcepts, "")
  
  # Armado de xml
  archXML <- newXMLDoc()
  dspl <- newXMLNode(name='dspl', doc=archXML, attrs=c(
    targetNamespace=targetNamespace), namespace=c(xsi=
      'http://www.w3.org/2001/XMLSchema-instance'))
  
  tmpattrs <- xmlAttrs(dspl)
  try(tmpattrs['xmlns'] <- 'http://schemas.google.com/dspl/2010', silent=T)
  
  # Definiciones dspl
  
  imports <- c('quantity', 'entity', 'geo', 'time', 'unit')
  
  sapply(imports,
         function(x) {
           newXMLNamespace(node=dspl, prefix=x,namespace=paste(
             'http://www.google.com/publicdata/dataset/google/',x,sep=''))
         })
  # Concepts import lines
  newXMLCommentNode('Concepts imports', parent=dspl)
  
  imports <- paste(
    "http://www.google.com/publicdata/dataset/google/",
    imports, sep = '')
  
  sapply(X = imports,
         FUN = function(x) newXMLNode(attrs=c(namespace=x), name = 'import',
                                      parent = dspl))
  # INFO
  newXMLCommentNode('Info lines', parent=dspl)
  info <- newXMLNode('info', parent = dspl)
  .addInfo('name', name, info, lang)
  .addInfo('description', description, info, lang)
  if (!is.na(url)) newXMLNode('url', newXMLNode('value', url), 
                              parent = info)
  
  # PROVIDER
  newXMLCommentNode('Data Provider', parent=dspl)
  provider <- newXMLNode('provider', parent = dspl)
  .addInfo('name', providerName, provider, lang)
  if (!is.na(providerURL)) newXMLNode('url', newXMLNode('value', providerURL), 
                                      parent = provider)
  
  # TOPICS
  if ("topicid" %in% colnames(varConcepts)) {
    test <- !all(is.na(varConcepts$topicid))
  }
  else {
    test <- F
  }
  if (test) {
    newXMLCommentNode('Topics definition', parent=dspl)
    topics <- newXMLNode('topics', parent=dspl)
    .addTopics('topic', varConcepts[c('topic', 'topicid')], topics, lang)
  }
    
  # CONCEPTS
  newXMLCommentNode('Concepts Definitions', parent=dspl)
  concepts <- newXMLNode('concepts', parent = dspl)
  .addConcepts(varConcepts,concepts, lang)
  
  # SLICES
  newXMLCommentNode('Slices Definitions', parent=dspl)
  slices <- newXMLNode('slices', parent = dspl)
  .addSlices(
    tableid  = vars[with(vars, !is.dim.tab),, FALSE][['slice']],
    sliceatt = vars[with(vars, !is.dim.tab),, FALSE], # subset(vars, is.dim.tab != T),
    parent   = slices
    )
  
  # TABLES
  newXMLCommentNode('Tables Definitios', parent=dspl)
  tables <- newXMLNode('tables', parent = dspl)
  .addTables(tableid=vars$slice,tableatt=vars,parent=tables, format=timeFormat)  
  
  # Building ouput
  .dimtabs  <- unique(vars[with(vars, is.dim.tab),,FALSE][['slice']]) #subset(vars, subset=is.dim.tab, select=slice))
  .slices   <- unique(vars[['slice']]) #subset(vars, select=slice))
  .concepts <- unique(vars[['label']]) #subset(vars, select=label))
  .dims <- unique(  
#     subset(
#       vars, 
#       subset=is.dim.tab & !(label %in% c('name', 'latitude','longitude','colour')),
#       select=label
#       )
    vars[ 
      with(vars, is.dim.tab & !(label %in% c('name', 'latitude','longitude','colour'))),
      colnames(vars)=='label',
      FALSE
    ]) 
  
  lapply(c(.dimtabs, .slices, .concepts, .dims), function(x) names(x) <- 'Name')
  
  pde.statistics <- matrix(
    c(
      NROW(.slices),
      NROW(.concepts),
      NROW(.dims)
    ), ncol=3)
  
  colnames(pde.statistics) <- c('slices','concepts','dimentions')
  
  result <- structure(.Data=
    list(
        saveXML(archXML, encoding="UTF-8"), vars, .dimtabs, .slices, .concepts,
        .dims, pde.statistics
    ), .Names=c('dspl', 'concepts.by.table', 'dimtabs', 'slices', 'concepts',
                  'dimentions','statistics'))
  class(result) <- c("dspl")
  
  # If an output file is specified, it writes it on it
  if (is.na(output)) {
    return(result)
  
  } else {
    path <- paste(temp.path,'/metadata.xml',sep='')
    print.dspl(x=result, path=path, replace=replace, quiet=T)
    
    # Zipping the files
    tozip <- list.files(temp.path, full.names=T, pattern="csv$|xml$")
    zip(output ,tozip,flags='-r9jm')
    
    message("Metadata created successfully at:\n", normalizePath(output))
  }
}
