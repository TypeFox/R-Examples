FromPositionToId <- function(react.id) {
  row <- which(LETTERS ==
                 gsub("([A-Z])[0-9]+", "\\1", react.id))
  col <- as.integer(gsub("[A-Z]([0-9]+)", "\\1", react.id))
  # as.character((row - 1) * 12 + col)
  (row - 1) * 12 + col
}

GetIds <- function(l) {
  unname(sapply(l, function(el) el$id))
}

# Gets concentrations (quantity) of each 
# dilution from XML for Roche
GetDilutionsRoche <- function(uniq.folder)
{
  cat("\nParsing Roche standards data...")
  if(!file.exists(paste0(uniq.folder,"/calculated_data.xml"))) {
    cat("NO SUCH FILE")
    return(NA)
  }
  rdml.doc <- xmlParse(paste0(uniq.folder,"/calculated_data.xml"))
  concs<-as.numeric(xpathSApply(
    rdml.doc,
    "//ns:absQuantDataSource/ns:standard",   
    namespaces = c(ns = "http://www.roche.ch/LC96AbsQuantCalculatedDataModel"),
    xmlValue))
  if(length(concs) == 0) {
    concs<-as.numeric(xpathSApply(
      rdml.doc,
      "//ns:relQuantDataSource/ns:standard",   
      namespaces = c(ns = "http://www.roche.ch/LC96RelQuantCalculatedDataModel"),
      xmlValue))
    concs.guids<-xpathSApply(
      rdml.doc,
      "//ns:relQuantDataSource/ns:standard/../ns:graphId",   
      namespaces = c(ns = "http://www.roche.ch/LC96RelQuantCalculatedDataModel"),
      xmlValue)
    names(concs) <- concs.guids
    concs <- sort(concs, decreasing=TRUE)
    positions <- 
      xpathSApply(
        rdml.doc, 
        paste0("//ns:standardPoints/ns:standardPoint/ns:position"), 
        xmlValue,
        namespaces = 
          c(ns = "http://www.roche.ch/LC96RelQuantCalculatedDataModel"))
    positions <- sapply(positions, FromPositionToId)
    dye.names <-xpathSApply(
      rdml.doc, 
      paste0("//ns:standardPoints/ns:standardPoint/ns:dyeName"), 
      xmlValue,
      namespaces = c(ns = "http://www.roche.ch/LC96RelQuantCalculatedDataModel"))
    positions.guids <- xpathSApply(
      rdml.doc, 
      paste0("//ns:standardPoints/ns:standardPoint/ns:graphIds/ns:guid"), 
      xmlValue,
      namespaces = c(ns = "http://www.roche.ch/LC96RelQuantCalculatedDataModel"))
  } else {
    concs.guids<-xpathSApply(
      rdml.doc,
      "//ns:absQuantDataSource/ns:standard/../ns:graphId",   
      namespaces = c(ns = "http://www.roche.ch/LC96AbsQuantCalculatedDataModel"),
      xmlValue)
    names(concs) <- concs.guids
    concs <- sort(concs, decreasing=TRUE)
    positions <- 
      xpathSApply(
        rdml.doc, 
        paste0("//ns:standardPoints/ns:standardPoint/ns:position"), 
        xmlValue,
        namespaces = 
          c(ns = "http://www.roche.ch/LC96AbsQuantCalculatedDataModel"))
    positions <- sapply(positions, FromPositionToId)
    dye.names <-xpathSApply(
      rdml.doc, 
      paste0("//ns:standardPoints/ns:standardPoint/ns:dyeName"), 
      xmlValue,
      namespaces = c(ns = "http://www.roche.ch/LC96AbsQuantCalculatedDataModel"))
    positions.guids <- xpathSApply(
      rdml.doc, 
      paste0("//ns:standardPoints/ns:standardPoint/ns:graphIds/ns:guid"), 
      xmlValue,
      namespaces = c(ns = "http://www.roche.ch/LC96AbsQuantCalculatedDataModel"))    
  }
  positions.table <- matrix(c(dye.names,
                              positions),
                            ncol = length(positions),
                            nrow = 2,
                            byrow = TRUE,
                            dimnames = list(c("dye.name","position"),
                                            positions.guids))
  positions.table <- positions.table[,
                                     order(match(colnames(positions.table), names(concs)))]
  positions.table <- rbind(positions.table, conc = concs)
  dyes <- unique(positions.table["dye.name",])
  dilutions <- lapply(dyes, function(dye) {
    dye.group.indecies <- which(positions.table["dye.name",] == dye)
    concs.by.dye <- concs[dye.group.indecies]
    names(concs.by.dye) <- positions.table["position",
                                           dye.group.indecies]
    concs.by.dye
  })
  if (length(dilutions) == 0) {
    cat("NONE")
    return(NULL)
  }
  names(dilutions) <- dyes
  cat("OK")
  return(dilutions)
}

GetConditionsRoche <- function(uniq.folder)
{
  cat("\nParsing Roche conditions data...")
  if(!file.exists(paste0(uniq.folder,"/app_data.xml"))) {
    cat("NO SUCH FILE")
    return(NA)
  }
  rdml.doc <- xmlParse(paste0(uniq.folder,"/app_data.xml"))
  nodes <- getNodeSet(rdml.doc,
                      "/ns:rocheLC96AppExtension/ns:experiment/ns:run/ns:react/ns:condition/..",
                      namespaces = c(ns = "http://www.roche.ch/LC96AppExtensionSchema"))
  
  reacts <- sapply(nodes,
                   function(node)  xmlAttrs(node, "id"))
  conditions <- sapply(nodes,
                       function(node)  xmlValue(node[["condition"]]))
  if (length(conditions) == 0) {
    cat("NONE")
    return(NULL)
  }
  names(conditions) <- reacts
  cat("OK")
  return(conditions)
}

GetRefGenesRoche <- function(uniq.folder)
{
  cat("\nParsing Roche reference genes data...")
  if(!file.exists(paste0(uniq.folder,"/module_data.xml"))) {
    cat("NO SUCH FILE")
    return(NA)
  }
  rdml.doc <- xmlParse(paste0(uniq.folder,"/module_data.xml"))
  
  ref <- getNodeSet(
    rdml.doc,
    "//ns:geneSettings/ns:relQuantGeneSettings",
    namespaces = c(ns = "http://www.roche.ch/LC96RelQuantGeneralDataModel"))
  
  if (length(ref) == 0) {
    cat("NONE")
    return(NULL)
  }
  
  cat("OK")
  return(ref)
}


#' Creates new instance of \code{RDML} class object
#' 
#' This function has been designed to import data from RDML v1.1 and v1.2 format files
#' 
#' @section Warning: Although the format RDML claimed as data exchange format, 
#'   the specific implementation of the format at devices from real 
#'   manufacturers differ significantly. Currently this function is checked 
#'   against RDML data from devices: \emph{Bio-Rad CFX96}, \emph{Roche 
#'   LightCycler 96} and \emph{Applied Biosystems StepOne}.
#' @param filename \code{string} -- path to RDML file
#' @param conditions.sep separator for condition defined at sample name
#' @author Konstantin A. Blagodatskikh <k.blag@@yandex.ru>, Stefan Roediger 
#'   <stefan.roediger@@hs-lausitz.de>, Michal Burdukiewicz 
#'   <michalburdukiewicz@@gmail.com>
#' @docType methods
#' @name new
#' @aliases RDML.new
#' @rdname new-method
#' @import XML
#' @include RDML.R
#' @examples
#' \dontrun{
#' ## Import from RDML file
#' PATH <- path.package("RDML")
#' filename <- paste(PATH, "/extdata/", "lc96_bACTXY.rdml", sep ="")
#' lc96 <- RDML$new(filename)
#' 
#' ## Some kind of overview for lc96
#' lc96$AsTable(name.pattern = sample[[react$sample$id]]$description)
#' lc96$AsDendrogram()
#' }
RDML$set("public", "initialize", function(filename,
                                          dateMade = NULL,
                                          dateUpdated = NULL,
                                          id = NULL,
                                          experimenter = NULL,
                                          documentation = NULL,
                                          dye = NULL,
                                          sample = NULL,
                                          target = NULL,
                                          thermalCyclingConditions = NULL,
                                          experiment = NULL,
                                          conditions.sep = NULL) {  
  if(missing(filename)) {
    assert_that(is.opt.string(dateMade))
    private$.dateMade <- dateMade
    assert_that(is.opt.string(dateUpdated))
    private$.dateUpdated <- dateUpdated
    assert_that(is.opt.list.type(id,
                                 rdmlIdType))
    private$.id <- id
    assert_that(is.opt.list.type(experimenter,
                                 experimenterType))
    private$.experimenter <- experimenter
    assert_that(is.opt.list.type(documentation,
                                 documentationType))                    
    private$.documentation <- documentation
    assert_that(is.opt.list.type(dye,
                                 dyeType))
    private$.dye <- dye
    assert_that(is.opt.list.type(sample,
                                 sampleType))
    private$.sample <- sample
    assert_that(is.opt.list.type(target,
                                 targetType))
    private$.target <- target
    assert_that(is.opt.list.type(thermalCyclingConditions,
                                 thermalCyclingConditionsType))
    private$.thermalCyclingConditions <- thermalCyclingConditions
    assert_that(is.opt.list.type(experiment,
                                 experimentType))
    private$.experiment <- experiment
    return()
  }
  assert_that(is.string(filename))
  # Unzips RDML to unique folder to get inner XML content.
  # Unique folder is needed to prevent file ovewriting
  # by parallel function usage.
  uniq.folder <- tempfile() #paste0(tempdir(), UUIDgenerate())
  cat(sprintf("Unzipping %s...", filename))
  unzipped.rdml <- unzip(filename, exdir = uniq.folder)
  dilutions.r <- NULL
  ref.genes.r <- NULL
  
  tryCatch({
    # Roche use more than one file at RDML zip.
    # One of the files store dilutions information.
    if(length(unzipped.rdml) > 1)
    {
      cat("\nParsing Roche(?) data...")
      rdml.doc <- xmlParse(paste0(uniq.folder,"/rdml_data.xml"))
      cat("OK")
      dilutions.r <- GetDilutionsRoche(uniq.folder)
      conditions.r <- GetConditionsRoche(uniq.folder)
      ref.genes.r <- GetRefGenesRoche(uniq.folder)
      # private$.dilutions <- dilutions.r
    }
    else
    {
      cat("\nParsing data...")
      rdml.doc <- xmlParse(unzipped.rdml)
      #     private$.dilutions <- GetDilutions(rdml.doc)
    }},
    error = function(e) { print(e) },
    finally = unlink(uniq.folder, recursive = TRUE)
  )
  ####
  
  rdml.root <- xmlRoot(rdml.doc)
  rdml.namespace <- c(rdml = "http://www.rdml.org")
  
  cat("\nGetting dateMade")
  private$.dateMade <- xmlValue(rdml.root[["dateMade"]])
  
  cat("\nGetting dateUpdated")
  private$.dateUpdated <- xmlValue(rdml.root[["dateUpdated"]])
  
  cat("\nGetting id")
  private$.id <- 
    llply(rdml.root["id"],
          function(id) {
            rdmlIdType$new(
              publisher = xmlValue(id[["publisher"]]),
              serialNumber = xmlValue(id[["serialNumber"]]),
              MD5Hash = xmlValue(id[["MD5Hash"]])
            )
          }
    ) %>% 
    with.names(quote(.$publisher))
  cat("\nGetting experementer")
  private$.experimenter <- {
    # experimenter.list <- 
    llply(rdml.root["experimenter"],
          function(experimenter)
            experimenterType$new(
              id = idType$new(xmlAttrs(experimenter, "id")),
              firstName = xmlValue(experimenter[["firstName"]]),
              lastName = xmlValue(experimenter[["lastName"]]),
              email = xmlValue(experimenter[["email"]]),
              labName = xmlValue(experimenter[["labName"]]),
              labAddress = xmlValue(experimenter[["labAddress"]])
            )
    ) %>% 
      with.names(quote(.$id$id))
  }
  
  cat("\nGetting documentation")
  private$.documentation <- {
    # documentation.list <- 
    llply(rdml.root["documentation"],
          function(el) 
            documentationType$new(
              id = xmlAttrs(el, "id"),
              text = xmlValue(el[["text"]])
            )) %>% 
      with.names(quote(.$id$id))
  }
  
  cat("\nGetting dye")
  private$.dye <-
    llply(rdml.root["dye"],
          function(el) {
            dyeType$new(
              id = idType$new(xmlAttrs(el, "id")),
              description = xmlValue(el[["description"]])
            )}) %>% 
    with.names(quote(.$id$id))
  
  
  cat("\nGetting sample")
  private$.sample <- 
    llply(rdml.root["sample"],
          function(sample) {
            type <- xmlValue(sample[["type"]])
            ######
            # remove Roche omitted ('ntp') samples
            if(type == "ntp")
              return(NULL)
            #####################
            # id <- xmlAttrs(sample, "id")
            sampleType$new(
              id = idType$new(xmlAttrs(sample, "id")),
              description = xmlValue(sample[["description"]]),
              documentation = 
                llply(sample["documentation"],
                      function(doc)
                        idReferencesType$new(xmlAttrs(doc, "id"))
                ),
              xRef = 
                llply(sample["xRef"],
                      function(xRef) 
                        xRefType$new(
                          name = xmlValue(xRef[["name"]]),
                          id = xmlValue(xRef[["id"]])
                        )),
              annotation = c(
                llply(sample["annotation"],
                      function(annotation)
                        annotationType$new(
                          property = xmlValue(annotation[["property"]]),
                          value = xmlValue(annotation[["value"]])
                        )),                  
                if (!is.null(conditions.sep))
                  annotationType$new(
                    property = "condition",
                    value = gsub(sprintf("^.*%s(.*)$",
                                         conditions.sep),
                                 "\\1", id))),
              type = sampleTypeType$new(type),
              interRunCalibrator = 
                as.logical(xmlValue(sample[["interRunCalibrator"]])),
              quantity = 
                quantityType$new(
                  value = as.numeric(
                    xmlValue(sample[["quantity"]][["value"]])),
                  unit = quantityUnitType$new(
                    sample[["quantity"]][["unit"]] %>% 
                    xmlValue)),
              calibratorSample = 
                as.logical(xmlValue(sample[["calibaratorSample"]])),
              cdnaSynthesisMethod = 
                cdnaSynthesisMethodType$new(
                  enzyme = xmlValue(sample[["cdnaSynthesisMethod"]][["enzyme"]]),
                  primingMethod =
                    primingMethodType$new(sample[["cdnaSynthesisMethod"]][["primingMethod"]],
                                          dnaseTreatment = as.logical(xmlValue(sample[["cdnaSynthesisMethod"]][["dnaseTreatment"]])),
                                          thermalCyclingConditions = 
                                            tryCatch(
                                              idReferencesType$new(
                                                xmlAttrs(sample[["cdnaSynthesisMethod"]][["thermalCyclingConditions"]],
                                                         "id")),
                                              error = function(e) NULL)
                    )),
              templateQuantity = 
                templateQuantityType$new(
                  conc = as.numeric(xmlValue(sample[["templateQuantity"]][["conc"]])),
                  nucleotide = nucleotideType$new(
                    xmlValue(sample[["templateQuantity"]][["nucleotide"]]))
                )
            )
            
          }) %>% 
    compact %>% 
    with.names(quote(.$id$id))
  
  cat("\nGetting target")
  private$.target <- 
    llply(rdml.root["target"],
          function(target) {
            targetType$new(
              id = idType$new({ 
                ifelse(length(private$.id) != 0 &&
                         private$.id[[1]]$publisher == "Roche Diagnostics",
                       {
                         id <- xmlAttrs(target, "id")
                         gsub("@(.+)$", "\\1", 
                              regmatches(id, gregexpr("@(.+)$",id))[[1]])
                       },
                       xmlAttrs(target, "id"))}),
              description = xmlValue(target[["description"]]),
              documentation = 
                llply(target["documentation"],
                      function(doc)
                        idReferencesType$new(xmlAttrs(doc, "id"))
                ),
              xRef = 
                llply(target["xRef"],
                      function(xRef) 
                        xRefType$new(
                          name = xmlValue(xRef[["name"]]),
                          id = xmlValue(xRef[["id"]])
                        )),
              type = targetTypeType$new(xmlValue(target[["type"]])),
              amplificationEfficiencyMethod = 
                xmlValue(target[["amplificationEfficiencyMethod"]]),
              amplificationEfficiency = 
                as.numeric(xmlValue(target[["amplificationEfficiency"]])),
              amplificationEfficiencySE = 
                as.numeric(xmlValue(target[["amplificationEfficiencySE"]])),
              detectionLimit = 
                as.numeric(xmlValue(target[["detectionLimit"]])),
              dyeId =
                tryCatch(
                  idReferencesType$new(xmlAttrs(target[["dyeId"]],"id")),
                  # StepOne stores dyeId as xmlValue 
                  error = function(e)
                    idReferencesType$new(xmlValue(target[["dyeId"]]))
                ),
              # dyeId = NA,
              
              sequences = sequencesType$new(
                forwardPrimer = 
                  oligoType$new(
                    threePrimeTag = 
                      xmlValue(target[["sequences"]][["forwardPrimer"]][["threePrimeTag"]]),
                    fivePrimeTag = 
                      xmlValue(target[["sequences"]][["forwardPrimer"]][["fivePrimeTag"]]),
                    sequence = 
                      xmlValue(target[["sequences"]][["forwardPrimer"]][["sequence"]])),
                reversePrimer = 
                  oligoType$new(
                    threePrimeTag = 
                      xmlValue(target[["sequences"]][["reversePrimer"]][["threePrimeTag"]]),
                    fivePrimeTag = 
                      xmlValue(target[["sequences"]][["reversePrimer"]][["fivePrimeTag"]]),
                    sequence = 
                      xmlValue(target[["sequences"]][["reversePrimer"]][["sequence"]])),
                probe1 = 
                  oligoType$new(
                    threePrimeTag = 
                      xmlValue(target[["sequences"]][["probe1"]][["threePrimeTag"]]),
                    fivePrimeTag = 
                      xmlValue(target[["sequences"]][["probe1"]][["fivePrimeTag"]]),
                    sequence = 
                      xmlValue(target[["sequences"]][["probe1"]][["sequence"]])),
                probe2 = 
                  oligoType$new(
                    threePrimeTag = 
                      xmlValue(target[["sequences"]][["probe2"]][["threePrimeTag"]]),
                    fivePrimeTag = 
                      xmlValue(target[["sequences"]][["probe2"]][["fivePrimeTag"]]),
                    sequence = 
                      xmlValue(target[["sequences"]][["probe2"]][["sequence"]])),
                amplicon = 
                  oligoType$new(
                    threePrimeTag = 
                      xmlValue(target[["sequences"]][["amplicon"]][["threePrimeTag"]]),
                    fivePrimeTag = 
                      xmlValue(target[["sequences"]][["amplicon"]][["fivePrimeTag"]]),
                    sequence = 
                      xmlValue(target[["sequences"]][["amplicon"]][["sequence"]]))
              ),
              commercialAssay = 
                commercialAssayType$new(
                  company = 
                    xmlValue(target[["commercialAssay"]][["company"]]),
                  orderNumber = 
                    xmlValue(target[["commercialAssay"]][["orderNumber"]]))
            )
          }
    ) %>% 
    with.names(quote(.$id$id))
  
  
  cat("\nGetting thermalCyclingConditions")
  private$.thermalCyclingConditions <- 
    llply(rdml.root["thermalCyclingConditions"],
          function(tcc) {
            thermalCyclingConditionsType$new(
              id = idType$new(xmlAttrs(tcc, "id")),
              description = xmlValue(tcc[["description"]]),
              documentation = 
                llply(tcc["documentation"],
                      function(doc)
                        idReferencesType$new(xmlAttrs(doc, "id"))
                ),
              lidTemperature = 
                as.numeric(xmlValue(tcc[["lidTemperature"]])),
              
              experimenter = llply(tcc["experimenter"],
                                   function(experimenter)
                                     idReferencesType$new(xmlAttrs(experimenter, "id"))
              ),
              
              step = llply(tcc["step"],
                           function(step) 
                             stepType$new(
                               nr = as.integer(xmlValue(step[["nr"]])),
                               description = xmlValue(step[["description"]]),
                               temperature = {
                                 if(is.null(step[["temperature"]][["temperature"]]))
                                   NULL
                                 else
                                   temperatureType$new(
                                     temperature = 
                                       as.numeric(xmlValue(step[["temperature"]][["temperature"]])),
                                     duration = 
                                       as.integer(xmlValue(step[["temperature"]][["duration"]])),
                                     temperatureChange = 
                                       as.numeric(xmlValue(step[["temperature"]][["temperatureChange"]])),
                                     durationChange = 
                                       as.integer(xmlValue(step[["temperature"]][["durationChange"]])),
                                     measure = measureType$new(
                                       xmlValue(step[["temperature"]][["measure"]])),
                                     ramp = 
                                       as.numeric(xmlValue(step[["temperature"]][["ramp"]]))
                                   )},
                               gradient = NULL,
                               #                                  gradientType$new(
                               #                                  highTemperature = 
                               #                                    as.numeric(xmlValue(step[["gradient"]][["highTemperature"]])),
                               #                                  lowTemperature = 
                               #                                    as.numeric(xmlValue(step[["gradient"]][["lowTemperature"]])),
                               #                                  duration = 
                               #                                    as.integer(xmlValue(step[["gradient"]][["duration"]])),
                               #                                  temperatureChange = 
                               #                                    as.numeric(xmlValue(step[["gradient"]][["temperatureChange"]])),
                               #                                  durationChange = 
                               #                                    as.integer(xmlValue(step[["gradient"]][["durationChange"]])),
                               #                                  measure = measureType$new(
                               #                                    xmlValue(step[["gradient"]][["measure"]])),
                               #                                  ramp = 
                               #                                    as.numeric(xmlValue(step[["gradient"]][["ramp"]]))
                               # ),
                               loop = NULL
                               #                                  {
                               #                                  if(is.null(xmlValue(step[["loop"]][["goto"]])))
                               #                                              NULL
                               #                                else
                               #                                              loopType$new(
                               #                                                goto = as.integer(xmlValue(step[["loop"]][["goto"]])),
                               #                                                # should be called "repeat" but this is reserved word
                               #                                                repeat.n = as.integer(xmlValue(step[["loop"]][["repeat"]])) 
                               #                                              )}
                               ,
                               pause = 
                               {
                                 if(is.null(xmlValue(step[["pause"]][["temperature"]])))
                                   NULL
                                 else
                                   pauseType$new(
                                     temperature = 
                                       as.numeric(xmlValue(step[["pause"]][["temperature"]]))
                                   )},
                               lidOpen = lidOpenType$new(xmlValue(step[["lidOpen"]][["lidOpenType"]]))
                             )
              )
            )
          }) %>% 
    with.names(quote(.$id$id))
  #     names(tcc.list) <- GetIds(tcc.list)
  #     tcc.list
  
  
  GetData <- function(data, experiment.id, run.id, react.id) {    
    tar.id <- xmlAttrs(data[["tar"]], "id")
    data.req <- paste0("/rdml:rdml/rdml:experiment[@id='",
                       experiment.id,
                       "']/rdml:run[@id='",                                                                         
                       run.id,
                       "']/rdml:react[@id='",
                       #                                                                          react.id[length(react.id)],
                       react.id,
                       "']/rdml:data/rdml:tar[@id='",
                       tar.id,
                       "']/..")                                                      
    dataType$new(
      tar = idReferencesType$new(
        ifelse(length(private$.id) != 0 &&
                 private$.id[[1]]$publisher == "Roche Diagnostics",
               gsub("@(.+)$", "\\1", 
                    regmatches(tar.id,gregexpr("@(.+)$",tar.id))[[1]])
               ,
               tar.id)), 
      cq = as.numeric(xmlValue(data[["cq"]])),
      excl = xmlValue(data[["excl"]]),
      adp = {                                                                                    
        cyc <- as.numeric(xpathSApply(rdml.doc,
                                      paste0(data.req,
                                             "/rdml:adp/rdml:cyc"),
                                      xmlValue,
                                      namespaces = c(rdml = "http://www.rdml.org")))
        tmp <- as.numeric(xpathSApply(rdml.doc,
                                      paste0(data.req,
                                             "/rdml:adp/rdml:tmp"),
                                      xmlValue,
                                      namespaces = c(rdml = "http://www.rdml.org")))                                                                        
        fluor <- as.numeric(xpathSApply(rdml.doc,
                                        paste0(data.req,
                                               "/rdml:adp/rdml:fluor"),
                                        xmlValue,
                                        namespaces = c(rdml = "http://www.rdml.org")))
        if(!is.null(fluor)) {
          if(length(tmp) != 0) {
            adpsType$new(matrix(c(cyc, tmp, fluor), 
                   byrow = FALSE,
                   ncol = 3,
                   dimnames = list(NULL,
                                   c("cyc", "tmp", "fluor"))))
          }
          else {
            adpsType$new(matrix(c(cyc, fluor), 
                   byrow = FALSE,
                   ncol = 2,
                   dimnames = list(NULL,
                                   c("cyc", "fluor"))))
          }
        } else {
          #           matrix(ncol = 2,
          #                  dimnames = list(NULL,
          #                                  c("cyc", "tmp", "fluor")))
          NULL
        }
      },
      mdp = {                                                             
        tmp <- as.numeric(xpathSApply(rdml.doc,
                                      paste0(data.req,
                                             "/rdml:mdp/rdml:tmp"),
                                      xmlValue,
                                      namespaces = c(rdml = "http://www.rdml.org")))
        fluor <- as.numeric(xpathSApply(rdml.doc,
                                        paste0(data.req,
                                               "/rdml:mdp/rdml:fluor"),
                                        xmlValue,
                                        namespaces = c(rdml = "http://www.rdml.org")))
        
        if(length(fluor) != 0 && !is.null(fluor)) {
#           matrix(c(tmp, fluor), 
#                                                byrow = FALSE,
#                                                ncol = 2,
#                                                dimnames = list(NULL,
#                                                                c("tmp", "fluor"))) %>% 
#             typeof %>% print
#           NULL
          mdpsType$new(matrix(c(tmp, fluor), 
                              byrow = FALSE,
                              ncol = 2,
                              dimnames = list(NULL,
                                              c("tmp", "fluor"))))
        }
        else
          #           matrix(ncol = 2,
          #                  dimnames = list(NULL,
          #                                  c("tmp", "fluor")))
          NULL
      },
      endPt = as.numeric(xmlValue(data[["endPt"]])),
      bgFluor = as.numeric(xmlValue(data[["bgFluor"]])),
      bgFluorSlp = as.numeric(xmlValue(data[["bgFluorSp"]])),
      quantFluor = as.numeric(xmlValue(data[["quantFluor"]]))
    )
  }
  
  GetReact <- function(react, experiment.id, run.id) {
    react.id <- xmlAttrs(react, "id")    
    react.id.corrected <- tryCatch(
      as.integer(react.id),
      warning = function(cond) {
        # if react.id is 'B1' not '13'
        # like in StepOne
        FromPositionToId(react.id)
      }    
    )
    #     cat(sprintf("\nreact: %i", react.id))
    sample <- xmlAttrs(react[["sample"]],"id")
    ######
    if(length(private$.id) != 0 && 
       private$.id[[1]]$publisher == "Roche Diagnostics") {
      # remove Roche omitted ('ntp') samples
      if(is.null(private$.sample[[sample]]))
        return(NULL)
      ## Better names for Roche    
      sample <- private$.sample[[xmlAttrs(react[["sample"]],"id")]]$description
    }
    #######    
    reactType$new(
      id = reactIdType$new(react.id.corrected), #sample.id
      #       # will be calculated at the end of init
      #       position = NA,
      sample = idReferencesType$new(sample),
      data = {
        llply(react["data"],
              function(data) GetData(data,
                                     experiment.id,
                                     run.id,
                                     react.id)
        ) 
      }
    )
  }  
  
  GetRun <- function(run, experiment.id) {
    run.id <-xmlAttrs(run, "id")
    cat(sprintf("\nrun: %s\n", run.id))
    runType$new(
      id = idType$new(run.id), #xmlAttrs(run, "id"),
      description = xmlValue(run[["description"]]),
      documentation = 
        llply(run["documentation"],
              function(doc)
                idReferencesType$new(xmlAttrs(doc, "id"))),
      experimenter = 
        llply(run["experimenter"],
              function(doc)
                idReferencesType$new(xmlAttrs(doc, "id"))),
      instrument = xmlValue(run[["instrument"]]),
      dataCollectionSoftware = dataCollectionSoftwareType$new(
        name = xmlValue(run[["dataCollectionSoftware"]][["name"]]),
        version = xmlValue(run[["dataCollectionSoftware"]][["version"]])
      ),
      backgroundDeterminationMethod = 
        xmlValue(run[["backgroundDeterminationMethod"]]),
      cqDetectionMethod = 
        cqDetectionMethodType$new(xmlValue(run[["cqDetectionMethod"]])),
      thermalCyclingConditions = 
        tryCatch(
          idReferencesType$new(
            xmlAttrs(run[["thermalCyclingConditions"]], "id")),
          error = function(e) NULL),
      pcrFormat = 
      {
        rows <- as.integer(xmlValue(run[["pcrFormat"]][["rows"]]))
        # check for absent of pcrFormat
        # like in StepOne
        if(!is.na(rows)) {
          pcrFormatType$new(
            rows = rows,
            columns = as.integer(xmlValue(run[["pcrFormat"]][["columns"]])),
            rowLabel = labelFormatType$new(
              xmlValue(run[["pcrFormat"]][["rowLabel"]])),
            columnLabel = labelFormatType$new(
              xmlValue(run[["pcrFormat"]][["columnLabel"]]))
          )
        } else {
          pcrFormatType$new(
            rows = 12,
            columns = 8,
            rowLabel = labelFormatType$new("ABC"),
            columnLabel = labelFormatType$new("123")
          )
        }
      },
      runDate = xmlValue(run[["runDate"]]),
      react =
        llply(run["react"],
              function(react) GetReact(react, 
                                       experiment.id,
                                       run.id),
              .parallel = FALSE,
              .progress = ifelse(interactive(),
                                 "text",
                                 "none")
        ) %>% 
        compact 
    )
    
  }                      
  
  GetExperiment <- function(experiment) {
    experiment.id <- xmlAttrs(experiment, "id")
    cat(sprintf("\nGetting experiment: %s", experiment.id))
    experimentType$new(
      id = idType$new(experiment.id),
      description = xmlValue(experiment[["description"]]),
      documentation = 
        llply(experiment["documentation"],
              function(doc)
                idReferencesType$new(xmlAttrs(doc, "id"))),
      run = 
        llply(experiment["run"],
              function(run) GetRun(run, experiment.id)
        )
    )
    
  }
  
  
  private$.experiment <- 
    llply(rdml.root["experiment"],
          function(experiment) GetExperiment(experiment)
    ) %>% 
    with.names(quote(.$id$id))
  
  # return()
  # private$.recalcPositions()
  
  if (length(private$.id) != 0 && private$.id[[1]]$publisher == "Roche Diagnostics") {    
    for(i in 1:length(private$.sample)) {
      private$.sample[[i]]$id <- idType$new(private$.sample[[i]]$description)
    }
    private$.sample <- with.names(private$.sample,
                                  quote(.$id$id))
    
    cat("Adding Roche ref genes\n")
    if(!is.null(ref.genes.r) && length(ref.genes.r) != 0) {
      for(ref.gene in ref.genes.r) {
        geneName <- xmlValue(ref.gene[["geneName"]])
        geneI <- grep(
          sprintf("^%s$", geneName),
          names(private$.target))
        private$.target[[geneI]]$type <-
          targetTypeType$new(
            ifelse(as.logical(xmlValue(ref.gene[["isReference"]])),
                   "ref",
                   "toi"))
      }
    }
    # return()
    tbl <- self$AsTable()
    cat("Adding Roche quantities\n")
    for(target in dilutions.r %>% names) {
      for(r.id in dilutions.r[[target]] %>% names) {
        sample.name <- filter(tbl, react.id == r.id)$sample[1]
        private$.sample[[sample.name]]$quantity <- 
          quantityType$new(
          value = unname(dilutions.r[[1]][r.id]),
          unit = quantityUnitType$new("other")
        )
        private$.sample[[sample.name]]$annotation <- 
          c(private$.sample[[sample.name]]$annotation,
          annotationType$new(
                  property = sprintf("Roche_quantity_at_%s_%s",
                                     target,
                                     r.id),
                  value = as.character(dilutions.r[[target]][r.id])))
      }
    }
    
    cat("Adding Roche conditions\n")
    for(r.id in conditions.r %>% names) {
      sample.name <- filter(tbl, react.id == r.id)$sample[1]
      private$.sample[[sample.name]]$annotation <- 
        c(private$.sample[[sample.name]]$annotation,
          annotationType$new(
                property = sprintf("Roche_condition_at_%s",r.id),
                value = conditions.r[r.id]))
    }
    
  }
}, 
overwrite = TRUE)