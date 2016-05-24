with.names <- function(l, id) {
  if (is.null(l))
    return(NULL)
  n <- list.select(l,
                   eval(id)) %>% 
    unlist
  if (is.null(n))
    return(l)
  names(l) <- n
  l
}

# rdmlBaseType ------------------------------------------------------------

#' Base R6 class for RDML package.
#' 
#' Most classes of RDML package inherit this class. It can't be directly 
#' accessed and serves only for inner package usage.
#' 
#' @section Initialization: \code{rdmlBaseType$new()}
#'   
#' @section Methods: \describe{\item{\code{.asXMLnodes(node.name)}}{Represents
#'   object as XML nodes. Should not be called directly. \code{node.name} --
#'   name of the root node for the generated XML
#'   tree.}\item{\code{print(...)}}{prints object}}
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
rdmlBaseType <-
  R6Class("rdmlBaseType",
          # class = FALSE,
          public = list(
            .asXMLnodes = function(node.name,
                                   namespaceDefinitions = NULL) {
              # has.attrs <- FALSE
              subnodes <- names(private)
              newXMLNode(
                name = node.name,
                namespaceDefinitions = namespaceDefinitions,
                attrs = {
                  get.attrs <- function() {
                    for(name in subnodes) {
                      if(class(private[[name]])[1] == "idType" ||
                         class(private[[name]])[1] == "reactIdType") {
                        subnodes <<- subnodes[subnodes != name]
                        # has.attrs <<- TRUE
                        return(list(id = private[[name]]$id))
                      }
                    }
                    NULL
                  }
                  get.attrs()
                },
                .children = {
                  llply(subnodes,
                        function(name) {
                          subnode.name <- gsub("^\\.(.*)$",
                                               "\\1", name)
                          switch(
                            typeof(private[[name]]),
                            closure = NULL,
                            list = 
                              llply(private[[name]],
                                    function(sublist)
                                      sublist$.asXMLnodes(subnode.name))
                            ,
                            environment = {
                              switch(class(private[[name]])[2],
                                     enumType = {
                                       if (is.null(private[[name]]$value) ||
                                           is.na(private[[name]]$value))
                                         NULL
                                       else
                                         newXMLNode(name = subnode.name,
                                                    text = private[[name]]$value)
                                     },
                                     idType = 
                                       newXMLNode(name = subnode.name,
                                                  attrs = 
                                                    list(id = private[[name]]$id)),
                                     private[[name]]$.asXMLnodes(subnode.name)
                              )},
                            {
                              if (is.null(private[[name]]) ||
                                  is.na(private[[name]]))
                                NULL
                              else
                                newXMLNode(name = subnode.name,
                                           text = private[[name]])
                            })
                        }) %>% 
                    compact
                }
              )
            },
            print = function(...) {
              sapply(names(private)[-which(names(private) == 
                                             "deep_clone")],
                     function(name) {
                       sprintf(
                         "\t%s: %s",
                         gsub("^\\.(.*)$",
                              "\\1", name),
                         switch(
                           typeof(private[[name]]),
                           closure = NULL,
                           list = sprintf("[%s]",
                                          names(private[[name]]) %>% 
                                            paste(collapse = ", ")),
                           environment = {
                             sprintf("~ %s",
                                     class(private[[name]])[1])
                           },
                           NULL = "",
                           {
                             if (class(private[[name]]) == "matrix")
                               sprintf("%s fluorescence data points",
                                       nrow(private[[name]]))
                             else
                               sprintf("%s", private[[name]])
                           }))
                     }) %>% 
                paste(sep = "\n", collapse = "\n") %>% 
                cat
            }
          ),
          private = list(
            deep_clone = function(name, value) {
              if (value %>% is.null)
                return(NULL)
              if (value %>%
                  class %>% 
                  tail(1) != "R6") {
                if (is.list(value)) {
                  llply(value,
                        function(el) el$clone(deep = TRUE))
                } else {
                  value
                }
              } else {
                value$clone(deep = TRUE)
              }
            }
          )
  )

# rdmlIdType ------------------------------------------------------------

#' rdmlIdType R6 class.
#' 
#' This element can be used to assign a publisher and id to the RDML file.\cr 
#' Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{rdmlIdType$new(publisher, serialNumber,
#'   MD5Hash = NULL)}
#'   
#' @section Fields: \describe{  
#'   \item{\code{publisher}}{\link[assertthat]{is.string}. RDML file publisher.}
#'   \item{\code{serialNumber}}{\link[assertthat]{is.string}. Serial number.} 
#'   \item{\code{MD5Hash}}{\link[assertthat]{is.string}. An MD5Hash calculated
#'   over the complete file after removing all rdmlIDTypes and all whitespaces
#'   between elements.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
rdmlIdType <- 
  R6Class("rdmlIdType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(publisher,
                                  serialNumber,
                                  MD5Hash = NULL) {
              assert_that(is.string(publisher))
              assert_that(is.string(serialNumber))
              assert_that(is.opt.string(MD5Hash))
              private$.publisher <- publisher
              private$.serialNumber <- serialNumber
              private$.MD5Hash <- MD5Hash
            }
          ),
          private = list(
            .publisher = NULL,
            .serialNumber = NULL,
            .MD5Hash = NULL
          ),
          active = list(
            publisher = function(publisher) {
              if (missing(publisher))
                return(private$.publisher)
              assert_that(is.string(publisher))
              private$.publisher <- publisher
            },
            serialNumber = function(serialNumber) {
              if (missing(serialNumber))
                return(private$.serialNumber)
              assert_that(is.string(serialNumber))
              private$.serialNumber <- serialNumber
            },
            MD5Hash = function(MD5Hash) {
              if (missing(MD5Hash))
                return(private$.MD5Hash)
              assert_that(is.opt.string(MD5Hash))
              private$.MD5Hash <- MD5Hash
            }
          ))

# idType ------------------------------------------------------------

#' idType R6 class.
#' 
#' Contains identificator for varius RDML types.Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{idType$new(id)}
#'
#'   @section Fields: \describe{     
#' \item{\code{id}}{\link[assertthat]{is.string}. Identificator.}
#' }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
idType <- 
  R6Class("idType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(id) {
              assert_that(is.string(id))
              private$.id <- id
            }#,
            #             print = function(...) {
            #               cat(private$.id)
            #             }
          ),
          private = list(
            .id = NULL
          ),
          active = list(
            id = function(id) {
              if (missing(id))
                return(private$.id)
              assert_that(is.string(id))
              private$.id <- id
            }
          ))

# reactIdType ------------------------------------------------------------

#' reactIdType R6 class.
#' 
#' Contains identificator for reactType.Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{reactIdType$new(id)}
#'
#'   @section Fields: \describe{     
#' \item{\code{id}}{\link[assertthat]{is.count}. Identificator.}
#' }
#' 
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
reactIdType <- 
  R6Class("reactIdType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(id) {
              assert_that(is.count(id))
              private$.id <- id
            }#,
            #             print = function(...) {
            #               cat(private$.id)
            #             }
          ),
          private = list(
            .id = NULL
          ),
          active = list(
            id = function(id) {
              if (missing(id))
                return(private$.id)
              assert_that(is.count(id))
              private$.id <- id
            }
          ))

# idReferencesType ------------------------------------------------------------

#' idReferencesType R6 class.
#' 
#' Contains id of another RDML object.Inherits: \link{idType}.
#' 
#' @section Initialization: \code{idReferencesType$new(id)}
#' 
#' @section Fields: \describe{  
#' \item{\code{id}}{\link[assertthat]{is.string}. Identificator.}
#' }
#'         
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
idReferencesType <- 
  R6Class("idReferencesType",
          # class = FALSE,
          inherit = idType)

# experimenterType ------------------------------------------------------------

#' experimenterType R6 class.
#' 
#' Contact details of the experimenter.Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{experimenterType$new(id, firstName, lastName,
#'   email = NULL, labName = NULL, labAddress = NULL)}
#'  
#'   @section Fields: \describe{   
#' \item{\code{id}}{\link{idType}. Identificator.}
#' \item{\code{firstName}}{\link[assertthat]{is.string}. First name.}
#' \item{\code{lastName}}{\link[assertthat]{is.string}. Last name.}
#' \item{\code{email}}{\link[assertthat]{is.string}. Email.}
#' \item{\code{labName}}{\link[assertthat]{is.string}. Lab name.}
#' \item{\code{labAddress}}{\link[assertthat]{is.string}. Lab address.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
experimenterType <- 
  R6Class("experimenterType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(id,
                                  firstName,
                                  lastName,
                                  email = NULL,
                                  labName = NULL,
                                  labAddress = NULL) {
              assert_that(is.type(id, idType))
              assert_that(is.string(firstName))
              assert_that(is.string(lastName))
              assert_that(is.opt.string(email))
              assert_that(is.opt.string(labName))
              assert_that(is.opt.string(labAddress))
              private$.id <- id
              private$.firstName <- firstName
              private$.lastName <- lastName
              private$.email <- email
              private$.labName <- labName
              private$.labAddress <- labAddress
            }
          ),
          private = list(
            .id = NULL,
            .firstName = NULL,
            .lastName = NULL,
            .email = NULL,
            .labName = NULL,
            .labAddress = NULL
          ),
          active = list(
            id = function(id) {
              if (missing(id))
                return(private$.id)
              assert_that(is.type(id, idType))
              private$.id <- id
            },
            firstName = function(firstName) {
              if (missing(firstName))
                return(private$.firstName)
              assert_that(is.string(firstName))
              private$.firstName <- firstName
            },
            lastName = function(lastName) {
              if (missing(lastName))
                return(private$.lastName)
              assert_that(is.string(lastName))
              private$.lastName <- lastName
            },
            email = function(email) {
              if (missing(email))
                return(private$.email)
              assert_that(is.opt.string(email))
              private$.email <- email
            },
            labName = function(labName) {
              if (missing(labName))
                return(private$.labName)
              assert_that(is.opt.string(labName))
              private$.labName <- labName
            },
            labAddress = function(labAddress) {
              if (missing(labAddress))
                return(private$.labAddress)
              assert_that(is.opt.string(labAddress))
              private$.labAddress <- labAddress
            }
          ))

# documentationType ------------------------------------------------------------

#' documentationType R6 class.
#' 
#' These elements should be used if the same description applies to many
#' samples, targets or experiments.Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{documentationType$new(id, text = NULL)}
#'   
#'   @section Fields: \describe{  
#' \item{\code{id}}{\link{idType}. Identificator.}
#' \item{\code{text}}{\link[assertthat]{is.string}. Text.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
documentationType <- 
  R6Class("documentationType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(id,
                                  text = NULL) {
              assert_that(is.type(id, idType))
              assert_that(is.opt.string(text))
              private$.id <- id
              private$.text <- text
            }
          ),
          private = list(
            .id = NULL,
            .text = NULL
          ),
          active = list(
            id = function(id) {
              if (missing(id))
                return(private$.id)
              assert_that(is.type(id, idType))
              private$.id <- id
            },
            text = function(text) {
              if (missing(text))
                return(private$.text)
              assert_that(is.opt.string(text))
              private$.text <- text
            }
          ))

# dyeType ------------------------------------------------------------

#' dyeType R6 class.
#' 
#' Information on a dye.Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{dyeType$new(id, description = NULL)}
#'   
#'   @section Fields: \describe{  
#' \item{\code{id}}{\link{idType}. Identificator.}
#' \item{\code{description}}{ \link[assertthat]{is.string}. Description.
#'   }}
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
dyeType <- 
  R6Class("dyeType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(id,
                                  description = NULL) {
              assert_that(is.type(id, idType))
              assert_that(is.opt.string(description))
              private$.id <- id
              private$.description <- description
            }
          ),
          private = list(
            .id = NULL,
            .description = NULL
          ),
          active = list(
            id = function(id) {
              if (missing(id))
                return(private$.id)
              assert_that(is.type(id, idType))
              private$.id <- id
            },
            description = function(description) {
              if (missing(description))
                return(private$.description)
              assert_that(is.opt.string(description))
              private$.description <- description
            }
          ))

# xRefType ------------------------------------------------------------

#' xRefType R6 class.
#' 
#' Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{xRefType$new(name = NULL, id = NULL)}
#'   
#'   @section Fields: \describe{  
#'   \item{\code{name}}{\link[assertthat]{is.string}. Reference to an external
#'   database, } for example "GenBank". 
#'   \item{\code{id}}{\link[assertthat]{is.string}. The ID of the entry within
#'   the external database, for example "AJ832138".}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
xRefType <- 
  R6Class("xRefType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(name = NULL,
                                  id = NULL) {
              assert_that(is.opt.string(name))
              assert_that(is.opt.string(id))
              private$.name <- name
              private$.id <- id
            }
          ),
          private = list(
            .name = NULL,
            .id = NULL
          ),
          active = list(
            name = function(name) {
              if (missing(name))
                return(private$.name)
              assert_that(is.opt.string(name))
              private$.name <- name
            },
            id = function(id) {
              if (missing(id))
                return(private$.id)
              assert_that(is.opt.string(id))
              private$.id <- id
            }
          ))

# annotationType ------------------------------------------------------------

#' annotationType R6 class.
#' 
#' These elements should be used to annotate samples by setting a property and a
#' value. A property could be sex, the value M or F.Inherits:
#' \link{rdmlBaseType}.
#' 
#' @section Fields: \describe{ \item{property}{\link[assertthat]{is.string}.
#'   Property} \item{value}{\link[assertthat]{is.string}. Value} }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
annotationType <- 
  R6Class("annotationType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(property,
                                  value) {
              assert_that(is.string(property))
              assert_that(is.string(value))
              private$.property <- property
              private$.value <- value
            }
          ),
          private = list(
            .property = NULL,
            .value = NULL
          ),
          active = list(
            property = function(property) {
              if (missing(property))
                return(private$.property)
              assert_that(is.string(property))
              private$.property <- property
            },
            value = function(value) {
              if (missing(value))
                return(private$.value)
              assert_that(is.string(value))
              private$.value <- value
            }
          ))

# quantityType ------------------------------------------------------------

#' quantityType R6 class.
#' 
#' A quantity is always defined by its value and its unit.Inherits:
#' \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{quantityType$new(value, unit)}
#' 
#'   @section Fields: \describe{    
#' \item{\code{value}}{\link[base]{is.double}. Value.}
#' \item{\code{unit}}{\link{quantityUnitType}. Unit.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
quantityType <- 
  R6Class("quantityType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(value,
                                  unit) {
              assert_that(is.float(value))
              assert_that(is.type(unit,
                                  quantityUnitType))
              private$.value <- value
              private$.unit <- unit
            }
          ),
          private = list(
            .value = NULL,
            .unit = NULL
          ),
          active = list(
            value = function(value) {
              if (missing(value))
                return(private$.value)
              assert_that(is.float(value))
              private$.value <- value
            },
            unit = function(unit) {
              if (missing(unit))
                return(private$.unit)
              assert_that(is.type(unit,
                                  quantityUnitType))
              private$.unit <- unit
            }
          ))

# cdnaSynthesisMethodType ------------------------------------------------------------

#' cdnaSynthesisMethodType R6 class.
#' 
#' Description of the cDNA synthesis method.Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{cdnaSynthesisMethodType$new(enzyme = NULL, 
#'   primingMethod = NULL, dnaseTreatment = NULL, thermalCyclingConditions = 
#'   NULL)}
#'   
#'   @section Fields: \describe{  
#'   \item{\code{enzyme}}{\link[base]{is.double}. Enzyme used for reverse
#'   transcription.} \item{\code{primingMethod}}{\link{primingMethodType}.} 
#'   \item{\code{dnaseTreatment}}{\link[assertthat]{is.flag} True if RNA was
#'   DNAse treated prior cDNA synthesis.}
#'   \item{\code{thermalCyclingConditions}}{\link{idReferencesType}.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
cdnaSynthesisMethodType <- 
  R6Class("cdnaSynthesisMethodType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(enzyme = NULL,
                                  primingMethod = NULL,
                                  dnaseTreatment = NULL,
                                  thermalCyclingConditions = NULL) {
              assert_that(is.opt.string(enzyme))
              assert_that(is.opt.type(primingMethod,
                                      primingMethodType))
              assert_that(is.opt.flag(dnaseTreatment))
              assert_that(is.opt.type(thermalCyclingConditions,
                                      idType))
              private$.enzyme <- enzyme
              private$.primingMethod <- primingMethod
              private$.dnaseTreatment <- dnaseTreatment
              private$.thermalCyclingConditions <- thermalCyclingConditions
            }
          ),
          private = list(
            .enzyme = NULL,
            .primingMethod = NULL,
            .dnaseTreatment = NULL,
            .thermalCyclingConditions = NULL
          ),
          active = list(
            enzyme = function(enzyme) {
              if (missing(enzyme))
                return(private$.enzyme)
              assert_that(is.opt.string(enzyme))
              private$.enzyme <- enzyme
            },
            primingMethod = function(primingMethod) {
              if (missing(primingMethod))
                return(private$.primingMethod)
              assert_that(is.opt.type(primingMethod,
                                      primingMethodType))
              private$.primingMethod <- primingMethod
            },
            dnaseTreatment = function(dnaseTreatment) {
              if (missing(dnaseTreatment))
                return(private$.dnaseTreatment)
              assert_that(is.opt.flag(dnaseTreatment))
              private$.dnaseTreatment <- dnaseTreatment
            },
            thermalCyclingConditions = function(thermalCyclingConditions) {
              if (missing(thermalCyclingConditions))
                return(private$.thermalCyclingConditions)
              assert_that(is.opt.type(thermalCyclingConditions,
                                      idType))
              private$.thermalCyclingConditions <- thermalCyclingConditions
            }
          ))

# templateQuantityType ------------------------------------------------------------

#' templateQuantityType R6 class.
#' 
#' Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{templateQuantityType$new(conc, nucleotide)}
#'   
#'   @section Fields: \describe{  
#' \item{\code{conc}}{\link[base]{is.double}. Concentration of the template in nanogram}
#'   per microliter in the final reaction mix.
#' \item{\code{nucleotide}}{\link{nucleotideType}.}
#'   }
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
templateQuantityType <- 
  R6Class("templateQuantityType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(conc,
                                  nucleotide) {
              assert_that(is.float(conc))
              assert_that(is.type(nucleotide,
                                  nucleotideType))
              private$.conc <- conc
              private$.nucleotide <- nucleotide
            }
          ),
          private = list(
            .conc = NULL,
            .nucleotide = NULL
          ),
          active = list(
            conc = function(conc) {
              if (missing(conc))
                return(private$.conc)
              assert_that(is.float(conc))
              private$.conc <- conc
            },
            nucleotide = function(nucleotide) {
              if (missing(nucleotide))
                return(private$.nucleotide)
              assert_that(is.type(nucleotide,
                                  nucleotideType))
              private$.nucleotide <- nucleotide
            }
          ))

# enumType ------------------------------------------------------------

#' enumType R6 class.
#' 
#' Generic class for creating objects thet can take limited list of values.\cr
#' Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{enumType$new(value)}
#'   @section Fields: \describe{  
#' \item{\code{value}}{\link[assertthat]{is.string}. Value.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
enumType <- 
  R6Class("enumType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(value, ...) {
              assert_that(is.enum(value, private$.levels))
              private$.value <- value
            },
            print = function(...) {
              cat(private$.value)
            }
          ),
          private = list(
            .value = NULL
          ),
          active = list(
            value = function(value) {
              if (missing(value))
                return(private$.value)
              assert_that(is.enum(value, private$.levels))
              private$.value <- value
            }
          ))

# sampleTypeType ------------------------------------------------------------

#' sampleTypeType R6 class.
#' 
#' Can take values:
#' \describe{
#' \item{unkn}{unknown sample}
#' \item{ntc}{non template control}
#' \item{nac}{no amplification control}
#' \item{std}{standard sample}
#' \item{ntp}{no target present}
#' \item{nrt}{minusRT}
#' \item{pos}{positive control}
#' \item{opt}{optical calibrator sample}}
#' 
#' Inherits: \link{enumType}.
#' 
#' @section Initialization: \code{sampleTypeType$new(value)}
#'  
#'   @section Fields: \describe{   
#' \item{\code{value}}{\link[assertthat]{is.string}. Value.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
sampleTypeType <- 
  R6Class("sampleTypeType",
          # class = FALSE,
          inherit = enumType,
          private = list(
            .levels = c("unkn", "ntc", "nac",
                        "std", "ntp", "nrt",
                        "pos", "opt")
          )
  )

# quantityUnitType ------------------------------------------------------------

#' quantityUnitType R6 class.
#' 
#' The unit the quantity. Can take values:
#' \describe{
#' \item{cop}{copies per microliter  }
#' \item{fold}{fold change  }
#' \item{dil}{dilution (10 would mean 1:10 dilution)  }
#' \item{nMol}{nanomol per microliter }
#' \item{ng}{nanogram per microliter }
#' \item{other}{other unit (must be linear, no exponents or logarithms allowed) }
#' }
#' 
#' Inherits: \link{enumType}.
#' 
#' @section Initialization: \code{quantityUnitType$new(value)}
#'   
#'   @section Fields: \describe{  
#' \item{\code{value}}{\link[assertthat]{is.string}. Value.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
quantityUnitType <- 
  R6Class("quantityUnitType",
          # class = FALSE,
          inherit = enumType,
          private = list(
            .levels = c("cop", "fold", "dil", "ng", "nMol", "other")
          )
  )

# primingMethodType ------------------------------------------------------------

#' primingMethodType R6 class.
#' 
#' The primers used to reverse transcribe the RNA to cDNA. Can take values:
#' \describe{
#' \item{oligo-dt}{}
#' \item{random}{}
#' \item{target-specific}{}
#' \item{oligo-dt and random}{}
#' \item{other}{} 
#' }
#' 
#' Inherits: \link{enumType}.
#' 
#' @section Initialization: \code{primingMethodType$new(value)}
#'   
#'   @section Fields: \describe{  
#' \item{\code{value}}{\link[assertthat]{is.string}. Value.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
primingMethodType <- 
  R6Class("primingMethodType",
          # class = FALSE,
          inherit = enumType,
          private = list(
            .levels = c("oligo-dt",
                        "random",
                        "target-specific",
                        "oligo-dt and random",
                        "other")
          )
  )

# nucleotideType ------------------------------------------------------------

#' nucleotideType R6 class.
#' 
#' Can take values:
#' \describe{
#' \item{DNA}{}
#' \item{genomic-DNA}{}
#' \item{cDNA}{}
#' \item{RNA}{}
#' }
#' 
#' Inherits: \link{enumType}.
#' 
#' @section Initialization: \code{nucleotideType$new(value)}
#'   
#'   @section Fields: \describe{  
#' \item{\code{value}}{\link[assertthat]{is.string}. Value.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
nucleotideType <- 
  R6Class("nucleotideType",
          # class = FALSE,
          inherit = enumType,
          private = list(
            .levels = c("DNA",
                        "genomic-DNA",
                        "cDNA",
                        "RNA")
          )
  )

# sampleType ------------------------------------------------------------

#'sampleType R6 class.
#'
#'A sample is a defined template solution. Dilutions of the same material differ
#'in concentration and are considered different samples. A technical replicate 
#'samples should contain the same name (reactions are performed on the same 
#'material), and biological replicates should contain different names (the 
#'nucleic acids derived from the different biological replicates are not the 
#'same). Serial dilutions in a standard curve must have a different name.\cr 
#'Inherits: \link{rdmlBaseType}.
#'
#'@section Initialization: \code{sampleType$new(id, description = NULL, 
#'  documentation = NULL, xRef =  NULL, annotation = NULL, type = 
#'  sampleTypeType$new("unkn"), interRunCalibrator = FALSE, quantity = NULL, 
#'  calibratorSample = FALSE, cdnaSynthesisMethod = NULL, templateQuantity = 
#'  NULL)}
#'  
#'  @section Fields: \describe{  
#'  \item{\code{id}}{\link{idType}. Concentration of the template in nanogram
#'  per microliter in the final reaction mix. }
#'  \item{\code{description}}{\link[assertthat]{is.string}.} 
#'  \item{\code{documentation}}{\code{list} of \link{idReferencesType}.} 
#'  \item{\code{xRef}}{\code{list} of \link{xRefType}.} 
#'  \item{\code{annotation}}{\code{list} of \link{annotationType}.} 
#'  \item{\code{type}}{\link{sampleTypeType}.} 
#'  \item{\code{interRunCalibrator}}{\link[assertthat]{is.flag}. True if this
#'  sample is used as inter run calibrator. }
#'  \item{\code{quantity}}{\link{quantityType}. Quantity - The reference
#'  quantity of this sample. It should be only used if the sample is part of a
#'  standard curve. The provided value will be used to quantify unknown samples
#'  in absolute quantification assays. Only the use of true numbers is valid
#'  like 1, 10, 100, 1000 or 1, 0.1, 0.01, 0.001. The use of exponents is not
#'  valid like 1, 2, 3, 4 or -1, -2, -3, -4 because it will not be interpreted
#'  as 10E1, 10E2, 10E3, 10E4 or 10E-1, 10E-2, 10E-3, 10E-4. }
#'  \item{\code{calibratorSample}}{\link[assertthat]{is.flag}. True if this
#'  sample is used as calibrator sample. }
#'  \item{\code{cdnaSynthesisMethod}}{\link{cdnaSynthesisMethodType}.} 
#'  \item{\code{templateQuantity}}{\link{templateQuantityType}.}
#'  }
#'  
#'@docType class
#'@format An \code{\link{R6Class}} generator object.
#'@export
sampleType <- 
  R6Class("sampleType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(id,
                                  description = NULL,
                                  documentation = NULL,
                                  xRef = NULL,
                                  annotation = NULL,
                                  type = sampleTypeType$new("unkn"),
                                  interRunCalibrator = FALSE,
                                  quantity = NULL,
                                  calibratorSample = FALSE,
                                  cdnaSynthesisMethod = NULL,
                                  templateQuantity = NULL) {
              assert_that(is.type(id, idType))
              assert_that(is.opt.string(description))
              assert_that(is.opt.list.type(documentation,
                                           idReferencesType))
              assert_that(is.opt.list.type(xRef,
                                           xRefType))
              assert_that(is.opt.list.type(annotation,
                                           annotationType))
              assert_that(is.type(type,
                                  sampleTypeType))
              assert_that(is.opt.flag(interRunCalibrator))
              assert_that(is.opt.type(quantity,
                                      quantityType))
              assert_that(is.opt.flag(calibratorSample))
              assert_that(is.opt.type(cdnaSynthesisMethod,
                                      cdnaSynthesisMethodType))
              assert_that(is.opt.type(templateQuantity,
                                      templateQuantityType))
              
              private$.id <- id
              private$.description <- description
              private$.documentation <- documentation
              private$.xRef <- xRef
              private$.annotation <- with.names(annotation,
                                                quote(.$property))
              private$.type <- type
              private$.interRunCalibrator <- interRunCalibrator
              private$.quantity <- quantity
              private$.calibratorSample <- calibratorSample
              private$.cdnaSynthesisMethod <- cdnaSynthesisMethod
              private$.templateQuantity <- templateQuantity
            }
          ),
          private = list(
            .id = NULL,
            .description = NULL,
            .documentation = NULL,
            .xRef = NULL,
            .annotation = NULL,
            .type = NULL,
            .interRunCalibrator = NULL,
            .quantity = NULL,
            .calibratorSample = NULL,
            .cdnaSynthesisMethod = NULL,
            .templateQuantity = NULL
          ),
          active = list(
            id = function(id) {
              if (missing(id))
                return(private$.id)
              assert_that(is.type(id, idType))
              private$.id <- id
            },
            description = function(description) {
              if (missing(description))
                return(private$.description)
              assert_that(is.opt.string(description))
              private$.description <- description
            },
            documentation = function(documentation) {
              if (missing(documentation))
                return(private$.documentation)
              assert_that(is.opt.list.type(documentation,
                                           idReferencesType))
              private$.documentation <- documentation
            },
            xRef = function(xRef) {
              if (missing(xRef))
                return(private$.xRef)
              assert_that(is.opt.list.type(xRef,
                                           xRefType))
              private$.xRef <- xRef
            },
            annotation = function(annotation) {
              if (missing(annotation))
                return(private$.annotation)
              assert_that(is.opt.list.type(annotation,
                                           annotationType))
              private$.annotation <- with.names(annotation,
                                                quote(.$property))
            },
            type = function(type) {
              if (missing(type))
                return(private$.type)
              assert_that(is.type(type,
                                  sampleTypeType))
              private$.type <- type
            },
            interRunCalibrator = function(interRunCalibrator) {
              if (missing(interRunCalibrator))
                return(private$.interRunCalibrator)
              assert_that(is.opt.flag(interRunCalibrator))
              private$.interRunCalibrator <- interRunCalibrator
            },
            quantity = function(quantity) {
              if (missing(quantity))
                return(private$.quantity)
              assert_that(is.opt.type(quantity,
                                      quantityType))
              private$.quantity <- quantity
            },
            calibratorSample = function(calibratorSample) {
              if (missing(calibratorSample))
                return(private$.calibratorSample)
              assert_that(is.opt.flag(calibratorSample))
              private$.calibratorSample <- calibratorSample
            },
            cdnaSynthesisMethod = function(cdnaSynthesisMethod) {
              if (missing(cdnaSynthesisMethod))
                return(private$.cdnaSynthesisMethod)
              assert_that(is.opt.type(cdnaSynthesisMethod,
                                      cdnaSynthesisMethodType))
              private$.cdnaSynthesisMethod <- cdnaSynthesisMethod
            },
            templateQuantity = function(templateQuantity) {
              if (missing(templateQuantity))
                return(private$.templateQuantity)
              assert_that(is.opt.type(templateQuantity,
                                      templateQuantityType))
              private$.templateQuantity <- templateQuantity
            }
          ))

# oligoType ------------------------------------------------------------

#' oligoType R6 class.
#' 
#' Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{oligoType$new(threePrimeTag = NULL, 
#'   fivePrimeTag = NULL, sequence)}
#'   
#'   @section Fields: \describe{  
#'   \item{\code{threePrimeTag}}{\link[assertthat]{is.string}. Description of
#'   three prime modification (if present). }
#'   \item{\code{fivePrimeTag}}{\link[assertthat]{is.string}. Description of
#'   five prime modification (if present).}
#'   \item{\code{sequence}}{\link[assertthat]{is.string}.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
oligoType <- 
  R6Class("oligoType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(threePrimeTag = NULL,
                                  fivePrimeTag = NULL,
                                  sequence) {
              assert_that(is.opt.string(threePrimeTag))
              assert_that(is.opt.string(fivePrimeTag))
              assert_that(is.string(sequence))
              private$.threePrimeTag <- threePrimeTag
              private$.fivePrimeTag <- fivePrimeTag
              private$.sequence <- sequence
            }
          ),
          private = list(
            .threePrimeTag = NULL,
            .fivePrimeTag = NULL,
            .sequence = NULL
          ),
          active = list(
            threePrimeTag = function(threePrimeTag) {
              if (missing(threePrimeTag))
                return(private$.threePrimeTag)
              assert_that(is.opt.string(threePrimeTag))
              private$.threePrimeTag <- threePrimeTag
            },
            fivePrimeTag = function(fivePrimeTag) {
              if (missing(fivePrimeTag))
                return(private$.fivePrimeTag)
              assert_that(is.opt.string(fivePrimeTag))
              private$.fivePrimeTag <- fivePrimeTag
            },
            sequence = function(sequence) {
              if (missing(sequence))
                return(private$.sequence)
              assert_that(is.string(sequence))
              private$.sequence <- sequence
            }
          ))

# sequencesType ------------------------------------------------------------

#' sequencesType R6 class.
#' 
#' Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{sequencesType$new(forwardPrimer = NULL, reversePrimer = NULL, probe1 = NULL, probe2 = NULL, amplicon = NULL)}
#'   
#'   @section Fields: \describe{  
#' \item{\code{forwardPrimer}}{\link{oligoType}.}
#' \item{\code{reversePrimer}}{\link{oligoType}.}
#' \item{\code{probe1}}{\link{oligoType}.}
#' \item{\code{probe2}}{\link{oligoType}.}
#' \item{\code{amplicon}}{\link{oligoType}.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
sequencesType <- 
  R6Class("sequencesType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(forwardPrimer = NULL,
                                  reversePrimer = NULL,
                                  probe1 = NULL,
                                  probe2 = NULL,
                                  amplicon = NULL) {
              assert_that(is.opt.type(forwardPrimer,
                                      oligoType))
              assert_that(is.opt.type(reversePrimer,
                                      oligoType))
              assert_that(is.opt.type(probe1,
                                      oligoType))
              assert_that(is.opt.type(probe2,
                                      oligoType))
              assert_that(is.opt.type(amplicon,
                                      oligoType))
              private$.forwardPrimer <- forwardPrimer
              private$.reversePrimer <- reversePrimer
              private$.probe1 <- probe1
              private$.probe2 <- probe2
              private$.amplicon <- amplicon
            }
          ),
          private = list(
            .forwardPrimer = NULL,
            .reversePrimer = NULL,
            .probe1 = NULL,
            .probe2 = NULL,
            .amplicon = NULL
          ),
          active = list(
            forwardPrimer = function(forwardPrimer) {
              if (missing(forwardPrimer))
                return(private$.forwardPrimer)
              assert_that(is.opt.type(forwardPrimer,
                                      oligoType))
              private$.forwardPrimer <- forwardPrimer
            },
            reversePrimer = function(reversePrimer) {
              if (missing(reversePrimer))
                return(private$.reversePrimer)
              assert_that(is.opt.type(reversePrimer,
                                      oligoType))
              private$.reversePrimer <- reversePrimer
            },
            probe1 = function(probe1) {
              if (missing(probe1))
                return(private$.probe1)
              assert_that(is.opt.type(probe1,
                                      oligoType))
              private$.probe1 <- probe1
            },
            probe2 = function(probe2) {
              if (missing(probe2))
                return(private$.probe2)
              assert_that(is.opt.type(probe2,
                                      oligoType))
              private$.probe2 <- probe2
            },
            amplicon = function(amplicon) {
              if (missing(amplicon))
                return(private$.amplicon)
              assert_that(is.opt.type(amplicon,
                                      oligoType))
              private$.amplicon <- amplicon
            }
          ))

# commercialAssayType ------------------------------------------------------------

#' commercialAssayType R6 class.
#' 
#' For some commercial assays, the primer sequences may be unknown. This element 
#' allows to describe commercial assays.Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{commercialAssayType$new(company, orderNumber)}
#'   
#'   @section Fields: \describe{  
#' \item{\code{company}}{\link[assertthat]{is.string}.}
#' \item{\code{orderNumber}}{\link[assertthat]{is.string}.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
commercialAssayType <- 
  R6Class("commercialAssayType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(company,
                                  orderNumber) {
              assert_that(is.string(company))
              assert_that(is.string(orderNumber))
              private$.company <- company
              private$.orderNumber <- orderNumber
            }
          ),
          private = list(
            .company = NULL,
            .orderNumber = NULL
          ),
          active = list(
            company = function(company) {
              if (missing(company))
                return(private$.company)
              assert_that(is.string(company))
              private$.company <- company
            },
            orderNumber = function(orderNumber) {
              if (missing(orderNumber))
                return(private$.orderNumber)
              assert_that(is.string(orderNumber))
              private$.orderNumber <- orderNumber
            }
          ))

# targetTypeType ------------------------------------------------------------

#' targetTypeType R6 class.
#' 
#' Can take values:
#' \describe{
#' \item{ref}{reference target}
#' \item{toi}{target of interest}
#' } 
#' Inherits: \link{enumType}.
#' 
#' @section Initialization: \code{targetTypeType$new(value)}
#' 
#'   @section Fields: \describe{  
#' \item{\code{value}}{\link[assertthat]{is.string}.}
#' }
#'  
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
targetTypeType <- 
  R6Class("targetTypeType",
          # class = FALSE,
          inherit = enumType,
          private = list(
            .levels = c("ref",
                        "toi")
          )
  )

# targetType ------------------------------------------------------------

#' targetType R6 class.
#' 
#' A target is a defined PCR reaction. PCR reactions for the same gene which 
#' differ in primer sequences are considered different targets.Inherits: 
#' \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{targetType$new(id, description = NULL, 
#'   documentation = NULL, xRef = NULL, type, amplificationEfficiencyMethod = 
#'   NULL, amplificationEfficiency = NULL, amplificationEfficiencySE = NULL, 
#'   detectionLimit = NULL, dyeId, sequences = NULL, commercialAssay = NULL)}
#'   
#' @section Fields: \describe{ \item{\code{id}}{\link{idType}.} 
#'   \item{\code{description}}{\link[assertthat]{is.string}.} 
#'   \item{\code{documentation}}{\code{list} of \link{idReferencesType}.} 
#'   \item{\code{xRef}}{\code{list} of \link{xRefType}.} 
#'   \item{\code{type}}{\link{targetTypeType}.} 
#'   \item{\code{amplificationEfficiencyMethod}}{\link[assertthat]{is.string}.} 
#'   \item{\code{amplificationEfficiency}}{\link[base]{double}.} 
#'   \item{\code{amplificationEfficiencySE}}{\link[base]{double}.} 
#'   \item{\code{detectionLimit}}{\link[base]{double}.} 
#'   \item{\code{dyeId}}{\link{idReferencesType}.} 
#'   \item{\code{sequences}}{\link{sequencesType}.} 
#'   \item{\code{commercialAssay}}{\link{commercialAssayType}.} }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
targetType <- 
  R6Class("targetType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(id,
                                  description = NULL,
                                  documentation = NULL,
                                  xRef = NULL,
                                  type,
                                  amplificationEfficiencyMethod = NULL,
                                  amplificationEfficiency = NULL,
                                  amplificationEfficiencySE = NULL,
                                  detectionLimit = NULL,
                                  dyeId,
                                  sequences = NULL,
                                  commercialAssay = NULL) {
              assert_that(is.type(id, idType))
              assert_that(is.opt.string(description))
              assert_that(is.opt.list.type(documentation,
                                           idReferencesType))
              assert_that(is.opt.list.type(xRef,
                                           xRefType))
              assert_that(is.type(type,
                                  targetTypeType))
              assert_that(is.opt.string(amplificationEfficiencyMethod))
              assert_that(is.opt.double(amplificationEfficiency))
              assert_that(is.opt.double(amplificationEfficiencySE))
              assert_that(is.opt.double(detectionLimit))
              assert_that(is.type(dyeId,
                                  idReferencesType))
              assert_that(is.opt.type(sequences,
                                      sequencesType))
              assert_that(is.opt.type(commercialAssay,
                                      commercialAssayType))
              
              private$.id <- id
              private$.description <- description
              private$.documentation <- documentation
              private$.xRef <- xRef
              private$.type <- type
              private$.amplificationEfficiencyMethod <- amplificationEfficiencyMethod
              private$.amplificationEfficiency <- amplificationEfficiency
              private$.amplificationEfficiencySE <- amplificationEfficiencySE
              private$.detectionLimit <- detectionLimit
              private$.dyeId <- dyeId
              private$.sequences <- sequences
              private$.commercialAssay <- commercialAssay
            }
          ),
          private = list(
            .id = NULL,
            .description = NULL,
            .documentation = NULL,
            .xRef = NULL,
            .type = NULL,
            .amplificationEfficiencyMethod = NULL,
            .amplificationEfficiency = NULL,
            .amplificationEfficiencySE = NULL,
            .detectionLimit = NULL,
            .dyeId = NULL,
            .sequences = NULL,
            .commercialAssay = NULL
          ),
          active = list(
            id = function(id) {
              if (missing(id))
                return(private$.id)
              assert_that(is.type(id, idType))
              private$.id <- id
            },
            description = function(description) {
              if (missing(description))
                return(private$.description)
              assert_that(is.opt.string(description))
              private$.description <- description
            },
            documentation = function(documentation) {
              if (missing(documentation))
                return(private$.documentation)
              assert_that(is.opt.list.type(documentation,
                                           idReferencesType))
              private$.documentation <- documentation
            },
            xRef = function(xRef) {
              if (missing(xRef))
                return(private$.xRef)
              assert_that(is.opt.list.type(xRef,
                                           xRefType))
              private$.xRef <- xRef
            },
            type = function(type) {
              if (missing(type))
                return(private$.type)
              assert_that(is.type(type,
                                  targetTypeType))
              private$.type <- type
            },
            amplificationEfficiencyMethod = 
              function(amplificationEfficiencyMethod) {
                if (missing(amplificationEfficiencyMethod))
                  return(private$.amplificationEfficiencyMethod)
                assert_that(is.opt.string(amplificationEfficiencyMethod))
                private$.amplificationEfficiencyMethod <- amplificationEfficiencyMethod
              },
            amplificationEfficiency = function(amplificationEfficiency) {
              if (missing(amplificationEfficiency))
                return(private$.amplificationEfficiency)
              assert_that(is.opt.double(amplificationEfficiency))
              private$.amplificationEfficiency <- amplificationEfficiency
            },
            amplificationEfficiencySE = function(amplificationEfficiencySE) {
              if (missing(amplificationEfficiencySE))
                return(private$.amplificationEfficiencySE)
              assert_that(is.opt.double(amplificationEfficiencySE))
              private$.amplificationEfficiencySE <- amplificationEfficiencySE
            },
            detectionLimit = function(detectionLimit) {
              if (missing(detectionLimit))
                return(private$.detectionLimit)
              assert_that(is.opt.double(detectionLimit))
              private$.detectionLimit <- detectionLimit
            },
            dyeId = function(dyeId) {
              if (missing(dyeId))
                return(private$.dyeId)
              assert_that(is.type(dyeId,
                                  idReferencesType))
              private$.dyeId <- dyeId
            },
            sequences = function(sequences) {
              if (missing(sequences))
                return(private$.sequences)
              assert_that(is.opt.type(sequences,
                                      sequencesType))
              private$.sequences <- sequences
            },
            commercialAssay = function(commercialAssay) {
              if (missing(commercialAssay))
                return(private$.commercialAssay)
              assert_that(is.opt.type(commercialAssay,
                                      commercialAssayType))
              private$.commercialAssay <- commercialAssay
            }
          ))

# # dpAmpCurveType ------------------------------------------------------------
# dpAmpCurveType <- 
#   R6Class("dpAmpCurveType",
#           # class = FALSE,
#           inherit = rdmlBaseType,
#           public = list(
#             initialize = function(cyc,
#                                   tmp = NULL,
#                                   fluor) {
#               assert_that(is.float(cyc))
#               assert_that(is.opt.double(tmp))
#               assert_that(is.float(fluor))
#               private$.cyc <- cyc
#               private$.tmp <- tmp
#               private$.fluor <- fluor
#             },
#             AsVector = function() {
#               c(cyc = private$.cyc,
#                 {
#                   if (!is.null(private$.tmp))
#                     c(tmp = private$.tmp)
#                   },
#                 fluor = private$.fluor)
#             }
#           ),
#           private = list(
#             .cyc = NULL,
#             .tmp = NULL,
#             .fluor = NULL
#           ),
#           active = list(
#             cyc = function(cyc) {
#               if (missing(cyc))
#                 return(private$.cyc)
#               assert_that(is.float(cyc))
#               private$.cyc <- cyc
#             },
#             tmp = function(tmp) {
#               if (missing(tmp))
#                 return(private$.tmp)
#               assert_that(is.opt.double(tmp))
#               private$.tmp <- tmp
#             },
#             fluor = function(fluor) {
#               if (missing(fluor))
#                 return(private$.fluor)
#               assert_that(is.float(fluor))
#               private$.fluor <- fluor
#             }
#           ))
# 

# # dpMeltingCurveType ------------------------------------------------------------
# dpMeltingCurveType <- 
#   R6Class("dpMeltingCurveType",
#           # class = FALSE,
#           inherit = rdmlBaseType,
#           public = list(
#             initialize = function(tmp,
#                                   fluor) {
#               assert_that(is.float(tmp))
#               assert_that(is.float(fluor))
#               private$.tmp <- tmp
#               private$.fluor <- fluor
#             },
#             AsVector = function() {
#               c(cyc = private$.tmp,
#                 fluor = private$.fluor)
#             }
#           ),
#           private = list(
#             .cyc = NULL,
#             .tmp = NULL,
#             .fluor = NULL
#           ),
#           active = list(
#             tmp = function(tmp) {
#               if (missing(tmp))
#                 return(private$.tmp)
#               assert_that(is.float(tmp))
#               private$.tmp <- tmp
#             },
#             fluor = function(fluor) {
#               if (missing(fluor))
#                 return(private$.fluor)
#               assert_that(is.float(fluor))
#               private$.fluor <- fluor
#             }
#           ))

# adpsType ------------------------------------------------------------

#' adpsType R6 class.
#' 
#' Contains of amplification data points \code{matrix} -- single data points 
#' measured during amplification. \code{Matrix} columns: \describe{ 
#' \item{cyc}{(every point must be unique) Cycle - The PCR cycle at which data 
#' point was collected.} \item{tmp}{(optional) Temperature - The temperature in 
#' degrees Celsius at the time of measurement.} \item{fluor}{Fluorescence - The 
#' fluorescence intensity measured without any correction. The fluorescence 
#' intensity must not be baseline corrected.}} Inherits: 
#' \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{adpsType$new(fpoints)}
#' 
#' @section Fields: \describe{    
#'   \item{\code{fpoints}}{\link[base]{matrix}. Matrix with amplification data
#'   points.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
adpsType <- 
  R6Class("adpsType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(fpoints) {
              assert_that(is.double.matrix(fpoints))
              assert_that(has.only.names(fpoints,
                                         c("cyc", "tmp", "fluor")))
              private$.fpoints <- fpoints
            },
            .asXMLnodes = function(node.name) {
              #               newXMLNode(
              #                 name = node.name,
              #                 .children = {
              alply(private$.fpoints,
                    1,
                    function(fpoints.row) {
                      newXMLNode(name = "adp",
                                 .children = 
                                   list(
                                     newXMLNode(
                                       name = "cyc",
                                       text = fpoints.row["cyc"]),
                                     {
                                       if(!is.na(fpoints.row["tmp"]))
                                         newXMLNode(
                                           name = "tmp",
                                           text = fpoints.row["tmp"])
                                     },
                                     newXMLNode(
                                       name = "fluor",
                                       text = fpoints.row["fluor"])
                                   ))
                    })
              # })
            }
          ),
          private = list(
            .fpoints = NULL
          ),
          active = list(
            fpoints = function(fpoints) {
              if (missing(fpoints))
                return(private$.fpoints)
              assert_that(is.double.matrix(fpoints))
              assert_that(has.only.names(fpoints,
                                         c("cyc", "tmp", "fluor")))
              private$.fpoints <- fpoints
            }
          ))

# mdpsType ------------------------------------------------------------

#' mdpsType R6 class.
#' 
#' Contains of melting data points \code{matrix} -- single data points measured
#' during amplification. \code{Matrix} columns: \describe{ \item{tmp}{(every
#' point must be unique) Temperature - The temperature in degrees Celsius at the
#' time of measurement.} \item{fluor}{Fluorescence - The fluorescence intensity 
#' measured without any correction. The fluorescence intensity must not be 
#' baseline corrected.}} Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{mdpsType$new(fpoints)}
#'   
#'   @section Fields: \describe{  
#' \item{\code{fpoints}}{\link[base]{matrix}. Matrix with amplification data points.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
mdpsType <- 
  R6Class("mdpsType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(fpoints) {
              assert_that(is.double.matrix(fpoints))
              assert_that(has.only.names(fpoints,
                                         c("tmp", "fluor")))
              private$.fpoints <- fpoints
            },
            .asXMLnodes = function(node.name) {
              alply(private$.fpoints,
                    1,
                    function(fpoints.row) {
                      newXMLNode(name = "mdp",
                                 .children = 
                                   list(
                                     newXMLNode(
                                       name = "tmp",
                                       text = fpoints.row["tmp"]),
                                     newXMLNode(
                                       name = "fluor",
                                       text = fpoints.row["fluor"])
                                   ))
                    })
            }
          ),
          private = list(
            .fpoints = NULL
          ),
          active = list(
            fpoints = function(fpoints) {
              if (missing(fpoints))
                return(private$.fpoints)
              assert_that(is.double.matrix(fpoints))
              assert_that(has.only.names(fpoints,
                                         c("tmp", "fluor")))
              private$.fpoints <- fpoints
            }
          ))

# dataType ------------------------------------------------------------

#' dataType R6 class.
#' 
#' Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{dataType$new(tar, cq = NULL, excl = NULL, adp 
#'   = NULL, mdp = NULL, endPt = NULL, bgFluor = NULL, bgFluorSlp = NULL, 
#'   quantFluor = NULL)}
#'   
#' @section Fields: \describe{ \item{\code{tar}}{\link{idReferencesType}. 
#'   TargetID - A reference to a target.} \item{\code{cq}}{\link[base]{double}. 
#'   Quantification cycle - The calculated fractional PCR cycle used for 
#'   downstream quantification. Negative values are used to express following 
#'   conditions: Not Available: -1.0 } 
#'   \item{\code{excl}}{\link[assertthat]{is.string}. Excluded - If present, 
#'   this entry should not be evaluated. Do not set this element to false if 
#'   this entry is valid, leave the entire element out instead. It may contain a
#'   string with reason for exclusion. Several reasons for exclusion should be 
#'   seperated by semicolons ";".} \item{\code{adp}}{\link{adpsType}.} 
#'   \item{\code{mdp}}{\link{mdpsType}.} 
#'   \item{\code{endPt}}{\link[base]{double}. End point - Result of an endpoint 
#'   } measurement. \item{\code{bgFluor}}{\link[base]{double}. Background 
#'   fluorescence - The y-intercept of the baseline trend based on the estimated
#'   background fluorescence. } \item{\code{bgFluorSlp}}{\link[base]{double}. 
#'   Background fluorescence slope - The slope of the baseline trend based on 
#'   the estimated background fluorescence. The element should be absent to 
#'   indicate a slope of 0.0; If this element is present without the bgFluor 
#'   element it should be ignored. } 
#'   \item{\code{quantFluor}}{\link[base]{double}. Quantification flourescence -
#'   The fluorescence value corresponding to the treshold line.} }
#'   
#' @section Methods: \describe{\item{\code{AsDataFrame(dp.type = 
#'   "adp")}}{Represents amplification (\code{dp.type = "adp"}) or melting 
#'   (\code{dp.type = "mdp"}) data points as \code{data.frame}}}
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
dataType <- 
  R6Class("dataType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(tar,
                                  cq = NULL,
                                  excl = NULL,
                                  adp = NULL,
                                  mdp = NULL,
                                  endPt = NULL,
                                  bgFluor = NULL,
                                  bgFluorSlp = NULL,
                                  quantFluor = NULL) {
              assert_that(is.type(tar,
                                  idReferencesType))
              assert_that(is.opt.double(cq))
              assert_that(is.opt.string(excl))
              #               assert_that(is.opt.list.type(adp,
              #                                            dpAmpCurveType))
              #               assert_that(is.opt.list.type(mdp,
              #                                            dpMeltingCurveType))
              assert_that(is.opt.type(adp,
                                      adpsType))
              assert_that(is.opt.type(mdp,
                                      mdpsType))
              assert_that(is.opt.double(endPt))
              assert_that(is.opt.double(bgFluor))
              assert_that(is.opt.double(bgFluorSlp))
              assert_that(is.opt.double(quantFluor))
              
              private$.tar <- tar
              private$.cq <- cq
              private$.excl <- excl
              private$.adp <- adp
              private$.mdp <- mdp
              private$.endPt <- endPt 
              private$.bgFluor <- bgFluor
              private$.bgFluorSlp <- bgFluorSlp
              private$.quantFluor <- quantFluor
            },
            AsDataFrame = function(dp.type = "adp") {
              self[[dp.type]]$fpoints %>% 
                as.data.frame %>% 
                select(ifelse(dp.type == "adp",
                              cyc,
                              tmp),
                       fluor)
            }
          ),
          private = list(
            .tar = NULL,
            .cq = NULL,
            .excl = NULL,
            .adp = NULL,
            .mdp = NULL,
            .endPt = NULL,
            .bgFluor = NULL,
            .bgFluorSlp = NULL,
            .quantFluor = NULL
          ),
          active = list(
            tar = function(tar) {
              if (missing(tar))
                return(private$.tar)
              assert_that(is.type(tar,
                                  idReferencesType))
              private$.tar <- tar
            },
            cq = function(cq) {
              if (missing(cq))
                return(private$.cq)
              assert_that(is.opt.double(cq))
              private$.cq <- cq
            },
            excl = function(excl) {
              if (missing(excl))
                return(private$.excl)
              assert_that(is.opt.string(excl))
              private$.excl <- excl
            },
            adp = function(adp) {
              if (missing(adp))
                return(private$.adp)
              assert_that(is.opt.type(adp,
                                      adpsType))
              private$.adp <- adp
            },
            mdp = function(mdp) {
              if (missing(mdp))
                return(private$.mdp)
              assert_that(is.opt.type(mdp,
                                      mdpsType))
              private$.mdp <- mdp
            },
            endPt = function(endPt) {
              if (missing(endPt))
                return(private$.endPt)
              assert_that(is.opt.double(endPt))
              private$.endPt <- endPt
            },
            bgFluor = function(bgFluor) {
              if (missing(bgFluor))
                return(private$.bgFluor)
              assert_that(is.opt.double(bgFluor))
              private$.bgFluor <- bgFluor
            },
            bgFluorSlp = function(bgFluorSlp) {
              if (missing(bgFluorSlp))
                return(private$.bgFluorSlp)
              assert_that(is.opt.double(bgFluorSlp))
              private$.bgFluorSlp <- bgFluorSlp
            },
            quantFluor = function(quantFluor) {
              if (missing(quantFluor))
                return(private$.quantFluor)
              assert_that(is.opt.double(quantFluor))
              private$.quantFluor <- quantFluor
            }
          ))

# reactType ------------------------------------------------------------

#' reactType R6 class.
#' 
#' A reaction is an independent chemical reaction corresponding for example to a
#' well in a 96 well plate, a capillary in a rotor, a through-hole on an array, 
#' etc. Inherits: \link{rdmlBaseType}.
#' 
#' The ID of this reaction
#' 
#' Schemas : \itemize{ \item rotor : assign IDs according to the position of the
#' sample on the rotor (1 for the 1st sample, 2 for the 2nd, ...) \item plate 
#' (96/384/1536 well) : the IDs are assigned in a row-first/column-second 
#' manner. For each row, the samples are numbered according to the increasing 
#' column number. At the end of a row, the numbering starts at the first column 
#' of the next row. An example for this type of plate can be found below : 
#' \tabular{lllll}{ \tab 1  \tab 2  \tab 3 \tab ... \cr A   \tab 1  \tab 2  \tab
#' 3 \tab     \cr B   \tab 13 \tab 14 \tab   \tab     \cr ... \tab    \tab \tab 
#' \tab    } or \tabular{lllll}{ \tab 1  \tab 2  \tab 3 \tab ... \cr 1 \tab 1 
#' \tab 2  \tab 3 \tab     \cr 2   \tab 13 \tab 14 \tab   \tab     \cr ... \tab 
#' \tab    \tab   \tab    }
#' 
#' \item multi-array plate (BioTrove) : the IDs are assigned in a 
#' row-first/column-second manner, ignoring the organisation of sub-arrays. For 
#' each row, the samples are numbered according to the increasing column number.
#' At the end of a row, the the next row. An example for this type of plate can 
#' be found below : todo... }
#' 
#' @section Initialization: \code{reactType$new(id, sample, data = NULL)}
#'   
#'   @section Fields: \describe{
#'   \item{\code{id}}{\link{reactIdType}. See 'Details'.} 
#'   \item{\code{sample}}{\link{idReferencesType}. SampleID - A reference to a
#'   sample.} 
#'   \item{\code{data}}{\code{list} of \link{dataType}.}
#'   }
#'   
#' @section Methods: \describe{\item{\code{AsDataFrame(dp.type = 
#'   "adp")}}{Represents amplification (\code{dp.type = "adp"}) or melting 
#'   (\code{dp.type = "mdp"}) data points of all targets as one 
#'   \code{data.frame}} \item{\code{Position(pcrformat)}}{Converts \code{react 
#'   id} to thew human readable form (i.e. '13' -> 'B1'). \code{pcrFormat} is 
#'   \code{pcrFormatType}. Only for 'ABC' '123' format!}}
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
reactType <- 
  R6Class("reactType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(id,
                                  sample,
                                  data = NULL) {
              assert_that(is.type(id,
                                  reactIdType))
              assert_that(is.type(sample,
                                  idReferencesType))
              assert_that(is.opt.list.type(data,
                                           dataType))
              private$.id <- id
              private$.sample <- sample
              private$.data <- with.names(data,
                                          quote(.$tar$id))
            },
            AsDataFrame = function(dp.type = "adp",
                                   long.table = FALSE) {
              assert_that(is.string(dp.type))
              out <-
                private$.data %>% 
                ldply(function(data)
                  data$AsDataFrame(dp.type) %>% 
                    cbind(tar = data$tar$id),
                  .id = "tar")
              
              if (long.table == FALSE) out <- out %>% spread(tar,
                                                             fluor) 
              out
            },
            Position = function(pcrFormat) {
              assert_that(is.type(pcrFormat,
                                  pcrFormatType))
              stopifnot(pcrFormat$rowLabel$value == "ABC",
                        pcrFormat$columnLabel$value == "123")
              sprintf("%s%02i",
                      LETTERS[(private$.id$id - 1) %/% pcrFormat$columns + 1],
                      as.integer((private$.id$id - 1) %% pcrFormat$columns + 1))
            }
          ),
          private = list(
            .id = NULL,
            .sample = NULL,
            .data = NULL
          ),
          active = list(
            id = function(id) {
              if (missing(id))
                return(private$.id)
              assert_that(is.type(id,
                                  reactIdType))
              private$.id <- id
            },
            sample = function(sample) {
              if (missing(sample))
                return(private$.sample)
              assert_that(is.type(sample,
                                  idReferencesType))
              private$.sample <- sample
            },
            data = function(data) {
              if (missing(data))
                return(private$.data)
              assert_that(is.opt.list.type(data,
                                           dataType))
              private$.data <- with.names(data,
                                          quote(.$tar$id))
            }
          ))

# dataCollectionSoftwareType ------------------------------------------------------------

#' dataCollectionSoftwareType R6 class.
#' 
#' Software name and version used to collect and analyze the data.Inherits: 
#' \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{dataCollectionSoftwareType$new(name, version)}
#'
#'   @section Fields: \describe{   
#' \item{\code{name}}{\link[assertthat]{is.string}.}
#' \item{\code{version}}{\link[assertthat]{is.string}.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
dataCollectionSoftwareType <- 
  R6Class("dataCollectionSoftwareType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(name,
                                  version) {
              assert_that(is.string(name))
              assert_that(is.string(version))
              private$.name <- name
              private$.version <- version
            }
          ),
          private = list(
            .name = NULL,
            .version = NULL
          ),
          active = list(
            name = function(name) {
              if (missing(name))
                return(private$.name)
              assert_that(is.string(name))
              private$.name <- name
            },
            version = function(version) {
              if (missing(version))
                return(private$.version)
              assert_that(is.string(version))
              private$.version <- version
            }
          ))

# cqDetectionMethodType ------------------------------------------------------------

#' cqDetectionMethodType R6 class.
#' 
#' The method used to determine the Cq value.
#' Can take values:
#' \describe{
#' \item{automated threshold and baseline settings}{}
#' \item{manual threshold and baseline settings}{}
#' \item{second derivative maximum}{}
#' \item{other}{}
#' }  
#' Inherits: \link{enumType}.
#' 
#' @section Initialization: \code{cqDetectionMethodType$new(value)}
#'   
#'   @section Fields: \describe{
#' \item{\code{value}}{\link[assertthat]{is.string}.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
cqDetectionMethodType <- 
  R6Class("cqDetectionMethodType",
          # class = FALSE,
          inherit = enumType,
          private = list(
            .levels = c("automated threshold and baseline settings",
                        "manual threshold and baseline settings",
                        "second derivative maximum",
                        "other")
          )
  )

# labelFormatType ------------------------------------------------------------

#' labelFormatType R6 class.
#' 
#' Label used for \link{pcrFormatType}.
#' Can take values:
#' \describe{
#' \item{ABC}{}
#' \item{123}{}
#' \item{A1a1}{}
#' }  
#' Inherits: \link{enumType}.
#' 
#' @section Initialization: \code{labelFormatType$new(value)}
#'   
#'   @section Fields: \describe{
#' \item{\code{value}}{\link[assertthat]{is.string}.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
labelFormatType <- 
  R6Class("labelFormatType",
          # class = FALSE,
          inherit = enumType,
          private = list(
            .levels = c("ABC",
                        "123",
                        "A1a1")
          )
  )

# pcrFormatType ------------------------------------------------------------

#' pcrFormatType R6 class.
#' 
#' The format of the run - This allows the software to display the data 
#' according to the qPCR instrument run format.Inherits: 
#' \link{rdmlBaseType}.
#' 
#' Rotor formats always have 1 column; rows correspond to the number of places 
#' in the rotor. Values for common formats are: \tabular{lllll}{ Format \tab
#' rows \tab columns \tab rowLabel \tab columnLabel \cr single-well \tab 1   
#' \tab 1       \tab 123      \tab 123         \cr 48-well plate \tab 6    \tab
#' 8       \tab ABC      \tab 123         \cr 96-well plate \tab 8    \tab 12   
#' \tab ABC      \tab 123         \cr 384-well plate \tab 16   \tab 24      \tab
#' ABC      \tab 123         \cr 1536-well plate \tab 32   \tab 48      \tab ABC
#' \tab 123         \cr 3072-well array \tab 32   \tab 96      \tab A1a1    
#' \tab A1a1        \cr 5184-well chip \tab 72   \tab 72      \tab ABC      \tab
#' 123         \cr 32-well rotor \tab 32   \tab 1       \tab 123      \tab 123  
#' \cr 72-well rotor \tab 72   \tab 1       \tab 123      \tab 123         \cr
#' 100-well rotor \tab 100  \tab 1       \tab 123      \tab 123         \cr free
#' format \tab -1   \tab 1       \tab 123      \tab 123 } If rows are -1 then
#' the software should not try to reconstruct a plate and just display all react
#' data in list (1 column) form. If columns is 1 then the software should not
#' display a column label.
#' 
#' @section Initialization: \code{pcrFormatType$new(rows, columns, rowLabel, columnLabel)}
#'   
#'   @section Fields: \describe{
#' \item{\code{rows}}{\link[assertthat]{is.count}.}
#' \item{\code{columns}}{\link[assertthat]{is.count}.}
#' \item{\code{rowLabel}}{\link{labelFormatType}.}
#' \item{\code{columnLabel}}{\link{labelFormatType}.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
pcrFormatType <- 
  R6Class("pcrFormatType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(rows,
                                  columns,
                                  rowLabel,
                                  columnLabel) {
              assert_that(is.count(rows))
              assert_that(is.count(columns))
              assert_that(is.type(rowLabel,
                                  labelFormatType))
              assert_that(is.type(columnLabel,
                                  labelFormatType))
              private$.rows <- rows
              private$.columns <- columns
              private$.rowLabel <- rowLabel
              private$.columnLabel <- columnLabel
            }
          ),
          private = list(
            .rows = NULL,
            .columns = NULL,
            .rowLabel = NULL,
            .columnLabel = NULL
          ),
          active = list(
            rows = function(rows) {
              if (missing(rows))
                return(private$.rows)
              assert_that(is.count(rows))
              private$.rows <- rows
            },
            columns = function(columns) {
              if (missing(columns))
                return(private$.columns)
              assert_that(is.count(columns))
              private$.columns <- columns
            },
            rowLabel = function(rowLabel) {
              if (missing(rowLabel))
                return(private$.rowLabel)
              assert_that(is.type(rowLabel,
                                  labelFormatType))
              private$.rowLabel <- rowLabel
            },
            columnLabel = function(columnLabel) {
              if (missing(columnLabel))
                return(private$.columnLabel)
              assert_that(is.type(columnLabel,
                                  labelFormatType))
              private$.columnLabel <- columnLabel
            }
          ))

# runType ------------------------------------------------------------

#' runType R6 class.
#' 
#' A run is a set of reactions performed in one "run", for example one plate, 
#' one rotor, one array, one chip.Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{runType$new(id, description = NULL, 
#'   documentation = NULL, experimenter = NULL, instrument = NULL, 
#'   dataCollectionSoftware = NULL, backgroundDeterminationMethod = NULL, 
#'   cqDetectionMethod = NULL, thermalCyclingConditions = NULL, pcrFormat, 
#'   runDate = NULL, react = NULL)}
#'   
#' @section Fields: \describe{ \item{\code{id}}{\link{idType}.} 
#'   \item{\code{description}}{\link[assertthat]{is.string}.} 
#'   \item{\code{documentation}}{\code{list} of \link{idReferencesType}.} 
#'   \item{\code{experimenter}}{\code{list} of \link{idReferencesType}.} 
#'   \item{\code{instrument}}{\link[assertthat]{is.string}. Description of the 
#'   instrument used to aquire the data.} 
#'   \item{\code{dataCollectionSoftware}}{\link{dataCollectionSoftwareType}. 
#'   Description of the software used to analyze/collect the data.} 
#'   \item{\code{backgroundDeterminationMethod}}{\link[base]{double}. 
#'   Description of method used to determine the background. } 
#'   \item{\code{cqDetectionMethod}}{\link[base]{double}. Description of method 
#'   used to calculate the quantification cycle. } 
#'   \item{\code{thermalCyclingConditions}}{\link[base]{double}. The program 
#'   used to aquire the data.} \item{\code{pcrFormat}}{\link{adpsType}.} 
#'   \item{\code{runDate}}{\link{adpsType}. Date and time stamp when the data 
#'   was aquired.} \item{\code{react}}{\code{list} of \link{adpsType}.} }
#'   
#' @section Methods: \describe{\item{\code{AsDataFrame(dp.type = 
#'   "adp")}}{Represents amplification (\code{dp.type = "adp"}) or melting 
#'   (\code{dp.type = "mdp"}) data points as \code{data.frame}}}
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
runType <- 
  R6Class("runType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(id,
                                  description = NULL,
                                  documentation = NULL,
                                  experimenter = NULL,
                                  instrument = NULL,
                                  dataCollectionSoftware = NULL,
                                  backgroundDeterminationMethod = NULL,
                                  cqDetectionMethod = NULL,
                                  thermalCyclingConditions = NULL,
                                  pcrFormat,
                                  runDate = NULL,
                                  react = NULL) {
              assert_that(is.type(id,
                                  idType))
              assert_that(is.opt.string(description))
              assert_that(is.opt.list.type(documentation,
                                       idReferencesType))
              assert_that(is.opt.list.type(experimenter,
                                       idReferencesType))
              assert_that(is.opt.string(instrument))
              assert_that(is.opt.type(dataCollectionSoftware,
                                      dataCollectionSoftwareType))
              assert_that(is.opt.string(backgroundDeterminationMethod))
              assert_that(is.opt.type(cqDetectionMethod,
                                      cqDetectionMethodType))
              assert_that(is.opt.type(thermalCyclingConditions,
                                      idReferencesType))
              assert_that(is.type(pcrFormat,
                                  pcrFormatType))
              assert_that(is.opt.string(runDate)) # date time
              assert_that(is.opt.list.type(react,
                                       reactType))
              
              private$.id <- id
              private$.description <- description
              private$.documentation <- documentation
              private$.experimenter <- experimenter
              private$.instrument <- instrument
              private$.dataCollectionSoftware <- dataCollectionSoftware
              private$.backgroundDeterminationMethod <- backgroundDeterminationMethod
              private$.cqDetectionMethod <- cqDetectionMethod
              private$.thermalCyclingConditions <- thermalCyclingConditions
              private$.pcrFormat <- pcrFormat
              private$.runDate <- runDate
              private$.react <- with.names(react,
                                           quote(.$id$id))
            },
            AsDataFrame = function(dp.type = "adp",
                                   long.table = FALSE) {
              assert_that(is.string(dp.type))
              out <- 
                private$.react %>% 
                ldply(function(react)
                  
                  react$AsDataFrame(
                    dp.type = dp.type,
                    long.table = TRUE) %>% 
                    cbind(.,
                          sname = 
                            sprintf("%s_%s",
                                    react$id$id,
                                    react$sample$id)
                    ),
                  .id = ifelse(dp.type == "adp",
                               "cyc",
                               "tmp"))
              
              if (long.table == FALSE) out <- out %>% 
                tidyr::unite(sname_tar, sname, tar, sep = "_") %>% 
                tidyr::spread(sname_tar, fluor) #%>% 
              #                   arrange(ifelse(dp.type == "adp",
              #                                  cyc,
              #                                  tmp))
              out
            }
          ),
          private = list(
            .id = NULL,
            .description = NULL,
            .documentation = NULL,
            .experimenter = NULL,
            .instrument = NULL,
            .dataCollectionSoftware = NULL,
            .backgroundDeterminationMethod = NULL,
            .cqDetectionMethod = NULL,
            .thermalCyclingConditions = NULL,
            .pcrFormat = NULL,
            .runDate = NULL,
            .react = NULL
          ),
          active = list(
            id = function(id) {
              if (missing(id))
                return(private$.id)
              assert_that(is.type(id,
                                  idType))
              private$.id <- id
            },
            description = function(description) {
              if (missing(description))
                return(private$.description)
              assert_that(is.opt.string(description))
              private$.description <- description
            },
            documentation = function(documentation) {
              if (missing(documentation))
                return(private$.documentation)
              assert_that(is.opt.list.type(documentation,
                                       idReferencesType))
              private$.documentation <- documentation
            },
            experimenter = function(experimenter) {
              if (missing(experimenter))
                return(private$.experimenter)
              assert_that(is.opt.list.type(experimenter,
                                       idReferencesType))
              private$.experimenter <- experimenter
            },
            instrument = function(instrument) {
              if (missing(instrument))
                return(private$.instrument)
              assert_that(is.opt.string(instrument))
              private$.instrument <- instrument
            },
            dataCollectionSoftware = function(dataCollectionSoftware) {
              if (missing(dataCollectionSoftware))
                return(private$.dataCollectionSoftware)
              assert_that(is.opt.type(dataCollectionSoftware,
                                      dataCollectionSoftwareType))
              private$.dataCollectionSoftware <- dataCollectionSoftware
            },
            backgroundDeterminationMethod = function(backgroundDeterminationMethod) {
              if (missing(backgroundDeterminationMethod))
                return(private$.backgroundDeterminationMethod)
              assert_that(is.opt.string(backgroundDeterminationMethod))
              private$.backgroundDeterminationMethod <- backgroundDeterminationMethod
            },
            cqDetectionMethod = function(cqDetectionMethod) {
              if (missing(cqDetectionMethod))
                return(private$.cqDetectionMethod)
              assert_that(is.opt.type(cqDetectionMethod,
                                      cqDetectionMethodType))
              private$.cqDetectionMethod <- cqDetectionMethod
            },
            thermalCyclingConditions = function(thermalCyclingConditions) {
              if (missing(thermalCyclingConditions))
                return(private$.thermalCyclingConditions)
              assert_that(is.opt.type(thermalCyclingConditions,
                                      idReferencesType))
              private$.thermalCyclingConditions <- thermalCyclingConditions
            },
            pcrFormat = function(pcrFormat) {
              if (missing(pcrFormat))
                return(private$.pcrFormat)
              assert_that(is.type(pcrFormat,
                                  pcrFormatType))
              private$.pcrFormat <- pcrFormat
            },
            runDate = function(runDate) {
              if (missing(runDate))
                return(private$.runDate)
              assert_that(is.opt.string(runDate)) # date time
              private$.runDate <- runDate
            },
            react = function(react) {
              if (missing(react))
                return(private$.react)
              assert_that(is.opt.list.type(react,
                                       reactType))
              private$.react <- with.names(react,
                                           quote(.$id$id))
            }
          ))

# experimentType ------------------------------------------------------------

#' experimentType R6 class.
#' 
#' An experiment can contain several runs (\link{runType}).Inherits: 
#' \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{experimentType$new(id, description = NULL,
#'   documentation = NULL, run = NULL)}
#'
#'   @section Fields: \describe{   
#' \item{\code{id}}{\link{idType}.}
#' \item{\code{description}}{\link[assertthat]{is.string}.}
#' \item{\code{documentation}}{\code{list} of \link{idReferencesType}.}
#' \item{\code{run}}{\code{list} of \link{runType}.}
#' }
#'   
#' @section Methods: \describe{\item{\code{AsDataFrame(dp.type = "adp", 
#'   long.table = FALSE)}}{Represents amplification (\code{dp.type = "adp"}) or 
#'   melting (\code{dp.type = "mdp"}) data points as \code{data.frame}. 
#'   \code{long.table = TRUE} means that fluorescence data for all runs and 
#'   reacts will be at one collumn.}}
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
experimentType <- 
  R6Class("experimentType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(id,
                                  description = NULL,
                                  documentation = NULL,
                                  run = NULL) {
              assert_that(is.type(id,
                                  idType))
              assert_that(is.opt.string(description))
              assert_that(is.opt.list.type(documentation,
                                           idReferencesType))
              assert_that(is.opt.list.type(run,
                                           runType))
              
              private$.id <- id
              private$.description <- description
              private$.documentation <- documentation
              private$.run <- with.names(run,
                                         quote(.$id$id))
            },
            AsDataFrame = function(dp.type = "adp",
                                   long.table = FALSE) {
              assert_that(is.string(dp.type))
              assert_that(is.flag(long.table))
              out <-
                private$.run %>% 
                ldply(function(run)
                  run$AsDataFrame(
                    dp.type = dp.type,
                    long.table = TRUE) %>% 
                    cbind(.,  run = run$id$id,
                          row.names = NULL),
                  .id = ifelse(dp.type == "adp",
                               "cyc",
                               "tmp"))
              if (long.table == FALSE) out <- out %>%
                tidyr::unite(run_sname_tar, run, sname, tar, sep = "_") %>% 
                tidyr::spread(run_sname_tar, fluor)
              out
            }
          ),
          private = list(
            .id = NULL,
            .description = NULL,
            .documentation = NULL,
            .run = NULL
          ),
          active = list(
            id = function(id) {
              if (missing(id))
                return(private$.id)
              assert_that(is.type(id, idType))
              private$.id <- id
            },
            description = function(description) {
              if (missing(description))
                return(private$.description)
              assert_that(is.opt.string(description))
              private$.description <- description
            },
            documentation = function(documentation) {
              if (missing(documentation))
                return(private$.documentation)
              assert_that(is.opt.list.type(documentation,
                                           idReferencesType))
              private$.documentation <- documentation
            },
            run = function(run) {
              if (missing(run))
                return(private$.run)
              assert_that(is.opt.list.type(run,
                                           runType))
              private$.run <- with.names(run,
                                         quote(.$id$id))
            }
          ))

# lidOpenType ------------------------------------------------------------

#' lidOpenType R6 class.
#' 
#' This step waits for the user to open the lid and continues afterwards. It 
#' allows to stop the program and to wait for the user to add for example 
#' enzymes and continue the program afterwards. The temperature of the previous 
#' step is maintained.Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{lidOpenType$new()}
#'   
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
lidOpenType <- 
  R6Class("lidOpenType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(lidOpen) {
            }
          ),
          private = list(
            .lidOpen = NA
          ))

# pauseType ------------------------------------------------------------

#' pauseType R6 class.
#' 
#' This step allows to pause at a certain temperature. It is typically the last 
#' step in an amplification protocol.Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{pauseType$new(temperature)}
#'   
#' @section Fields: \describe{ \item{\code{temperature}}{\link[base]{numeric}.
#'   The temperature in degrees Celsius to maintain during the pause.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
pauseType <- 
  R6Class("pauseType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(temperature) {
              assert_that(is.float(temperature))
              private$.temperature <- temperature
            }
          ),
          private = list(
            .temperature = NULL
          ),
          active = list(
            temperature = function(temperature) {
              if (missing(temperature))
                return(private$.temperature)
              assert_that(is.float(temperature))
              private$.temperature <- temperature
            }
          ))


# loopType ------------------------------------------------------------

#' loopType R6 class.
#' 
#' This step allows to form a loop or to exclude some steps. It allows to jump 
#' to a certain "goto" step for "repeat" times. If the "goto" step is higher 
#' than the step of the loop, "repeat" must be "0".Inherits: 
#' \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{loopType$new(goto, repeat.n)}
#'   
#' @section Fields: \describe{ \item{\code{goto}}{\link[base]{numeric}.  The
#'   step to go to to form the loop.}
#' \item{\code{repeat.n}}{\link[base]{numeric}. Determines how often the loop is 
#'   repeated. The first run through the loop is counted as 0, the last loop is 
#'   "repeat" - 1. The loop is run through exactly "repeat" times.}}
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
loopType <- 
  R6Class("loopType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(goto,
                                  repeat.n) {
              assert_that(is.count(goto))
              assert_that(is.count(repeat.n))
              private$.goto <- goto
              private$.repeat <- repeat.n
            }
          ),
          private = list(
            .goto = NULL,
            .repeat = NULL
          ),
          active = list(
            goto = function(goto) {
              if (missing(goto))
                return(private$.goto)
              assert_that(is.float(goto))
              private$.goto <- goto
            },
            repeat.n = function(repeat.n) {
              if (missing(repeat.n))
                return(private$.repeat)
              assert_that(is.float(repeat.n))
              private$.repeat <- repeat.n
            }
          ))

# measureType ------------------------------------------------------------

#' measureType R6 class.
#' 
#' Can take values:
#' \describe{
#' \item{real time}{}
#' \item{meltcurve}{}
#' }  
#' Inherits: \link{enumType}.
#' 
#' @section Initialization: \code{measureType$new(value)}
#' 
#'   @section Fields: \describe{ 
#' \item{\code{value}}{\link[assertthat]{is.string}.}
#' }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
measureType <- 
  R6Class("measureType",
          # class = FALSE,
          inherit = enumType,
          private = list(
            .levels = c("real time",
                        "meltcurve")
          )
  )

# baseTemperatureType ------------------------------------------------------------

#' baseTemperatureType R6 class.
#' 
#' Parent class for inner usage. Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{baseTemperatureType$new(duration, 
#'   temperatureChange = NULL,  durationChange = NULL, measure = NULL, ramp = 
#'   NULL)}
#'   
#' @section Fields: \describe{ 
#'   \item{\code{duration}}{\link[assertthat]{is.count}. The duration of this
#'   step in } seconds. \item{\code{temperatureChange}}{\link[base]{double}. The
#'   change of the temperature from one cycle to the next: actual temperature
#'   = temperature + (temperatureChange * cycle counter)}
#'   \item{\code{durationChange}}{\link[assertthat]{is.count}. The change of the
#'   duration from one cycle to the next: actual duration = duration +
#'   (durationChange * cycle counter)} \item{\code{measure}}{\link{measureType}.
#'   Indicates to make a measurement and store it as meltcurve or real-time
#'   data.} \item{\code{ramp}}{\link[base]{double}. The allowed temperature
#'   change from one step to the next in degrees Celsius per second. No value
#'   means maximal change rate.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
baseTemperatureType <- 
  R6Class("baseTemperatureType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(duration,
                                  temperatureChange = NULL,
                                  durationChange = NULL,
                                  measure = NULL,
                                  ramp = NULL) {
              assert_that(is.count(duration))
              assert_that(is.opt.double(temperatureChange))
              assert_that(is.opt.count(durationChange))
              assert_that(is.opt.type(measure,
                                      measureType))
              assert_that(is.opt.double(ramp))
              
              private$.duration <- duration
              private$.temperatureChange <- temperatureChange
              private$.durationChange <- durationChange
              private$.measure <- measure
              private$.ramp <- ramp
            }
          ),
          private = list(
            .duration = NULL,
            .temperatureChange = NULL,
            .durationChange = NULL,
            .measure = NULL,
            .ramp = NULL
          ),
          active = list(
            duration = function(duration) {
              if (missing(duration))
                return(private$.duration)
              assert_that(is.count(duration))
              private$.duration <- duration
            },
            temperatureChange = function(temperatureChange) {
              if (missing(temperatureChange))
                return(private$.temperatureChange)
              assert_that(is.opt.double(temperatureChange))
              private$.temperatureChange <- temperatureChange
            },
            durationChange = function(durationChange) {
              if (missing(durationChange))
                return(private$.durationChange)
              assert_that(is.opt.count(durationChange))
              private$.durationChange <- durationChange
            },
            measure = function(measure) {
              if (missing(measure))
                return(private$.measure)
              assert_that(is.opt.type(measure,
                                      measureType))
              private$.measure <- measure
            },
            ramp = function(ramp) {
              if (missing(ramp))
                return(private$.ramp)
              assert_that(is.opt.double(ramp))
              private$.ramp <- ramp
            }
          ))

# temperatureType ------------------------------------------------------------

#' temperatureType R6 class.
#' 
#' This step keeps a constant temperature on the heat block. Inherits: 
#' \link{baseTemperatureType}.
#' 
#' @section Initialization: \code{temperatureType$new(temperature, ...)}
#'   
#' @section Fields: \describe{ \item{\code{temperature}}{\link[base]{double}.
#'   The temperature of the step in  degrees Celsius.}
#' \item{\code{...}}{ Params of parent class \link{baseTemperatureType}.}
#' } 
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
temperatureType <- 
  R6Class("temperatureType",
          # class = FALSE,
          inherit = baseTemperatureType,
          public = list(
            initialize = function(temperature, ...) {
              assert_that(is.float(temperature))
              private$.temperature <- temperature
              super$initialize(...)
            }
          ),
          private = list(
            .temperature = NULL
          ),
          active = list(
            temperature = function(temperature) {
              if (missing(temperature))
                return(private$.temperature)
              assert_that(is.float(temperature))
              private$.temperature <- temperature
            }
          ))

# gradientType ------------------------------------------------------------

#' gradientType R6 class.
#' 
#' This step forms a temperature gradient across the PCR block. Inherits: 
#' \link{baseTemperatureType}.
#' 
#' @section Initialization: \code{gradientType$new(highTemperature, 
#'   lowTemperature, ...)}
#'   
#' @section Fields: \describe{ 
#'   \item{\code{highTemperature}}{\link[base]{double}. The high temperature of
#'   thegradient in degrees Celsius.}
#'   \item{\code{lowTemperature}}{\link[base]{double}. The low temperature of
#'   the gradient in degrees Celsius.}
#' \item{\code{...}}{ Params of parent class \link{baseTemperatureType}. }}
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
gradientType <- 
  R6Class("gradientType",
          # class = FALSE,
          inherit = baseTemperatureType,
          public = list(
            initialize = function(highTemperature,
                                  lowTemperature,
                                  ...) {
              assert_that(is.float(highTemperature))
              assert_that(is.float(lowTemperature))
              private$.highTemperature <- highTemperature
              private$.lowTemperature <- lowTemperature
              super$initialize(...)
            }
          ),
          private = list(
            .highTemperature = NULL,
            .lowTemperature = NULL
          ),
          active = list(
            highTemperature = function(highTemperature) {
              if (missing(highTemperature))
                return(private$.highTemperature)
              assert_that(is.float(highTemperature))
              private$.highTemperature <- highTemperature
            },
            lowTemperature = function(lowTemperature) {
              if (missing(lowTemperature))
                return(private$.lowTemperature)
              assert_that(is.float(lowTemperature))
              private$.lowTemperature <- lowTemperature
            }
          ))


# stepType ------------------------------------------------------------

#' stepType R6 class.
#' 
#' Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{stepType$new(nr, description = NULL, 
#'   temperature = NULL, gradient = NULL, loop = NULL, pause = NULL, lidOpen = 
#'   NULL)}
#'   
#' @section Fields: \describe{ \item{\code{nr}}{\link[assertthat]{is.count}. The
#'   incremental number of the step. First step should be nr = 1 and then
#'   increment each step by + 1. }
#'   \item{\code{description}}{\link[assertthat]{is.string}.} 
#'   \item{\code{temperature}}{\link{temperatureType}.} 
#'   \item{\code{gradient}}{\link{gradientType}.} 
#'   \item{\code{loop}}{\link{loopType}.} \item{\code{pause}}{\link{pauseType}.}
#'   \item{\code{lidOpen}}{\link{lidOpenType}.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
stepType <- 
  R6Class("stepType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(nr,
                                  description = NULL,
                                  temperature = NULL,
                                  gradient = NULL,
                                  loop = NULL,
                                  pause = NULL,
                                  lidOpen = NULL) {
              assert_that(is.count(nr))
              assert_that(is.opt.string(description))
              assert_that(is.opt.type(temperature,
                                      temperatureType))
              assert_that(is.opt.type(gradient,
                                      gradientType))
              assert_that(is.opt.type(loop,
                                      loopType))
              assert_that(is.opt.type(pause,
                                      pauseType))
              assert_that(is.opt.type(lidOpen,
                                      lidOpenType))
              
              private$.nr <- nr
              private$.description <- description
              private$.temperature <- temperature
              private$.gradient <- gradient
              private$.loop <- loop
              private$.pause <- pause
              private$.lidOpen <- lidOpen
            }
          ),
          private = list(
            .nr = NULL,
            .description = NULL,
            .temperature = NULL,
            .gradient = NULL,
            .loop = NULL,
            .pause = NULL,
            .lidOpen = NULL
          ),
          active = list(
            nr = function(nr) {
              if (missing(nr))
                return(private$.nr)
              assert_that(is.count(nr))
              private$.nr <- nr
            },
            description = function(description) {
              if (missing(description))
                return(private$.description)
              assert_that(is.opt.string(description))
              private$.description <- description
            },
            temperature = function(temperature) {
              if (missing(temperature))
                return(private$.temperature)
              assert_that(is.opt.type(temperature,
                                      temperatureType))
              private$.temperature <- temperature
            },
            gradient = function(gradient) {
              if (missing(gradient))
                return(private$.gradient)
              assert_that(is.opt.type(gradient,
                                      gradientType))
              private$.gradient <- gradient
            },
            loop = function(loop) {
              if (missing(loop))
                return(private$.loop)
              assert_that(is.opt.type(loop,
                                      loopType))
              private$.loop <- loop
            },
            pause = function(pause) {
              if (missing(pause))
                return(private$.pause)
              assert_that(is.opt.type(pause,
                                      pauseType))
              private$.pause <- pause
            },
            lidOpen = function(lidOpen) {
              if (missing(lidOpen))
                return(private$.lidOpen)
              assert_that(is.opt.type(lidOpen,
                                      lidOpenType))
              private$.lidOpen <- lidOpen
            }
          ))

# thermalCyclingConditionsType ------------------------------------------------------------

#' thermalCyclingConditionsType R6 class.
#' 
#' A cycling program for PCR or to amplify cDNA. Inherits: \link{rdmlBaseType}.
#' 
#' @section Initialization: \code{thermalCyclingConditionsType$new(id, 
#'   description = NULL, documentation = NULL, lidTemperature = NULL, 
#'   experimenter = NULL, step)}
#'   
#' @section Fields: \describe{ \item{\code{id}}{\link{idType}.} 
#'   \item{\code{description}}{\link[assertthat]{is.string}.} 
#'   \item{\code{documentation}}{\code{list} of \link{idReferencesType}.} 
#'   \item{\code{lidTemperature}}{\link[base]{double}. The temperature in
#'   degrees Celsius of the lid during cycling. }
#'   \item{\code{experimenter}}{\code{list} of \link{idReferencesType}.
#'   Reference to the person who made or uses this protocol. }
#'   \item{\code{step}}{\code{list} of \link{stepType}. The steps a protocol
#'   runs through to amplify DNA.}
#'   }
#'   
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @export
thermalCyclingConditionsType <- 
  R6Class("thermalCyclingConditionsType",
          # class = FALSE,
          inherit = rdmlBaseType,
          public = list(
            initialize = function(id,
                                  description = NULL,
                                  documentation = NULL,
                                  lidTemperature = NULL,
                                  experimenter = NULL,
                                  step) {
              assert_that(is.type(id,
                                  idType))
              assert_that(is.opt.string(description))
              assert_that(is.opt.list.type(documentation,
                                           idReferencesType))
              assert_that(is.opt.double(lidTemperature))
              assert_that(is.opt.list.type(experimenter,
                                           idReferencesType))
              assert_that(is.opt.list.type(step,
                                           stepType))
              
              private$.id <- id
              private$.description <- description
              private$.documentation <- documentation
              private$.lidTemperature <- lidTemperature
              private$.experimenter <- experimenter
              private$.step <- step
              
            }
          ),
          private = list(
            .id = NULL,
            .description = NULL,
            .documentation = NULL,
            .lidTemperature = NULL,
            .experimenter = NULL,
            .step = NULL
          ),
          active = list(
            id = function(id) {
              if (missing(id))
                return(private$.id)
              assert_that(is.type(id, idType))
              private$.id <- id
            },
            description = function(description) {
              if (missing(description))
                return(private$.description)
              assert_that(is.opt.string(description))
              private$.description <- description
            },
            documentation = function(documentation) {
              if (missing(documentation))
                return(private$.documentation)
              assert_that(is.opt.list.type(documentation,
                                           idReferencesType))
              private$.documentation <- documentation
            },
            lidTemperature = function(lidTemperature) {
              if (missing(lidTemperature))
                return(private$.lidTemperature)
              assert_that(is.opt.double(lidTemperature))
              private$.lidTemperature <- lidTemperature
            },
            experimenter = function(experimenter) {
              if (missing(experimenter))
                return(private$.experimenter)
              assert_that(is.opt.list.type(experimenter,
                                           idReferencesType))
              private$.experimenter <- experimenter
            },
            step = function(step) {
              if (missing(step))
                return(private$.step)
              assert_that(is.opt.list.type(step,
                                           idReferencesType))
              private$.step <- step
            }
          ))
