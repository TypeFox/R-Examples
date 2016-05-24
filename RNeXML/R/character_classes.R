#' @include classes.R

####################################################


setClass("char",
         slots = c(states = "character"),
         contains = "IDTagged")
setMethod("fromNeXML", 
          signature("char", "XMLInternalElementNode"),
          function(obj, from){
            obj <- callNextMethod()
            if(!is.na(xmlAttrs(from)["states"]))
              obj@states <- xmlAttrs(from)["states"]
            obj
          })
setMethod("toNeXML", 
          signature("char", "XMLInternalElementNode"),
          function(object, parent){
            parent <- callNextMethod()
            if(length(object@states) > 0)
              addAttributes(parent, "states" = object@states)
            parent
          })
setAs("char", "XMLInternalNode",
      function(from) toNeXML(from, newXMLNode("char")))
setAs("char", "XMLInternalElementNode",
      function(from) toNeXML(from, newXMLNode("char")))
setAs("XMLInternalElementNode", "char",
      function(from) fromNeXML(new("char"), from))




###############################################

setClass("ListOfrow", slots = c(names="character"), contains="list")
setClass("obsmatrix",
         slots = c(row="ListOfrow"),
         contains = "Annotated")
setMethod("fromNeXML", 
          signature("obsmatrix", "XMLInternalElementNode"),
          function(obj, from){
            obj <- callNextMethod()
            kids <- xmlChildren(from)
            if(length(kids) > 0)
              obj@row <- new("ListOfrow", 
                              lapply(kids[names(kids) == "row"], 
                                     as, "row"))
            obj
          })
setMethod("toNeXML", 
          signature("obsmatrix", "XMLInternalElementNode"),
          function(object, parent){
            parent <- callNextMethod()
            addChildren(parent, kids = object@row)
            parent
          })
setAs("obsmatrix", "XMLInternalNode",
      function(from) toNeXML(from, newXMLNode("matrix")))
setAs("obsmatrix", "XMLInternalElementNode",
      function(from) toNeXML(from, newXMLNode("matrix")))
setAs("XMLInternalElementNode", "obsmatrix",
      function(from) fromNeXML(new("obsmatrix"), from))





######################################################

setClass("ListOfcell", slots = c(names="character"), contains="list")
setClass("ListOfseq", slots = c(names="character"), contains="list")

setClass("row",
         slots = c(cell = "ListOfcell",
                        seq = "ListOfseq"),
         contains = "OptionalTaxonLinked")
setMethod("fromNeXML", 
          signature("row", "XMLInternalElementNode"),
          function(obj, from){
            obj <- callNextMethod()
            kids <- xmlChildren(from)
            if(length(kids) > 0){
              if("cell" %in% names(kids))
              obj@cell <- new("ListOfcell", 
                lapply(kids[names(kids) == "cell"], as, "cell"))
              if("seq" %in% names(kids))
              obj@seq <- new("ListOfseq", 
                lapply(kids[names(kids) == "seq"], as, "seq"))
            }
            obj
          })
setMethod("toNeXML", 
          signature("row", "XMLInternalElementNode"),
          function(object, parent){
            parent <- callNextMethod()
            addChildren(parent, kids = object@cell)
            addChildren(parent, kids = object@seq)
            parent
          })
setAs("row", "XMLInternalNode",
      function(from) toNeXML(from, newXMLNode("row")))
setAs("row", "XMLInternalElementNode",
      function(from) toNeXML(from, newXMLNode("row")))
setAs("XMLInternalElementNode", "row",
      function(from) fromNeXML(new("row"), from))

#######################################################
setClass("ListOfstate", slots = c(names="character"), contains="list")
setClass("ListOfpolymorphic_state_set", slots = c(names="character"), contains="list")
setClass("ListOfuncertain_state_set", slots = c(names="character"), contains="list")

setClass("states",
         slots = c(state="ListOfstate", 
                   polymorphic_state_set="ListOfpolymorphic_state_set",
                   uncertain_state_set="ListOfuncertain_state_set"),
         contains = "IDTagged")
setMethod("fromNeXML", 
          signature("states", "XMLInternalElementNode"),
          function(obj, from){
            obj <- callNextMethod()
            kids <- xmlChildren(from)
            if(length(kids) > 0){
              obj@state <- new("ListOfstate", 
                              lapply(kids[names(kids) == "state"], 
                                     as, "state"))
              obj@polymorphic_state_set <- new("ListOfpolymorphic_state_set", 
                             lapply(kids[names(kids) == "polymorphic_state_set"], 
                                    as, "polymorphic_state_set"))
              obj@uncertain_state_set <- new("ListOfuncertain_state_set", 
                                               lapply(kids[names(kids) == "uncertain_state_set"], 
                                                      as, "uncertain_state_set"))
            }
            obj
          })
setMethod("toNeXML", 
          signature("states", "XMLInternalElementNode"),
          function(object, parent){
            parent <- callNextMethod()
            addChildren(parent, kids = object@state)
            parent
          })
setAs("states", "XMLInternalNode",
      function(from) toNeXML(from, newXMLNode("states")))
setAs("states", "XMLInternalElementNode",
      function(from) toNeXML(from, newXMLNode("states")))
setAs("XMLInternalElementNode", "states",
      function(from) fromNeXML(new("states"), from))




####################################################### 
## technically symbol is positive integer http://nexml.org/doc/schema-1/characters/standard/#StandardToken
setClass("state",
         slots = c(symbol = "integer"), 
         contains = "IDTagged")
setMethod("fromNeXML", 
          signature("state", "XMLInternalElementNode"),
          function(obj, from){
            obj <- callNextMethod()
            obj@symbol <- as.integer(xmlAttrs(from)["symbol"])
            obj
          })
setMethod("toNeXML", 
          signature("state", "XMLInternalElementNode"),
          function(object, parent){
            parent <- callNextMethod()
            addAttributes(parent, "symbol" = object@symbol)
            parent
          })
setAs("state", "XMLInternalNode",
      function(from) toNeXML(from, newXMLNode("state")))
setAs("state", "XMLInternalElementNode",
      function(from) toNeXML(from, newXMLNode("state")))
setAs("XMLInternalElementNode", "state",
      function(from) suppressWarnings(fromNeXML(new("state"), from)))

################################################

setClass("ListOfmember", slots = c(names="character"), contains="list")

setClass("uncertain_state_set", 
         slots = c(member = "ListOfmember"),
         contains="state")
setMethod("fromNeXML", 
          signature("uncertain_state_set", "XMLInternalElementNode"),
          function(obj, from){
            obj <- callNextMethod()
            kids <- xmlChildren(from)
            if(length(kids) > 0)
              obj@member <- new("ListOfmember", 
                              lapply(kids[names(kids) == "member"], 
                                     as, "member"))
            obj
          })
setMethod("toNeXML", 
          signature("uncertain_state_set", "XMLInternalElementNode"),
          function(object, parent){
            parent <- callNextMethod()
            addChildren(parent, kids = object@member)
            parent
          })
setAs("uncertain_state_set", "XMLInternalNode",
      function(from) toNeXML(from, newXMLNode("uncertain_state_set")))
setAs("uncertain_state_set", "XMLInternalElementNode",
      function(from) toNeXML(from, newXMLNode("uncertain_state_set")))
setAs("XMLInternalElementNode", "uncertain_state_set",
      function(from) fromNeXML(new("uncertain_state_set"), from))

################################################

setClass("polymorphic_state_set", contains="uncertain_state_set")
setMethod("fromNeXML", 
          signature("polymorphic_state_set", "XMLInternalElementNode"),
          function(obj, from){
            obj <- callNextMethod()
            kids <- xmlChildren(from)
            if(length(kids) > 0)
              obj@member <- new("ListOfmember", 
                              lapply(kids[names(kids) == "member"], 
                                     as, "member"))
            obj
          })
setMethod("toNeXML", 
          signature("polymorphic_state_set", "XMLInternalElementNode"),
          function(object, parent){
            parent <- callNextMethod()
            addChildren(parent, kids = object@member)
            parent
          })
setAs("polymorphic_state_set", "XMLInternalNode",
      function(from) toNeXML(from, newXMLNode("polymorphic_state_set")))
setAs("polymorphic_state_set", "XMLInternalElementNode",
      function(from) toNeXML(from, newXMLNode("polymorphic_state_set")))
setAs("XMLInternalElementNode", "polymorphic_state_set",
      function(from) fromNeXML(new("polymorphic_state_set"), from))


#####################

setClass("cell",
         slots = c(char="character", 
                        state= "character"),
         contains="Base")
setMethod("fromNeXML", 
          signature("cell", "XMLInternalElementNode"),
          function(obj, from){
            obj <- callNextMethod()
            obj@char <- xmlAttrs(from)["char"]
            obj@state <- xmlAttrs(from)["state"]
            obj
          })
setMethod("toNeXML", 
          signature("cell", "XMLInternalElementNode"),
          function(object, parent){
            parent <- callNextMethod()
            addAttributes(parent, "char" = object@char)
            addAttributes(parent, "state" = object@state)
            parent
          })
setAs("cell", "XMLInternalNode",
      function(from) toNeXML(from, newXMLNode("cell")))
setAs("cell", "XMLInternalElementNode",
      function(from) toNeXML(from, newXMLNode("cell")))
setAs("XMLInternalElementNode", "cell",
      function(from) fromNeXML(new("cell"), from))

#########################

setClass("member", 
         slots = c(state="character"),
         contains="Base")
setMethod("fromNeXML", 
          signature("member", "XMLInternalElementNode"),
          function(obj, from){
            obj <- callNextMethod()
            obj@state <- xmlAttrs(from)["state"]
            obj
          })
setMethod("toNeXML", 
          signature("member", "XMLInternalElementNode"),
          function(object, parent){
            parent <- callNextMethod()
            addAttributes(parent, "state" = object@state)
            parent
          })
setAs("member", "XMLInternalNode",
      function(from) toNeXML(from, newXMLNode("member")))
setAs("member", "XMLInternalElementNode",
      function(from) toNeXML(from, newXMLNode("member")))
setAs("XMLInternalElementNode", "member",
      function(from) fromNeXML(new("member"), from))


########################

setClass("seq", 
         slots = c(seq = "character"),
         contains="Base")
setMethod("fromNeXML", 
          signature("seq", "XMLInternalElementNode"),
          function(obj, from){
            obj <- callNextMethod()
            obj@seq <- xmlValue(from)
            obj
          }
)
setMethod("toNeXML", 
          signature("seq", "XMLInternalElementNode"),
          function(object, parent){
            parent <- callNextMethod()
            addChildren(parent, object@seq)
            parent
          })
setAs("seq", "XMLInternalNode",
      function(from) toNeXML(from, newXMLNode("seq")))
setAs("seq", "XMLInternalElementNode",
      function(from) toNeXML(from, newXMLNode("seq")))
setAs("XMLInternalElementNode", "seq",
      function(from) fromNeXML(new("seq"), from))

#########################################

setClass("ListOfchar", slots = c(names="character"), contains="list")
setClass("ListOfstates", slots = c(names="character"), contains="list")

setClass("format", 
         slots = c(states = "ListOfstates", ## FIXME Should be ListOfstates
                        char = "ListOfchar"),
         contains = "Annotated")
setMethod("fromNeXML", 
          signature("format", "XMLInternalElementNode"),
          function(obj, from){
            obj <- callNextMethod()
            kids <- xmlChildren(from)
            if(length(kids) > 0){
              if("char" %in% names(kids))
                obj@char <- new("ListOfchar", 
                                lapply(kids[names(kids) == "char"], 
                                       as, "char"))
              if("states" %in% names(kids))
                obj@states <- new("ListOfstates", 
                                lapply(kids[names(kids) == "states"], 
                                       as, "states"))
            }
            obj
          })
setMethod("toNeXML", 
          signature("format", "XMLInternalElementNode"),
          function(object, parent){
            parent <- callNextMethod()
            if(!isEmpty(object@char))
              addChildren(parent, kids = object@char)
            if(length(object@states) > 0)
              addChildren(parent, kids = object@states)
            parent
          })
setAs("format", "XMLInternalNode",
      function(from) toNeXML(from, newXMLNode("format")))
setAs("format", "XMLInternalElementNode",
      function(from) toNeXML(from, newXMLNode("format")))
setAs("XMLInternalElementNode", "format",
      function(from) fromNeXML(new("format"), from))





####################################################
setClass("characters",
         slots = c(format = "format",
                        matrix = "obsmatrix"),
        contains = "TaxaLinked")
setMethod("fromNeXML", 
          signature("characters", "XMLInternalElementNode"),
          function(obj, from){
            obj <- callNextMethod()
            obj@format <- as(from[["format"]], "format")
            obj@matrix <- as(from[["matrix"]], "obsmatrix")
            obj
          })
setMethod("toNeXML", 
          signature("characters", "XMLInternalElementNode"),
          function(object, parent){
            parent <- callNextMethod()
            parent <- addChildren(parent, format = object@format)
            parent <- addChildren(parent, matrix = object@matrix)
            parent
          })
setAs("characters", "XMLInternalNode",
      function(from) toNeXML(from, newXMLNode("characters")))
setAs("characters", "XMLInternalElementNode",
      function(from) toNeXML(from, newXMLNode("characters")))
setAs("XMLInternalElementNode", "characters",
      function(from) fromNeXML(new("characters"), from))


