## Transformations.

## The OAI-PMH request issues can either return the "raw" XML results,
## or a suitable transformation thereof, which in turn can be obtained
## by calling oaih_transform() on the raw XML.

## What should the XML results be transformed into?
##
## Earlier versions were based on the following design:
##
## OAI-PMH requests that describe a single "entity" (repository for
## Identify(), record for GetRecord()) give a list.  Requests that
## possibly describe multiple entities (all ListSomething() requests)
## give data frames.
## Elements of entities are mapped as follows.
## Elements which are strings occurring exactly once: character vectors.
## Elements which are strings occurring arbitrarily often: lists of
## character vectors (NULL if missing).
## Elements which are XML "nodes" occurring at most once: lists of XML
## nodes (or NULL if missing).
##
## This complies with the "usual" approach of arranging rectangular,
## case by variables data in data frames (and allows to distinguish the
## XML sequence text-type elements which must occur *exactly* once from
## others).  However, it arranges data by "columns" (variables) which is
## somewhat unnatural when processing data by "rows".  In particular, a
## row from the data frame corresponding to a ListRecords() request would
## be different from the list transformation of the corresponding
## GetRecord() request.  Or equivalently, the data frame cannot simply
## be built from the individual lists via rbind().
##
## We thus provide a design which, in the rectangular cases, treats rows
## and columns symmetrically by arranging data in a "list matrix" (a
## list with a dim attribute).  As matrix subscripting drops dimensions
## when a single row or column is selected, one gets a "simple" list
## (without a dim attribute) in these cases.  Transforming in a list
## context (ListSomething() requests, or lists of same-kind XML nodes)
## gives a list matrix.
##
## The specifics are somewhat tricky.
## One possible guideline for transformation could be the following:
## A list of character vectors can be left as is.
## For a single non-leaf XML node transform the children.
## For a list of different-kind XML nodes, split according to kind, and
## return a list of the transformed sublists.
## For a list of same-kind XML nodes, get their values if they are text
## nodes, and otherwise leave alone.
##
## A possible implementation:
##
## oaih_transform <-
## function(x, rbind = FALSE)
## {
##     xml_value_if_text_node <- function(n)
##         if(length(v <- xmlValue(n, recursive = FALSE))) v else n

##     if(inherits(x, "XMLNode")) {
##         kids <- xmlChildren(x)
##         if(!length(kids))
##             return(xml_value_if_text_node(x))
##         x <- kids
##     } else {
##         if(!is.list(x))
##             stop("Cannot transform")
##         else if(all(sapply(x, is.character)))
##             return(x)
##         else if(!all(sapply(x, inherits, "XMLNode")))
##             stop("Cannot transform")
##     }
##     names(x) <- NULL
##     nms <- sapply(x, xmlName)
##     if(length(unique(nms)) > 1L)
##         tapply(x, nms, sapply, xml_value_if_text_node)
##     else
##         sapply(x, xml_value_if_text_node)
## }
##
## (in fact, the sapply() at the end is not quite correct).
##
## The downside of this "agnostic" approach to transformation is:
## * Possible loss of convenience in cases where it is felt that a
##   different representation is more appropriate.  E.g., for a record
##   we have
##     header metadata? about*
##   where in turn header is
##     identifier datestamp setSpec*
##   (all character) and about can be "anything", and it seems useful to
##   "flatten out" the all-character nodes.  (One might be able to turn
##   this into a general principle, though.)
## * Elements which may be missing (minOccurs = 0 in the DTD) will be
##   missing, which at least makes it inconvenient to combine results.
##
## Hence, instead we use tailor-made explicit transformations for
## certains kinds of nodes, with the underlying idea that transforming a
## list of "same" nodes will result in a list matrix with rows
## corresponding to transforming a single node.  We thus arrange for the
## transformation methods to always work in a list context (i.e., to be
## vectorized) and return list matrices, and the transformation function
## to maintain the context.

oaih_transform <-
function(x)
{
    list_type_node_names <-
        c("ListIdentifiers", "ListMetadataFormats", "ListRecords",
          "ListSets")
    ## Internal nodes would not inherit from XMLNode ...
    xml_node_classes <- c("XMLNode", "XMLAbstractNode")
    if(inherits(x, xml_node_classes)) {
        if(xmlName(x) %in% list_type_node_names)
            .oaih_vtransform(xmlChildren(x))
        else
            .oaih_vtransform(list(x))[1L, ]
    }
    else {
        if(!is.list(x))
            stop("Can only transform lists and XML nodes")
        else if(all(sapply(x, is.character)))
            return(x)
        else if(!all(sapply(x, inherits, xml_node_classes)))
            stop("Can only transform lists of XML nodes or character vectors")
        .oaih_vtransform(x)
    }
}

## Only allow for lists of same-kind nodes for now.

.oaih_vtransform <-
function(x)
{
    if(length(nm <- unique(sapply(x, xmlName))) > 1L)
        stop("Cannot transform mixed-kind node lists")
    names(x) <- NULL
    ## <DM>
    ## registry implemented instead of:
    ## trafo <- sprintf(".oaih_vtransform_%s", nm)
    ## if(exists(trafo, mode = "function"))
    ##    trafo <- get(trafo, mode = "function")
    ## else
    ##    stop(gettextf("Cannot transform node kind '%s'", nm))
    ## </DM>
    if(is.null(trafo <- oaih_transform_methods_db[[nm]]))
        stop(gettextf("Cannot transform node kind '%s'", nm),
             domain = NA)
    trafo(x)
}

oaih_transform_methods_db <- new.env()

oaih_transform_methods_db$Identify <-
function(x)
{
    ## A list of <Identify> nodes.
    ## See http://www.openarchives.org/OAI/2.0/OAI-PMH.xsd for the
    ## definition of identifyType: this
    ## * MUST include one instance of:
    ##     repositoryName baseURL protocolVersion earliestDatestamp
    ##     deletedRecord granularity
    ## * MUST include one or more instances of
    ##     adminEmail
    ## * MAY include multiple instances of the following OPTIONAL
    ##   elements:
    ##     compression description
    ## i.e.,
    ##    repositoryName baseURL protocolVersion earliestDatestamp
    ##    adminEmail+ compression* description*
    ## where description can be anything and everything else is string
    ## type.
    Identify_vars <-
        c("repositoryName", "baseURL", "protocolVersion",
          "earliestDatestamp", "deletedRecord", "granularity",
          "adminEmail", "compression")
    cbind(.get_value_of_any_text_nodes(x, Identify_vars),
          description =
          lapply(x, .get_all_any_nodes, "description"))
}

## .oaih_vtransform_branding <-
## function(x)
## {
##     ## A list of <branding> nodes using a schema for collection branding
##     ## within OAI.
##     ## See http://www.openarchives.org/OAI/2.0/branding.xsd for the
##     ## definition of branding nodes: a sequence with elements
##     ##   collectionIcon? metadataRendering*
##     ## which in turn are sequences with elements
##     ##   url link? title? width? height?
##     ## (the first three string-type, the last two integer) and of base
##     ## URLs with attributes metadataNamespace and mimeType.
##     ## But what can this usefully be transformed into?
##     ## So let's leave it for now ...
## }

oaih_transform_methods_db$dc <-
function(x)
{
    ## Dublin Core:
    ## A list of <oai:dc> nodes.
    ## See http://www.openarchives.org/OAI/2.0/oai_dc.xsd for the
    ## definition of oai_dcType.
    ## Note that each of the 15 string-type variables can occur
    ## arbitrarly often in an oai_dc metadata node.
    dc_vars <-
        c("title", "creator", "subject", "description", "publisher",
          "contributor", "date", "type", "format", "identifier",
          "source", "language", "relation", "coverage", "rights")
    .get_value_of_any_text_nodes(x, dc_vars)
}

oaih_transform_methods_db$eprints <-
function(x)
{
    ## A list of <eprint> nodes, for use in the description section of
    ## an Identify() reply, defined by the e-print community.
    ## See http://www.openarchives.org/OAI/2.0/eprints.xsd for the
    ## definition of eprints:eprintsDescriptionType: a sequence with
    ## elements
    ##   content? metadataPolicy dataPolicy submissionPolicy comment*
    ## with the last a string and the others of TextURLType, consisting
    ## of an URL and text.
    ## Grr, so we must flatten this out ... ugly.
    eprint_textURL_vars <-
        c("content", "metadataPolicy", "dataPolicy", "submissionPolicy")
    out <- lapply(x,
                  function(n)
                  lapply(n[eprint_textURL_vars],
                         function(x)
                         c(.xml_value_if_not_null(x[["URL"]],
                                                  NA_character_),
                           .xml_value_if_not_null(x[["text"]],
                                                  NA_character_)
                           )))
    out <- matrix(as.list(do.call(rbind, lapply(out, unlist))),
                  nrow = length(x), ncol = 8L,
                  dimnames =
                  list(NULL,
                       c(t(outer(eprint_textURL_vars, c("URL", "text"),
                                 paste, sep = ".")))))
    cbind(out, .get_value_of_any_text_nodes(x, "comment"))
}

oaih_transform_methods_db$friends <-
function(x)
{
    ## A list of <friend> nodes, for use in the description section of
    ## an Identify() reply for repositories that list "friends"
    ## (confederate data providers).
    ## See http://www.openarchives.org/OAI/2.0/friends.xsd for the
    ## definition of friends:friendsType: a sequence with elements
    ##   baseURL*
    ## (string-type).
    .get_value_of_any_text_nodes(x, "baseURL")
}

oaih_transform_methods_db$header <-
function(x)
{
    ## A list of <header> nodes.
    ## See http://www.openarchives.org/OAI/2.0/OAI-PMH.xsd for the
    ## definition of headerType: a sequence with elements
    ##   identifier datestamp setSpec*
    ## and optionally a status attribute (all string-type).
    ## <FIXME>
    ## Add support for the status attribute.
    ## </FIXME>
    cbind(.get_value_of_unique_text_nodes(x, c("identifier",
                                               "datestamp")),
          setSpec =
          lapply(x,
                 function(kid)
                 as.character(sapply(kid[names(kid) == "setSpec"],
                                     .xml_value_in_utf8))))
}

oaih_transform_methods_db$metadataFormat <-
function(x)
{
    ## A list of <metadataFormat> nodes.
    ## See http://www.openarchives.org/OAI/2.0/OAI-PMH.xsd for the
    ## definition of metadataFormatType: a sequence with elements
    ##   metadataPrefix schema metadataNamespace
    ## (all string-type).
    metadataFormat_vars <-
        c("metadataPrefix", "schema", "metadataNamespace")
    .get_value_of_unique_text_nodes(x, metadataFormat_vars)
}

oaih_transform_methods_db$`oai-identifier` <-
function(x)
{
    ## A list of <oai-identifier> nodes, for use in the description
    ## section of an Identify() reply.
    ## See http://www.openarchives.org/OAI/2.0/oai-identifier.xsd for
    ## the definition of oai-identifierType: a sequence with elements
    ##   scheme repositoryIdentifier delimiter sampleIdentifier
    ## (all string-type).
    oai_identifier_vars <-
        c("scheme", "repositoryIdentifier", "delimiter",
          "sampleIdentifier")
    .get_value_of_unique_text_nodes(x, oai_identifier_vars)
}

oaih_transform_methods_db$rfc1807 <-
function(x)
{
    ## A list of <rfc1807> metadata nodes.
    ## See http://www.openarchives.org/OAI/rfc1807.xsd for the
    ## definition of fc1807:rfc1807Type: a sequence with the elements in
    ## rfc1807_vars below, where the first three occur exactly once and
    ## the others arbitrarily often (all text-type).
    rfc1807_vars <-
        c("bib-version", "id", "entry", "organization", "title", "type",
          "revision", "withdraw", "author", "corp-author", "contact",
          "date", "pages", "copyright", "handle", "other_access",
          "retrieval", "keyword", "cr-category", "period", "series",
          "monitoring", "funding", "contract", "grant", "language",
          "notes", "abstract")
    .get_value_of_any_text_nodes(x, rfc1807_vars)
}

oaih_transform_methods_db$record <-
function(x)
{
    ## A list of <record> nodes.
    ## See http://www.openarchives.org/OAI/2.0/OAI-PMH.xsd for the
    ## definition of recordType: a sequence with elements
    ##   header metadata? about*
    ## where metadata and about can be "anything".
    cbind(.oaih_vtransform(lapply(x, `[[`, "header")),
          metadata = lapply(x, .get_one_any_node, "metadata"),
          about = lapply(x, .get_all_any_nodes, "about"))
}

oaih_transform_methods_db$set <-
function(x)
{
    ## A list of <set> nodes.
    ## See http://www.openarchives.org/OAI/2.0/OAI-PMH.xsd for the
    ## definition of setType: a sequence with elements
    ##   setSpec setName setDescription*
    ## where setDescription can be "anything".
    cbind(.get_value_of_unique_text_nodes(x, c("setSpec", "setName")),
          setDescription =
          lapply(x, .get_all_any_nodes, "setDescription"))
}

### register functions
## transfuns <- registry()
## transfuns$set_field("Key", type = "character", is_key = TRUE)
## transfuns$set_field("FUN", type = "function")
## transfuns$set_entry("Identify", .oaih_vtransform_Identify)
## transfuns$set_entry("dc", .oaih_vtransform_dc)
## transfuns$set_entry("eprints", .oaih_vtransform_eprints)
## transfuns$set_entry("friends", .oaih_vtransform_friends)
## transfuns$set_entry("header", .oaih_vtransform_header)
## transfuns$set_entry("metadataFormat", .oaih_vtransform_metadataFormat)
## transfuns$set_entry("oai-identifier", `.oaih_vtransform_oai-identifier`)
## transfuns$set_entry("rfc1807", .oaih_vtransform_rfc1807)
## transfuns$set_entry("record", .oaih_vtransform_record)
## transfuns$set_entry("set", .oaih_vtransform_set)

