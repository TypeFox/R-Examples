#' A class for querying SPARQL end points
#' 
#' This class allows to run SELECT statement on SPARQL endpoints.
#' The resource parameter is interpreted as SPARQL statement.
#'
#' @seealso \code{\link{xsparql}}, \code{\link{dbpedia}}, \code{\link{enipedia}}, \code{\link{openei}}
#' 
#' @examples
#' getSlots("Xsparql")
#' 
#' @name Xsparql-class
#' @rdname Xsparql-class
#' @exportClass Xsparql
setClass(
    Class="Xsparql", 
    representation=representation(url="character", nspace="character", map.lst="list", stmt="character", resource="character"), 
    contains="Xdata",
    validity=function(object) {
    if((object@nspace != "") && length(object@nspace) %% 2 != 0)
        stop("invalid nspace parameter, character vector of even length or '' expected.")
    }
)

#' Constructor for Xsparql
#'
#' The function \code{xsparql} constructs a Xsparql object.
#'
#' @param resource  the resource the query represents
#' @param url       sparql end point
#' @param statement the SPARQL statement the query stands for
#' @param nspace    character vector with short name / namespace expansions
#' @param clss      class to create, default Xsparql
#' @param ...       named parameters for resources
#'
#' @return a Xsparql object
#' @export
xsparql <- function(resource, url, statement, nspace="", clss="Xsparql", ...) {
    res <- new(clss, resource=resource, url=url, nspace=nspace, stmt=statement, map.lst=list(...))
    return(res)
}

# Internal query method for SPARQL end points
#
# Internal function, use query(xsparql(), ...) instead.
#
# @param url        URL of SPARQL end point
# @param query      SPARQL statement
# @param typeconv   if TRUE (default), converts numbers and dates to R types
# @param verbose    if TRUE, print diagnostic messages. Defaults to getOption("verbose")
#
# @return data.frame object
# @author see SPARQL package
mySPARQL <- function (url, query, typeconv=TRUE, verbose=getOption("verbose")) {
  if (url=="") stop("missing url")
  if (query=="") stop("missing query")
  if(verbose) cat(paste(url, "?query=", URLencode(query), sep = ""),"\n")
    
  js <- RJSONIO::fromJSON(RCurl::getURL(
          paste(url, "?query=", URLencode(query), sep = ""), 
          httpheader = c(Accept = "application/sparql-results+json")
        )
  )
  if(length(js$results$bindings)==0) return(data.frame())
  attrs <- js$head$vars
  empty_row <- as.data.frame(
      matrix(rep(NA, length(attrs)),
      ncol=length(attrs),
      dimnames=list(NULL, attrs)
      )
  )
  
  one_row <- function(l) {
    row <- empty_row
    one_value <- function(name) {
      if (typeconv && "datatype" %in% names(l[[name]]))  
        if (grepl(".*integer$|.*float$|.*double$|.*int$", l[[name]][["datatype"]]))
          row[[name]] <<- as.numeric(l[[name]][["value"]])
        else if(grepl(".*date$", l[[name]][["datatype"]]))
          row[[name]] <<- as.Date(l[[name]][["value"]])
        else 
          row[[name]] <<- l[[name]][["value"]]
      else
        row[[name]] <<- l[[name]][["value"]]
    }
    sapply(intersect(attrs, names(l)), one_value)
    return(row)
  }
  
  res <- Reduce(rbind, lapply(js$results$bindings, one_row))
  rownames(res) <- NULL
  return(data.frame(res))
}

#' @param maxrows      (Xsparql) limit of lines to return (default NULL)
#' @param interactive  (Xsparql) if TRUE, display result in chunks (default FALSE)
#' @param typeconv     (Xsparql) if TRUE (default), convert numbers and dates
#' @rdname query-methods
#' @name query
#' @export
#' @docType methods
#' @aliases query query,Xsparql,character-method
setMethod(
    f="query",
    signature=c(self="Xsparql", resource="character"),
    definition=function(self, resource, maxrows=NULL, interactive=FALSE, typeconv=TRUE, verbose=getOption("verbose"), ...) {
        if(resource==self@resource) {
            if(verbose) cat("building sparql statement..\n")
            mapped <- list()
            arg.lst <- list(...)
            unused <- setdiff(names(arg.lst), names(self@map.lst))
            if(length(unused)>0) warning("unused argument(s) to query: '", paste(unused, collapse="', '"), "'")
            for(n in names(self@map.lst)) {
                arg.fct <- self@map.lst[[n]]
                arg <- arg.lst[[n]]
                if(!is.null(arg)) 
                    mapped[[n]] <- if(is.function(arg.fct)) arg.fct(arg) else as.character(arg)
                else
                    mapped[[n]] <- if(is.function(arg.fct)) arg.fct() else arg.fct
            }
            # browser()
            na.arg <- sapply(mapped, is.na)
            if(any(na.arg)) stop("query with missing required parameters: '", paste(names(which(na.arg)), collapse="', '"), "'")
            stmt <- strsubst(self@stmt, mapped)
        
            if(verbose) cat("querying '", stmt, "'\n")
            prefix <- c()
            if(length(self@nspace) %% 2 != 0) 
                for(i in seq(1,length(self@nspace)-1,2)) 
                    prefix <- c(prefix, paste("PREFIX ", self@nspace[[i]], ": ", self@nspace[[i+1]], sep=""))
            prefix <- paste(prefix, collapse="\n")
            query <- paste(prefix, stmt, sep="")
            if(is.null(maxrows)) {
                d <- mySPARQL(url=self@url, query=query, typeconv=typeconv, verbose=verbose)
            } else {
                offset <- 0
                limit <- maxrows
                d <- NULL
                repeat {
                    chunk <- mySPARQL(url=self@url, query=paste(query, "LIMIT", limit, "OFFSET", offset), typeconv=typeconv, verbose=verbose)
                    d <- if(is.null(d)) chunk$results else rbind(d,chunk$results)
                    offset <- offset + limit
                    if(is.null(chunk$results) || nrow(chunk$results) < limit) break
                    if(interactive) {
                        print(chunk$results)
                        input <- readline("hit <c> to continue, <enter> to exit.. ")
                        if(input=="") break
                    } else if(nrow(d) >= maxrows) break
                }
            }
            return(d)
        }
    }
)

#' @rdname queries-methods
#' @name queries
#' @export
#' @docType methods
#' @aliases queries queries,Xsparql-method
setMethod(
  f="queries",
  signature=c("Xsparql"),
  definition=function(self) c(self@resource, callNextMethod())
)

#' @rdname meta-methods
#' @name meta
#' @export
#' @docType methods
#' @aliases meta meta,Xsparql-method
setMethod(
    f="meta",
    signature="Xsparql",
    definition=function(self) {
        
        prefix <- c()
        if(length(self@nspace) %% 2 != 0) 
            for(i in seq(1,length(self@nspace)-1,2)) 
                prefix <- c(prefix, paste("PREFIX ", self@nspace[[i]], ": ", self@nspace[[i+1]], sep=""))
        prefix <- paste(prefix, collapse="\n")
        
        entity.count <- try(mySPARQL(
            url=self@url, 
            query=paste(prefix, "SELECT COUNT(distinct ?s) AS ?no { ?s a []  }", sep=""), 
            typeconv=TRUE
        )[1,1], silent=TRUE)
        if(inherits(entity.count, "try-error")) {warning("could not sparquery entitiy count data"); entity.count <- data.frame(no=NA)}
        counts <- try(mySPARQL(
            url=self@url, 
            query="SELECT (COUNT(*) AS ?triples) (COUNT(distinct ?p) AS ?predicates) (COUNT(distinct ?s) AS ?subjects) { ?s ?p ?o  }", 
            typeconv=TRUE
        ), silent=TRUE)
        if(inherits(counts, "try-error")) {warning("could not sparquery count data"); counts <- data.frame(triples=NA, predicates=NA, subjects=NA)}
        class.count <- try(mySPARQL(
            url=self@url, 
            query=paste(prefix, "SELECT COUNT(distinct ?o) AS ?no { ?s rdf:type ?o }", sep=""),
            typeconv=TRUE
        )[1,1], silent=TRUE)
        if(inherits(class.count, "try-error")) {warning("could not sparquery class count data"); entity.count <- data.frame(no=NA)}
        object.count <- try(mySPARQL(
            url=self@url, 
            query=paste(prefix, "SELECT (COUNT(DISTINCT ?o ) AS ?no) {  ?s ?p ?o  filter(!isLiteral(?o)) } ", sep=""),
            typeconv=TRUE
        )[1,1], silent=TRUE)
        if(inherits(object.count, "try-error")) {warning("could not sparquery object count data"); object.count <- data.frame(no=NA)}
        props <- try(mySPARQL(
            url=self@url, 
            query=paste(prefix, "SELECT  ?p (COUNT(?s) AS ?subjects) (COUNT(?o) AS ?objects) (COUNT(*) AS ?triples) { ?s ?p ?o } GROUP BY ?p", sep=""),
            typeconv=TRUE
        ), silent=TRUE)
        if(inherits(props, "try-error")) {warning("could not sparquery property data"); props <- data.frame() }

        d <- list(
            counts=c(
                triples=counts[1, "triples"], 
                predicates=counts[1, "predicates"],
                subjects=counts[1, "subjects"],
                entities=entity.count,
                objects=object.count
            ),
            properties=props
        )
        return(d)
    }
)



# getinfo.d.sparql <- function(x, tag, param=NULL, maxrows=100, interactive=TRUE) {
  # if(tag=="classes") {
    # if(is.null(param)) 
      # stmt <- "SELECT DISTINCT ?txt, ?class WHERE {?s a ?class . ?class rdfs:label ?txt .}"
    # else
      # stmt <- sprintf("SELECT DISTINCT ?txt, ?class WHERE {?s a ?class . ?class rdfs:label ?txt . ?txt bif:contains '%s'.}", param)
    # return(query(x,stmt, maxrows, interactive))
  # } else if(tag=="properties") {
    # if(is.null(param)) 
      # stmt <- "SELECT DISTINCT ?property WHERE {?s ?property ?o .}"
    # else
      # stmt <- sprintf("SELECT DISTINCT ?property WHERE {?s a %s; ?property ?o .}", param)
    # return(query(x,stmt, maxrows, interactive))  
  # }
# }
