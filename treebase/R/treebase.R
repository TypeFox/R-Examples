#' A function to pull in the phyologeny/phylogenies matching a search query
#' 
#' @param input a search query (character string)
#' @param by the kind of search; author, taxon, subject, study, etc
#' (see list of possible search terms, details)
#' @param returns should the fn return the tree or the character matrix?
#' @param exact_match force exact matching for author name, taxon, etc.  
#'   Otherwise does partial matching
#' @param max_trees Upper bound for the number of trees returned, good for
#' keeping possibly large initial queries fast
#' @param branch_lengths logical indicating whether should only return 
#'   trees that have branch lengths. 
#' @param curl the handle to the curl web utility for repeated calls, see
#'  the getCurlHandle() function in RCurl package for details.  
#' @param verbose logical indicating level of progress reporting
#' @param pause1 number of seconds to hesitate between requests
#' @param pause2 number of seconds to hesitate between individual files
#' @param attempts number of attempts to access a particular resource
#' @param only_metadata option to only return metadata about matching trees 
#' which lists study.id, tree.id, kind (gene,species,barcode) type (single, consensus)
#' number of taxa, and possible quality score.  
#' @return either a list of trees (multiphylo) or a list of character matrices
#' @keywords utility
#' @details Choose the search type.  Options are: \itemize{
#' \item{abstract        }{ search terms in the publication abstract}
#' \item{author          }{ match authors in the publication}
#' \item{subject         }{ match subject}
#' \item{doi             }{ the unique object identifier for the publication }
#' \item{ncbi            }{ NCBI identifier number for the taxon}
#' \item{kind.tree       }{ Kind of tree (Gene tree, species tree, barcode tree)  }
#' \item{type.tree       }{ type of tree (Consensus or Single)}
#' \item{ntax            }{ number of taxa in the matrix}
#' \item{quality         }{ A quality score for the tree, if it has been rated.  }
#' \item{study           }{ match words in the title of the study or publication}
#' \item{taxon           }{ taxon scientific name }
#' \item{id.study        }{ TreeBASE study ID}
#' \item{id.tree         }{ TreeBASE's unique tree identifier (Tr.id)}
#' \item{id.taxon        }{ taxon identifier number from TreeBase }
#' \item{tree            }{ The title for the tree}
#' \item{type.matrix     }{ Type of matrix }
#' \item{matrix          }{ Name given the the matrix }
#' \item{id.matrix       }{ TreeBASE's unique matrix identifier}
#' \item{nchar           }{ number of characters in the matrix}
#' } 
#' 
#' The package provides partial support for character matrices provided by TreeBASE.
#' At the time of writing, TreeBASE permits ambiguous DNA characters in these matrices,
#' such as `{CG}` indicating either a C or G, which is not supported by any R interpreter,
#' and thus may lead to errors.  
#'   for a description of all possible search options, see
#'   https://spreadsheets.google.com/pub?key=rL--O7pyhR8FcnnG5-ofAlw. 
#' 
#' @examples \dontrun{
#' ## defaults to return phylogeny 
#' Huelsenbeck <- search_treebase("Huelsenbeck", by="author")
#'
#' ## can ask for character matrices:
#' wingless <- search_treebase("2907", by="id.matrix", returns="matrix")
#'
#' ## Some nexus matrices don't meet read.nexus.data's strict requirements,
#' ## these aren't returned
#' H_matrices <- search_treebase("Huelsenbeck", by="author", returns="matrix")
#'
#' ## Use Booleans in search: and, or, not
#' ## Note that by must identify each entry type if a Boolean is given
#' HR_trees <- search_treebase("Ronquist or Hulesenbeck", by=c("author", "author"))
#'
#' ## We'll often use max_trees in the example so that they run quickly, 
#' ## notice the quotes for species.  
#' dolphins <- search_treebase('"Delphinus"', by="taxon", max_trees=5)
#' ## can do exact matches
#' humans <- search_treebase('"Homo sapiens"', by="taxon", exact_match=TRUE, max_trees=10)
#' ## all trees with 5 taxa
#' five <- search_treebase(5, by="ntax", max_trees = 10)
#' ## These are different, a tree id isn't a Study id.  we report both
#' studies <- search_treebase("2377", by="id.study")
#' tree <- search_treebase("2377", by="id.tree")
#' c("TreeID" = tree$Tr.id, "StudyID" = tree$S.id)
#' ## Only results with branch lengths
#' ## Has to grab all the trees first, then toss out ones without branch_lengths
#'Near <- search_treebase("Near", "author", branch_lengths=TRUE)
#'  }
#' @export
search_treebase <- function(input, by, returns = c("tree", "matrix"),   
                            exact_match = FALSE, max_trees = Inf,
                            branch_lengths = FALSE, curl = getCurlHandle(),
                            verbose = TRUE, pause1 = 2, pause2 = 1, attempts = 3,
                            only_metadata = FALSE){

  nterms <- length(by)
  search_term <- character(nterms)
  section <- character(nterms)

  for(i in 1:nterms){
    search_term[i] <- switch(by[i],
                       abstract="dcterms.abtract",
                       citation="dcterms.bibliographicCitation",
                       author = "dcterms.contributor",
                       subject = "dcterms.subject",
                       id.matrix = "tb.identifier.matrix",
                       id.matrix.tb1 = "tb.identifer.matrix.tb1",
                       ncbi = "tb.identifier.ncbi",
                       id.study = "tb.identifier.study",
                       id.study.tb1 = "tb.identifier.study.tb1",
                       id.taxon = "tb.identifer.taxon",
                       taxon.tb1 = "tb.identifier.taxon.tb1",
                       taxonVariant.tb1 = "tb.identifier.taxonVarient.tb1",
                       id.tree = "tb.identifier.tree",
                       ubio = "tb.identifier.ubio",
                       kind.tree = "tb.kind.tree",
                       nchar = "tb.nchar.matrix",
                       ntax = "tb.ntax.matrix", 
                       quality="tb.quality.tree",
                       matrix = "tb.title.matrix",
                       study = "tb.title.study",
                       taxon = "tb.title.taxon",
                       taxonLabel = "tb.title.taxonLabel",
                       taxonVariant = "tb.title.taxonVariant",
                       tree = "tb.title.tree",
                       type.matrix="tb.type.matrix",
                       type.tree = "tb.type.tree",
                       doi="prism.doi")

  # Section specifies what kind of resource the above id refers to, a tree, matrix, or study
  section[i] <- switch(by[i],
                       abstract= "study",
                       citation= "study",
                       author = "study",
                       subject = "study",
                       id.matrix ="matrix",  
                       id.matrix.tb1 = "matrix",
                       ncbi =  "taxon",
                       id.study = "study",
                       id.study.tb1 = "study",
                       id.taxon = "taxon",
                       taxon.tb1 = "taxon",
                       taxonVariant.tb1 = "taxon",
                       id.tree =  "tree",
                       ubio = "taxon",
                       kind.tree = "tree", 
                       nchar = "matrix",
                       ntax = "matrix",
                       quality = "tree",
                       matrix = "matrix",
                       study="study",
                       taxon = "taxon",
                       taxonLabel =  "taxon",
                       taxonVariant = "taxon", 
                       tree = "tree", 
                       type.matrix= "matrix",
                       type.tree = "tree",
                       doi="study")
  }
  if(!all(section == section[1]))
    stop("Multiple queries must belong to the same section (study/taxon/tree/matrix)")
  search_type <- paste(section[1], "/find?query=", sep="")

  search_term[1] <- paste(search_term[1], "=", sep="") 
  if(nterms > 1){
    for(i in 2:nterms){
      input <- sub("(and|or) ", paste("\\1%20", search_term[i], "=", sep=""), input)
    }
  }

  input <- gsub(" +", "%20\\1", input) # whitespace to html space symbol
  input <- gsub("\"", "%22", input) # html quote code at start
  input <- gsub("'", "%22", input) # html quote code at start
  if(by %in% c("doi")) # list of search types that need to be quoted
    input <- paste("%22", input,"%22", sep="")
  if(exact_match){
    search_term <- gsub("=", "==", search_term) # exact match uses (==) 
  }

  returns <- match.arg(returns)
  schema <- switch(returns, tree = "tree", matrix = "matrix")

# We'll always use rss1 as the machine-readable format 
# could attempt to open a webpage instead with html format to allow manual user search
  format <- "&format=rss1"

  # combine into a search query
  # Should eventually update to allow for multiple query terms with booleans
  query <- paste("http://purl.org/phylo/treebase/phylows/", search_type, 
                 search_term[1], input, format, "&recordSchema=", schema, sep="")
  ## display the constructed query to the user
  message(query)

  if(max_trees == Inf)
    max_trees <- "last()"

  out <- get_nex(query, max_trees = max_trees, returns = returns, curl = curl,
                 pause1 = pause1, pause2 = pause2, attempts = attempts,
                 only_metadata = only_metadata)

  if(schema == "tree" && only_metadata == FALSE){
    out <- drop_nontrees(out)
    if(branch_lengths){
      have <- have_branchlength(out)
      out <- out[have]
    }
#  class(out) <- "multiPhylo"
  }
  out
}

#' drop errors from the search
#' @param tr a list of phylogenetic trees returned by search_treebase
#' @return the list of phylogenetic trees returned successfully
#' @details primarily for the internal use of search_treebase, but may be useful
#' @export
drop_nontrees <- function(tr){
  tt <- tr[sapply(tr, function(x) is(x, "phylo"))]
  message(paste("dropped", length(tr)-length(tt), "objects"))
  tt
}




#' imports phylogenetic trees from treebase. internal function
#' @param query : a phylows formatted search, 
#'     see https://sourceforge.net/apps/mediawiki/treebase/index.php?title=API
#' @param max_trees limits the number of trees returned
#'  should be kept.  
#' @param returns should return the tree object or the matrix (of sequences)
#' @param curl the handle to the curl 
#' @param verbose a logical indicating if output should be printed to screen
#' @param pause1 number of seconds to hesitate between requests
#' @param pause2 number of seconds to hesitate between individual files
#' @param attempts number of attempts to access a particular resource
#' @return A list object containing all the trees matching the search 
#'    (of class phylo)
#' @import XML
#' @import RCurl
#' @import ape
#' @import utils
#' @import methods
#' @keywords internal
get_nex <- function(query, max_trees = "last()", returns = "tree", 
                    curl = getCurlHandle(), verbose = TRUE,
                    pause1 = 1, pause2 = 1, attempts = 5,
                    only_metadata = FALSE){
  n_trees <- 0
  ## Note the need for followlocation -- the actual url just resolves to a page that forwards us on
  page1 <- getURLContent(query, followlocation=TRUE, curl=curl)
  xml_hits <- xmlParse(page1)
  message("Query resolved, looking at each matching resource...")

  resources <- getNodeSet(xml_hits, paste("//rdf:li[position()<= ", max_trees, "]", sep=""))

  ## process some metadata
  metadata <- getNodeSet(xml_hits, "//@rdf:about/./..")
  metadata <- metadata[-1] # first value is for the search
  Studies <- sapply(metadata, function(x) xmlValue(x[["isDefinedBy"]]))
  Trees <- sapply(metadata, function(x) xmlValue(x[["link"]]))
  Study.ids <- gsub(".*TB2:S([1-9]+)+", "\\1", Studies)
  Tree.ids <- gsub(".*Tr([1-9]+)+", "\\1", Trees)
  kind <- sapply(metadata, function(x) xmlValue(x[["kind.tree"]]))
  quality <- sapply(metadata, function(x) xmlValue(x[["quality.tree"]]))
  type <- sapply(metadata, function(x) xmlValue(x[["type.tree"]]))
  ntax <- sapply(metadata, function(x) xmlValue(x[["ntax.tree"]]))

  #all_metadata <- sapply(metadata, xmlToList)

  message(paste(length(resources), "resources found matching query"))


  if(only_metadata){
    out <- vector("list", length=length(Trees))
  } else {
    out <- lapply(Trees, 
      try_recursive, returns=returns, curl=curl, pause1=pause1, 
      pause2=pause2, attempts=attempts)
  }
    out <- lapply(1:length(out), function(i){ 
      out[[i]]$S.id <- Study.ids[i]
      out[[i]]$Tr.id <- Tree.ids[i]
      out[[i]]$type <- type[i]
      out[[i]]$kind <- kind[i]
      out[[i]]$quality <- quality[i]
      out[[i]]$ntax <- ntax[i]
      out[[i]] })

  
  out
}

#' Simple function to identify which trees have branch lengths
#' @param trees a list of phylogenetic trees (ape/phylo format)
#' @return logical string indicating which have branch length data
#' @export
have_branchlength <- function(trees){
 sapply(trees, function(x) !is.null(x[["edge.length"]]))
}




## an internal function which descends through the pages to get the nexus resources
dig <- function(tree_url, returns="tree", curl=getCurlHandle(), pause1=1, pause2=1){
# Get the URL to the actual resource on that page 
  #thepage <- xmlAttrs(x, "rdf:resource")

  ## being patient will let the server get the resource ready
  Sys.sleep(pause1)
  target <- getURLContent(tree_url, followlocation=TRUE, curl=curl)
  seconddoc <- xmlParse(target) ## This fails if we rush

  message("Looking for nexus files...")

  ### Here we sometimes get 304 errors, and want to try again
  node <- xpathApply(seconddoc, "//x:item[x:title='Nexus file']", 
            namespaces=c(x="http://purl.org/rss/1.0/"),
            function(x){
              if(is.list(x))
                x = x[[1]]

              ## being patient will let the server get the resource ready
              Sys.sleep(pause2)
              con <- url(xmlValue(x[["link"]]))
              if(returns=="tree"){
                nex <- read.nexus(con) ## fails if we rush
                message("Tree read in successfully")
              } else if (returns=="matrix"){
                 nex <- read.nexus.data(xmlValue(x[["link"]]))
              }
            nex
            })
  node <- node[[1]] # will fail if above has errored
  node
  # One tree per seconddoc 
}


# Helper function to make multiple trys of dig, increasing the patience timing
try_thrice <- function(x,returns, curl, pause1, pause2, attempts){
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  retry <- function(e){
    message("Query failed, attempting a second call")
    tryCatch(dig(x,returns, curl, pause1+2, pause2+2), 
      error = retry2)
  }

  retry2 <- function(e){
    message("Query failed, attempting a third call")
    tryCatch(dig(x,returns, curl, pause1+5, pause2+5), 
      error = function(e){
        message("Resource cannot be reached")
        "Retry failed"})
  }
  withCallingHandlers(
       tryCatch(dig(x, returns, curl, pause1, pause2),
       error = retry))
}


# Helper function to make a specified number of attempts to access a resource

#' @keywords internal
try_recursive <- function(x,returns, curl, pause1, pause2, attempts=3){
  try <- 1
  while(try <= attempts){
    message(paste("Attempting try", try))
    out <- try(dig(x,returns, curl, pause1, pause2))
    try <- try + 1
    if(!is(out, "try-error"))
      try <- attempts+1
  }
  out
}


