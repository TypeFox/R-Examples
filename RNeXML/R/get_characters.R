#' Get character data.frame from nexml
#' 
#' @param nex a nexml object 
#' @param rownames_as_col option to return character matrix rownames (with taxon ids) as it's own column in the 
#' data.frame. Default is FALSE for compatibility with geiger and similar packages.
#' @param otu_id logical, default FALSE. return a column with the 
#'  otu id (for joining with otu metadata, etc)
#' @param otus_id logical, default FALSE. return a column with the 
#'  otus block id (for joining with otu metadata, etc)
#' @return the character matrix as a data.frame
#' @details RNeXML will attempt to return the matrix using the NeXML taxon (otu) labels to name the rows
#'  and the NeXML char labels to name the traits (columns).  If these are unavailable or not unique, the NeXML
#'  id values for the otus or traits will be used instead.
#' @importFrom tidyr spread
#' @importFrom dplyr left_join select_
#' @importFrom stringr str_replace
#' @importFrom stats setNames
#' @export
#' @examples
#' \dontrun{
#' # A simple example with a discrete and a continous trait
#' f <- system.file("examples", "comp_analysis.xml", package="RNeXML")
#' nex <- read.nexml(f)
#' get_characters(nex)
#' 
#' # A more complex example -- currently ignores sequence-type characters
#' f <- system.file("examples", "characters.xml", package="RNeXML")
#' nex <- read.nexml(f)
#' get_characters(nex)
#' }
get_characters <- function(nex, rownames_as_col=FALSE, otu_id = FALSE, otus_id = FALSE){
  
  drop = lazyeval::interp(~-matches(x), x="about|xsi.type|format")
  
  otus <- get_level(nex, "otus/otu") %>% 
    select_(drop) %>%
    optional_labels(id_col = "otu")
  
  char <- get_level(nex, "characters/format/char") %>% 
    select_(drop) %>%
    optional_labels(id_col = "char")
  
  ## Rows have otu information
  rows <- get_level(nex, "characters/matrix/row") %>% 
    dplyr::select_(.dots = c("otu", "row"))
  cells <- get_level(nex, "characters/matrix/row/cell") %>% 
    dplyr::select_(.dots = c("char", "state", "row")) %>% 
    dplyr::left_join(rows, by = "row")
    
  characters <- get_level(nex, "characters")
  
  ## States, including polymorphic states (or uncertain states)
  states <- get_level(nex, "characters/format/states/state") 
  
  ## Include polymorphic and uncertain states
  polymorph <- get_level(nex, "characters/format/states/polymorphic_state_set") 
  uncertain <- get_level(nex, "characters/format/states/uncertain_state_set") 
  if(dim(polymorph)[1] > 0)
    states <- dplyr::bind_rows(states, polymorph)
  if(dim(uncertain)[1] > 0)
    states <- dplyr::bind_rows(states, uncertain)
  states <- dplyr::select_(states, drop)


  if(dim(states)[1] > 0) 
    cells <- cells %>% 
    dplyr::left_join(states, by = c("state"))  %>% 
    dplyr::select_(.dots = c("char", "symbol", "otu", "state"))
  
  ## Join the matrices.  Note that we select unique column names after each join to avoid collisions
  cells %>% 
    dplyr::left_join(char, by = c("char")) %>% 
    dplyr::rename_(.dots = c("trait" = "label")) %>% 
    dplyr::left_join(otus, by = c("otu")) %>%
    dplyr::rename_(.dots = c("taxa" = "label")) %>%
    na_symbol_to_state() %>% 
    dplyr::select_(.dots = c("taxa", "symbol", "trait", "otu", "otus")) %>% 
    tidyr::spread("trait", "symbol") ->
    out
  

  ## Identify the class of each column and reset it appropriately
  cellclass <- function(x){
    x %>%
    stringr::str_replace(".*ContinuousCells", "numeric") %>%
    stringr::str_replace(".*StandardCells", "integer")
  }
  
  type <-
    get_level(nex, "characters/matrix/row/cell")  %>% 
    dplyr::select_(drop) %>%
    dplyr::left_join(characters, by = "characters") %>%
    dplyr::select_(.dots = c( "xsi.type", "char", "characters")) %>%
    dplyr::left_join(char, by = c("char")) %>%
    dplyr::select_(.dots = c("label", "xsi.type")) %>% 
    dplyr::distinct() %>%
    dplyr::mutate_(.dots = setNames(list(~cellclass(xsi.type)), "class"))
  
  for(i in dim(type)[1])
    class(out[[type$label[i]]]) <- type$class[i]
  
  
  ## drop unwanted columns if requested (default)
  if(!otu_id){
    out <- dplyr::select_(out, quote(-otu))
  }
  if(!otus_id){
    out <- dplyr::select_(out, quote(-otus))
  }
  if(!rownames_as_col){
    taxa <- out$taxa
    out <- dplyr::select_(out, quote(-taxa))
    out <- as.data.frame(out)  
    rownames(out) <- taxa
    out
  }
  
  out
}

## If 'label' column is missing, create it from 'id' column
## if label exists but has missing or non-unique values, also use ids instead
optional_labels <- function(df, id_col = "id"){
  who <- names(df)
  if(! "label" %in% who)
    df$label <- df[[id_col]]
  if(length(unique(df$label)) < length(df$label))
    df$label <- df[[id_col]]
  df
}

## Continuous traits have the values in "state" column, whereas 
## for discrete states we return the value of the "symbol" column 
na_symbol_to_state <- function(df){
  if(is.null(df$symbol))
    df$symbol <- NA
  df$symbol[is.na(df$symbol)] <- suppressWarnings(as.numeric(df$state[is.na(df$symbol)]))
  df
  }
