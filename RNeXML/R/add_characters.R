#################### Write character matices into S4 #####################


#' Add character data to a nexml object
#'
#' @param x character data, in which character traits labels are column names
#'        and taxon labels are row names.  x can be in matrix or data.frame
#'        format. 
#' @param nexml a nexml object, if appending character table to an existing
#'        nexml object.  If ommitted will initiate a new nexml object.  
#' @param append_to_existing_otus logical. If TRUE, will add any new taxa
#' (taxa not matching any existing otus block) to the existing (first)
#' otus block.  Otherwise (default), a new otus block is created, even
#' though it may contain duplicate taxa to those already present.  While
#' FALSE is the safe option, TRUE may be appropriate when building nexml
#' files from scratch with both characters and trees.  
#' @include classes.R
#' @examples
#' library("geiger")
#' data(geospiza)
#' geiger_nex <- add_characters(geospiza$dat)
#' @export 
add_characters <- function(x, 
                           nexml = new("nexml"), 
                           append_to_existing_otus=FALSE){

  # FIXME does it make sense to take a phylo object here as an option?  
  # If so, perhaps don't call the argument 'nexml'.  
  # If not, then we don't really need this conversion. (maybe a type-check instead).   
  nexml <- as(nexml, "nexml")

  ## Check types  & row names ##
  x <- format_characters(x) 
  j <- length(nexml@characters) ## add after any existing character matrices

  nexml <- add_character_nodes(nexml, x)
  for(i in 1:length(x)){

    new_taxa <- rownames(x[[i]]) 
    nexml <- add_otu(nexml, new_taxa, append=append_to_existing_otus)
    ## Add the otus id to the characters node
    otus_id <- nexml@otus[[length(nexml@otus)]]@id
    nexml@characters[[i+j]]@otus <- get_by_id(nexml@otus, otus_id)@id

    nexml <- add_char(nexml, x, i, j)
    nexml <- add_states(nexml, x, i, j)
  }
  for(i in 1:length(x)){
    nexml <- add_rows(nexml, x, i, j)  
  }

  nexml
}

add_character_nodes <- function(nexml, x){
  n <- length(x)
  cs_list <- lapply(1:n, function(i){
         uid <- nexml_id("cs")
         characters <- new("characters", 
             id = uid,
             about = paste0("#", uid))
         if(class(x[[i]][[1]]) == "numeric")  ## Should be numeric but not integer!
           type <- "ContinuousCells"
         else ## Should be integer!
           type <- "StandardCells"
         slot(characters, "xsi:type") <- type  
         characters
      })
  cs_list <- c(nexml@characters, cs_list)
  nexml@characters <-  new("ListOfcharacters", cs_list)
  nexml
}



otu_list <- function(to_add, prefix="ou"){
  lapply(to_add, 
         function(label){
         uid <- nexml_id(prefix)
         new("otu", label=label, id =uid, about = paste0("#", uid))})
}


add_otu <- function(nexml, new_taxa, append=FALSE){

  current_taxa <- get_taxa_list(nexml)

  if(length(current_taxa) == 0) { # No otus exist, create a new node 
    otus <- new_otus_block(nexml, new_taxa)    
    nexml@otus <- new("ListOfotus", c(nexml@otus, otus))

  } else {

    otu_pos <- lapply(current_taxa, 
                      function(current) 
                        match(new_taxa, current))

    if(any(is.na(unlist(otu_pos)))){ # We have missing taxa
      if(append){
        ## append to otus block `otus_id` ##
        otus_id <- 1 # position that matches the id string
        to_add <- new_taxa[sapply(otu_pos, is.na)]  
        nexml@otus[[otus_id]]@otu <- new("ListOfotu", 
                                   c(nexml@otus[[otus_id]]@otu, 
                                     otu_list(to_add, "ou_char"))) ## FIXME hack to make sure new ids are 'unique', 
      } else {
        ## Alternatively, do not append ## 
        otus <- new_otus_block(nexml, new_taxa)    
        nexml@otus <- new("ListOfotus", c(nexml@otus, otus))
      }
    } # else # all taxa matched, so we're all set
  }

  nexml
}

new_otus_block <- function(nexml, to_add){
    id <- nexml_id("os")
    new("ListOfotus", 
        list(new("otus",
                 id = id,
                 about = paste0("#", id),
                 otu = new("ListOfotu", 
                           otu_list(to_add))
                 )))
}




# Turns char names into char nodes 
add_char <- function(nexml, x, i = 1, j = 0){
  char_labels <- colnames(x[[i]])
  char_list <- 
    lapply(char_labels, function(lab){
      id <- nexml_id("cr")
      char <- new("char", 
                  id = id, 
                  about = paste0("#", id),
                  label = lab) 
      })
  nexml@characters[[i+j]]@format@char <- new("ListOfchar", char_list)
  nexml 
}



add_states <- function(nexml, x, i = 1, J = 0){
  # don't ctreate a states node if data is numeric
  if(all(sapply(x[[i]], is.numeric)))
    nexml
  else { 
    nchars <- length(x[[i]])
    char <- nexml@characters[[i+J]]@format@char
    states_list <- lapply(1:nchars, function(j){
      lab <- char[[j]]@label
      lvls <- levels(x[[i]][[lab]])
      id <- nexml_id("ss")
      states <- 
      new("states", 
          id = id,
          about = paste0("#", id),
          state = new("ListOfstate",  
          lapply(lvls, function(lvl){
            new("state", 
                id=nexml_id("s"),
                symbol = as.integer(as.factor(lvl)))
          }))
      )
    })
    nexml@characters[[i+J]]@format@states <- new("ListOfstates", states_list)
    # Add the states's id to char
    for(j in 1:nchars)
      nexml@characters[[i+J]]@format@char[[j]]@states <- states_list[[j]]@id
  nexml
  }
  nexml
}      



## Assumes that otu ids have already been added to the nexml 
add_rows <- function(nexml, x, i = 1, j = 0){

  X <- x[[i]]
  taxa <- rownames(X)
  char_labels <- colnames(X)

  ## get the relevant characters block and otus block
  cs <- nexml@characters[[i+j]]@id   
  os <- nexml@characters[[i+j]]@otus 

  otu_map <- get_otu_maps(nexml)[[os]]
  char_map <- get_char_maps(nexml)[[cs]] 
  state_map <- get_state_maps(nexml)[[cs]]

  reverse_otu_map <- reverse_map(otu_map) 
  reverse_char_map <- reverse_map(char_map) 
  reverse_state_map <- reverse_map(state_map) 

  mat <- 
    new("obsmatrix", 
        row = new("ListOfrow", 
    lapply(taxa, 
      function(taxon){
        id = nexml_id("rw")
        new("row",
            id = id,
            about = paste0("#", id),
            label = taxon,
            otu = reverse_otu_map[taxon],
            cell = new("ListOfcell",  
             lapply(char_labels, function(char){
                    state <- X[taxon,char] # unmapped
                    char_id <- reverse_char_map[[char]]
                    if(!is.null(state_map))
                      state <- reverse_state_map[[char_id]][state]
                    new("cell",
                        char = char_id,
                        state = as.character(state))
                       }))
           )
        }))
    )
  nexml@characters[[i+j]]@matrix <- mat
  nexml
}




## divide matrix into discrete and continuous trait matrices, if necessary
## then write each as separate <characters> nodes: 
## x should now be a list of data.frames of common type
format_characters <- function(x){

  if(is(x, "numeric"))
    x <- as.data.frame(x)

  ## Actually useful conversions ##
  ## Matrices are either all-numeric or all-character class, so no risk of mixed discrete and continous states.  
  if(is(x, "matrix")){
    x <- list(as.data.frame(x))
    
  ## Data.frames can mix discrete and continous states, so we need to seperate them
  } else if(is(x, "data.frame") && dim(x)[2] > 1) {
    x <- split_by_class(x)

  } else if(is(x, "data.frame") && dim(x)[2] == 1) {
    x <- list(x) 

  #### Ugh, this next bit isn't pretty.  Maybe we should  just hope lists are formatted correctly, e.g. come from get_character_list. 

  ## If we're getting a list with matrices, coerce them into data.frames and hope for the best.  
  ## If the list has data.frames, check that each one has consistent class type.  
  ## Otherwise, panic.  
  } else if(is(x, "list")) {
    for(i in 1:length(x)){
      if(is(x[[i]], "matrix"))
        x[[i]] <- as.data.frame(x[[i]])  ## A list of matrices we can make into a list of data.frames...
    }
  ## Someone didn't even try to read the documentation...
  } else { 
    stop("x must be a named numeric, matrix, data.frame, or list thereof")
  }
  ## Let's just hope folks read the documentation and have 
  ## row names as taxa and column names as be character traits.  
  ## Kinda hard to check that for sure?  

  ## return the updated object: a list of data.frames
  x
}

## Helper function for the above, contains the primary functionality
# divide a data.frame into a list of data.frames, in which each has only a unique column class
split_by_class <- function(x){
    col.classes <- sapply(x, class)
    if(all(sapply(col.classes, identical, col.classes[1])))
      x <- list(x)
    else {
      ## split into numerics and non-numerics 
      cts <- unname(which(col.classes=="numeric"))
      discrete <- unname(which(col.classes!="numeric"))
      x <- list(x[cts], x[discrete])
    }

  x
}


