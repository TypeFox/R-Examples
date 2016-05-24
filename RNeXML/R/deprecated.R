## Some of these methods are still in use by simmap and various add_ functions.
## Ideally all would be replaced by the equivalent get_level() version



## Conversions between matrix and characters node

#' Extract the character matrix
#'
#' @param nexml nexml object (e.g. from read.nexml)
#' @param rownames_as_col option to return character matrix rownames
#'  (with taxon ids) as it's own column in the data.frame. Default is FALSE
#'  for compatibility with geiger and similar packages.
#' @return the list of taxa
#' @examples
#' comp_analysis <- system.file("examples", "comp_analysis.xml", package="RNeXML")
#' nex <- nexml_read(comp_analysis)
#' get_characters_list(nex)
#' @export
get_characters_list <- function(nexml, rownames_as_col=FALSE){
  # extract mapping between otus and taxon ids 
  maps <- get_otu_maps(nexml)
  # loop over all character matrices 
  out <- lapply(nexml@characters, function(characters){
    
    # Make numeric data class numeric, discrete data class discrete
    type <- slot(characters, 'xsi:type') # check without namespace?
    if (grepl("ContinuousCells", type)) {
      dat <- extract_character_matrix(characters@matrix)
      dat <- otu_to_label(dat, maps[[characters@otus]])
      dat <- character_to_label(dat, characters@format)
      for(i in length(dat)) ## FIXME something more elegant, no?
        dat[[i]] <- as.numeric(dat[[i]])
    } else if (grepl("StandardCells", type)) {
      dat <- extract_character_matrix(characters@matrix)
      dat <- state_to_symbol(dat, characters@format)
      dat <- otu_to_label(dat, maps[[characters@otus]])
      dat <- character_to_label(dat, characters@format)
      for(i in length(dat))
        dat[[i]] <- factor(dat[[i]])
    } else {
      dat <- NULL 
    }
    if(rownames_as_col){
      dat <- cbind(taxa = rownames(dat), dat)
      rownames(dat) <- NULL
    }
    dat
  })
  # name the character matrices by their labels,
  # if available, otherwise, by id.    
  names(out) <- name_by_id_or_label(nexml@characters) 
  out
} 


# for lists only 
# identical_rownames <- function(x) all(sapply(lapply(x, rownames), identical, rownames(x[[1]])))


####  Subroutines (not exported)  ########### 


### The subroutines of "get_characters_list# function
## Fixme these could be adapted to use the get_*_maps functions


otu_to_label <- function(dat, otu_map){
  rownames(dat) <- otu_map[rownames(dat)]
  dat
}

character_to_label <- function(dat, format){ 
  ## Compute the mapping
  map <- map_chars_to_label(format)
  ## replace colnames with matching labels
  colnames(dat) <- map[colnames(dat)]
  dat
}

state_to_symbol <- function(dat, format){
  if(!isEmpty(format@states)){
    map_by_char <- map_state_to_symbol(format)
    for(n in names(dat)){
      symbol <- map_by_char[[n]] 
      dat[[n]] <- symbol[dat[[n]]]
    }
    dat 
  } else {
    dat ## Nothing to do if we don't have a states list
  }
}



### Subroutine for characters_to_label, 
###  Also subroutine for get_char_map . 

map_chars_to_label <- function(format){
  map <- sapply(format@char, function(char){
    if(length(char@label) > 0)
      label <- char@label
    else
      label <- char@id
    c(char@id, label)
  })
  out <- map[2,]
  names(out) <- map[1,]
  out
}




# Subroutine of the get_state_maps and state_to_symbol functions
#
# For each character, find the matching `states` set.
# For that set, map each state id to the state symbol 
map_state_to_symbol <- function(format){
  # loop over characters
  map <- lapply(format@char, function(char){ 
    # name the list with elements as `states` sets by their ids 
    states <- format@states 
    ids <- sapply(states, function(states) states@id)
    names(states) <- ids
    # Get the relevant states set matching the current character 
    map_states_to_symbols( states[[char@states]] )
  })  
  names(map) <- name_by_id(format@char)
  map
}
## Subroutine of the map_state_to_symbol function above  
map_states_to_symbols <- function(states){
  map <- sapply(states@state,
                function(state){
                  if(length(state@symbol) > 0)
                    symbol <- state@symbol
                  else
                    symbol <- state@id
                  c(state@id, symbol)
                })
  out <- map[2, ]
  names(out) <- map[1, ]
  out
}





#' @import reshape2
extract_character_matrix <- function(matrix){ 
  otu <- sapply(matrix@row, function(row) row@otu)
  # charnames <- unname(sapply(matrix@row[[1]]@cell, function(cell) cell@char))
  names(matrix@row) <- otu
  mat <- lapply(matrix@row, function(row){
    names(row@cell) <- unname(sapply(row@cell, function(b) b@char))
    # names(row@cell) <- charnames
    lapply(row@cell, function(cell) cell@state)
  })
  mat <- melt(mat)
  colnames(mat) <- c("state", "character", "otu")
  mat <- dcast(mat, otu ~ character, value.var = "state")
  
  # Move otus into rownames and drop the column 
  rownames(mat) <- mat[["otu"]]
  # mat <- mat[-1]
  mat 
}
