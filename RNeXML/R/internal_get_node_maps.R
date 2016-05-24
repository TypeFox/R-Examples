# get otus map
#
# @param nexml nexml object 
# @return a list showing the mapping between (internal) otu identifiers and labels (taxonomic names). List is named by the id of the otus block. 
# @details largely for internal use   
get_otu_maps <- function(nexml){
  otus <- as.list(nexml@otus)
  names(otus) <- name_by_id(otus)
  otu_maps <- 
    lapply(otus, function(otus){ # loop over all otus nodes  
    taxon <- sapply(otus@otu, function(otu){ # loop over each otu in the otus set
      if(length(otu@label) > 0) 
        label <- otu@label
      else
        label <- otu@id
      c(otu@id, label)
    })
    out <- taxon[2, ] #label 
    names(out) <- taxon[1, ] #id 
    out
  })
  otu_maps
}

get_char_maps <- function(nexml){
  map <- lapply(nexml@characters, function(characters)
         map_chars_to_label(characters@format))
  names(map) <- name_by_id(nexml@characters)
  map
}

get_state_maps <- function(nexml){
  map <- lapply(nexml@characters, function(characters){
    if(!isEmpty(characters@format@states))
      map_state_to_symbol(characters@format)
    else
     NULL
  })    
  names(map) <- name_by_id(nexml@characters)
  map
}


reverse_map <- function(map){
  out <- NULL
  if(is.list(map)){
    out <- lapply(map, 
                  function(x){
                    out <- names(x)
                    names(out) <- x
                    out})
   } else if(is.character(map)) {
    out <- names(map)
    names(out) <- map
    out
  }
  out   
}


