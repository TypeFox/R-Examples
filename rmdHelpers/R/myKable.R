myKable <-
function(x
                    , row.names = NA
                    , boldRowNames = TRUE
                    , boldColNames = TRUE
                    , ...){
  # Function to bold row.names and colnames
  # I still need to add explicit handling for things other than markdown
  if(boldRowNames){
    if(is.na(row.names)){
      # Handle defaults
      if(is.null(row.names(x))){
        # Do nothing, won't print
      } else if(identical(row.names(x), as.character(1:nrow(x)))){
        # Do nothing, won't print
      } else{
        row.names(x) <- paste0("**",row.names(x),"**")
      }
    } else if(!row.names){
      # Do nothing
    } else if(row.names){
      # Handle auto include row.names
      row.names(x) <- paste0("**",row.names(x),"**")
    }
  }
  
  if(boldColNames){
    colnames(x) <- paste0("**",colnames(x),"**")
  }

  
  # Send to kable
  kable(x, row.names = row.names, ...)
}
