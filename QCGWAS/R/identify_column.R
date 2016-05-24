identify_column <-
function(std_name, alt_names, header) { return(which(header %in% c(alt_names[ which(alt_names[ ,1]==std_name), 2], std_name))) }
