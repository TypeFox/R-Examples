build.profile <-
function(table, Profile=NULL, name=NULL, depth=NULL, col=NULL, skel_dim=NULL, skel_ab=NULL, type=NULL, root_ab=NULL, root_dim=NULL,
                          orientation=NULL) {
  required.names <- c('Profile','name','depth','skel_dim', 'skel_ab', 'col', 'type', 'root_ab', 'root_dim', 'orientation')
  output.table <- data.frame(matrix(ncol=length(required.names), nrow=dim(table)[1]))
  names(output.table) <- required.names
  names <- names(table)
  used.columns <- NULL
  for (a in 1:length(required.names)) {
    actual.name <- required.names[a]
    ##if the name
    if (length(grep(actual.name, names))==1) {
      output.table[,actual.name] <- table[,actual.name]
      column <- which(names(table)==actual.name)
      used.columns <- c(used.columns, column)
      next()
    }
    else {
      print(paste('looking for a vector with', actual.name))
      if (is.null(get(actual.name))) {
        print(toupper(paste(actual.name, 'is missing!')))
        next()
      } else  {
        output.table[,actual.name] <- get(actual.name)
        responses <- NULL
        for (i in 1:dim(table)[2]) {
          response <- ifelse(length(which(as.character(output.table[,actual.name])==as.character(table[,i])))==0, 'NO', 'YES')
          responses <- c(responses, response)
        }
        column <- which(responses=='YES')
        used.columns <- c(used.columns, column)
        print(paste(actual.name, 'found'))
      }
    }
  }
  remaining <- as.data.frame(table[,-used.columns])
  names(remaining) <- names(table)[-used.columns]
  output.table <- cbind(output.table, remaining)
  class(output.table) <- 'profile.data.frame'
  return(output.table)
}
