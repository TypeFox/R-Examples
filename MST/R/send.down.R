send.down <-
function(data, tree){

  #Converts a character into logical (e.g. "3 <= 8" equals TRUE)  
  charToLog<-function(char){
    eval(parse(text=char))
  }

  data$node <- 1  #Initialize node to 1
  cut.point <- as.vector(tree$cut)
  split.var <- as.numeric(as.vector(tree$var))
  operator <- tree$operator
  for(i in 1:nrow(tree)){
    in.node <- (data$node)==(tree$node[i])
    if(!is.na(split.var[i])){
      var.split <- data[,split.var[i]]
      cut <- cut.point[i]
      if(operator[i] %in% c("<=",">")){
        cut1 <- as.numeric(cut)
        splitNode<-as.logical(sapply(paste(var.split[in.node],operator[i],cut1),charToLog))
        data$node[in.node]<-paste0(data$node[in.node],(!splitNode)+1)
      }
      else {
        var.split <- as.character(var.split)
        cut1 <- unlist(strsplit(as.character(cut), split=","))
        if(operator[i]=="in"){splitNode <- is.element(var.split[in.node],cut1)
        } else {splitNode <- !is.element(var.split[in.node],cut1)}
        data$node[in.node]<-paste0(data$node[in.node],(!splitNode)+1)
      }
    }
  }
  return(data)
}
