
lastobs <- function(x, ...)  UseMethod("lastobs")

firstobs <- function(x, ...) UseMethod("firstobs")

lastobs.default <- function(x,...){
  ux<-unique(x)
  m <-match(ux,x)
  sort(sapply(ux, function(i) {max(which(i==x))}))
}

firstobs.default <- function(x,...){
  ux<-unique(x)
  m <-match(ux,x)
  sort(sapply(ux, function(i) {min(which(i==x))}))
}

lastobs.formula <- function(formula, data=parent.frame(), ...){
  rhs<-gsub(" +","",strsplit(paste(formula[2]),"\\+")[[1]][1])
  lastobs(data[,rhs])
}


firstobs.formula <- function(formula, data=parent.frame(), ...){
   rhs<-gsub(" +","",strsplit(paste(formula[2]),"\\+")[[1]][1])
   firstobs(data[,rhs])
}



# lastobs.formula <- function(x, data=parent.frame(), ...){
#    args<-list(...)
#    cl <- match('data.frame', sapply(args,class))   
#    if (is.na(cl))
#       stop("A dataframe must be given\n")
#    data  <- args[[cl]]
#    mcall <- match.call(expand.dots = FALSE)
#    ff  <- as.formula(eval.parent(mcall[[2]]))
#    rhs <- paste(ff)[2]
#    lastobs(data[,rhs])
# }


# firstobs.formula <- function(x, data=parent.frame(), ...){
#    args<-list(...)
#    cl <- match('data.frame', sapply(args,class))   
#    if (is.na(cl))
#       stop("A dataframe must be given\n")
#    data  <- args[[cl]]
#    mcall <- match.call(expand.dots = FALSE)
#    ff  <- as.formula(eval.parent(mcall[[2]]))
#    rhs <- paste(ff)[2]
#    firstobs(data[,rhs])
# }
