

### Attributes: ######################################################

extract_attr <- function(x) {
  UseMethod("extract_attr")
}

extract_attr.osm_parsed <- function(parsed) {
  nodeattr <- extract_attr(parsed[[1]])
  wayattr <- extract_attr(parsed[[2]])
  relationattr <- extract_attr(parsed[[3]])

  list(nodeattr = nodeattr, wayattr = wayattr, relationattr = relationattr)
}

extract_attr.node_parsed<- function(elparsed){
  ret <- lapply(elparsed$elements, xmlAttrs)

  if(length(ret) == 0) {
    return(data.frame(id = numeric(),
                      visible = character(),
                      timestamp = character(),  # How to create a POSIXlt(0)?
                      version = numeric(),
                      changeset = numeric(),
                      user = factor(),
                      uid = factor(),
                      lat = numeric(),
                      lon = numeric()))
  }
  
  r0 <- structure(rep(NA_character_, 9),
                  names = c("id", "visible", "timestamp", "version",
                            "changeset", "user", "uid", "lat", "lon"))
  ret <- lapply(ret,
                function(r) {
                  r2 <- r0
                  r2[names(r)] <- unname(r)
                  r2
                })
  
  # if(any(sapply(ret, length)!=9)){
  #   ret <- as.data.frame(do.call("smartbind", ret))
  # } else{
  ret <- as.data.frame(do.call("rbind", ret))
  # }

  ret$timestamp <- strptime(ret$timestamp, format="%Y-%m-%dT%H:%M:%S")
  ret$lat<- as.numeric(as.character(ret$lat))
  ret$lon<- as.numeric(as.character(ret$lon))
  ret$id<- as.numeric(as.character(ret$id))
  ret$version<- as.numeric(as.character(ret$version))
  ret$uid<- as.factor(as.character(ret$uid))
  ret$user<- as.factor(as.character(ret$user))
  ret$changeset<- as.numeric(as.character(ret$changeset))

  ret
}

extract_attr.way_parsed <- function(elparsed){
  ret <- lapply(elparsed$elements, xmlAttrs)

  if(length(ret)==0) {
    return(data.frame(id = numeric(),
                      visible = character(),
                      timestamp = character(),  # How to create a POSIXlt(0)?
                      version = numeric(),
                      changeset = numeric(),
                      user = factor(),
                      uid = factor()))
  }
  
  r0 <- structure(rep(NA_character_, 7),
                  names = c("id", "visible", "timestamp", "version",
                            "changeset", "user", "uid"))
  ret <- lapply(ret,
                function(r) {
                  r2 <- r0
                  r2[names(r)] <- unname(r)
                  r2
                })
  

  # if(any(sapply(ret, length)!=7)){
  #  ret<-as.data.frame(do.call("smartbind", ret))
  # } else{
  ret<-as.data.frame(do.call("rbind", ret))
  # }

  ret$timestamp <- strptime(ret$timestamp, format="%Y-%m-%dT%H:%M:%S")
  ret$id<- as.numeric(as.character(ret$id))
  ret$version<- as.numeric(as.character(ret$version))
  ret$uid<- as.factor(as.character(ret$uid))
  ret$user<- as.factor(as.character(ret$user))
  ret$changeset<- as.numeric(as.character(ret$changeset))

  ret
}

extract_attr.relation_parsed <- extract_attr.way_parsed




### Data: ############################################################

extract_data <- function(x){
  UseMethod("extract_data")
}

extract_data.osm_parsed <- function(parsed){
  nodedata <- extract_data(parsed[[1]])
  waydata <- extract_data(parsed[[2]])
  relationdata <- extract_data(parsed[[3]])
  list(nodedata=nodedata, waydata=waydata, relationdata=relationdata)
}

extract_data.node_parsed <- function(nparsed){
  values<-nparsed$elements[which(sapply(nparsed$elements,xmlSize)!=0)]
  ##auswahl der nodes MIT daten
  if(length(values)==0)
    return(data.frame(id=numeric(),
                      k=factor(),
                      v=factor()))

  ret <- do.call("rbind", lapply(values, xml2long, "data"))

  ret$id <- as.numeric(as.character(ret$id))
  ret$k <- as.factor(as.character(ret$k))
  ret$v <- as.factor(as.character(ret$v))

  ret
}

## extract_data.way_parsed <- function(wparsed){
##   XMLclone<- lapply(wparsed$elements, xmlClone)
##   XMLclone<- removeKids(XMLclone, "nd")
##   XMLclone<- XMLclone[which(sapply(XMLclone, xmlSize)!=0)]
##   if(length(XMLclone)==0)
##     return(data.frame(id=numeric(),
##                       k=factor(),
##                       v=factor()))
##
##   ret <- do.call("rbind", lapply(XMLclone, xml2long, "data"))
##
##   ret$id <- as.numeric(as.character(ret$id))
##   ret$k <- as.factor(as.character(ret$k))
##   ret$v <- as.factor(as.character(ret$v))
##
##   ret
## }
## NOTE: A very strange xmlClone segfault on Linux and MacOS systems
##       make the above version not usable. We were contacted by
##       Roger Bivand who found that problem. We then contacted
##       Duncan Temple Lang to discuss the problem, but he couldn't
##       reproduce the problems on _any_ of his machines.
extract_data.way_parsed <- function(wparsed){
  clone <- lapply(wparsed$elements, xmlToList)
  clone <- lapply(clone,
                  function(x) {
                    x[names(x) != "nd"]
                  })
  clone <- clone[sapply(clone, length) > 1]  # .attrs element

  if ( length(clone) == 0 ) {
    return(data.frame(id = numeric(),
                      k = factor(),
                      v = factor()))
  }

  ret2 <- lapply(clone, xml2long2)
  ret2 <- do.call(rbind, ret2)

  ret2$id <- as.numeric(ret2$id)
  ret2$k <- as.factor(ret2$k)
  ret2$v <- as.factor(ret2$v)

  attr(ret2$k, "names") <- NULL
  attr(ret2$v, "names") <- NULL

  ret2
}


## extract_data.relation_parsed <- function(rparsed){
##   XMLclone<- lapply(rparsed$elements, xmlClone)
##   XMLclone<- removeKids(XMLclone, "member")
##   XMLclone<- XMLclone[which(sapply(XMLclone,xmlSize)!=0)]
##   if(length(XMLclone)==0)
##     return(data.frame(id=numeric(),
##                       k=factor(),
##                       v=factor()))
##
##   ret <- do.call("rbind", lapply(XMLclone, xml2long, "data"))
##
##   ret$id <- as.numeric(as.character(ret$id))
##   ret$k <- as.factor(as.character(ret$k))
##   ret$v <- as.factor(as.character(ret$v))
##
##   ret
## }
extract_data.relation_parsed <- function(rparsed){
  clone <- lapply(rparsed$elements, xmlToList)
  clone <- lapply(clone,
                  function(x) {
                    x[names(x) != "member"]
                  })
  clone <- clone[sapply(clone, length) > 1]  # .attrs element

  if ( length(clone) == 0 ) {
    return(data.frame(id = numeric(),
                      k = factor(),
                      v = factor()))
  }

  ret2 <- lapply(clone, xml2long2)
  ret2 <- do.call(rbind, ret2)

  ret2$id <- as.numeric(ret2$id)
  ret2$k <- as.factor(ret2$k)
  ret2$v <- as.factor(ret2$v)

  attr(ret2$k, "names") <- NULL
  attr(ret2$v, "names") <- NULL

  ret2
}


xml2long2 <- function(x) {
  attrs <- x$.attrs
  x$.attrs <- NULL

  x <- do.call(rbind, x)
  rownames(x) <- NULL

  ret <- as.data.frame(x, stringsAsFactors = FALSE)
  ret$id <- unname(attrs["id"])

  ret <- ret[c(ncol(ret), seq(length=ncol(ret)-1))]

  ret
}


xml2long <- function(x, dfType){
  size<-xmlSize(x)
  ret<-data.frame(id=character(size))
  ret$id<-rep(xmlGetAttr(x, "id"), each=size)

  if(dfType=="data")
    colnames<-c("k", "v")
  if(dfType=="member")
    colnames<-names(xmlAttrs(xmlChildren(x)[[1]]))

  ret[colnames]<-character(size)
  for(i in 1:length(colnames))
    ret[colnames[i]]<-xmlSApply(x, xmlGetAttr, colnames[i])
  ret
}

removeKids <- function(XML, kidsname){
  ## Returns nodes without the specified children
  lapply(XML, function(x) removeChildren(x, kids=which(names(x)==kidsname)))
}



### Relation refernces: ##############################################

extract_ref <- function(x){
  UseMethod("extract_ref")
}

extract_ref.osm_parsed <- function(parsed){
  wayref <- extract_ref(parsed[[2]])
  relationref <- extract_ref(parsed[[3]])
  list(wayref=wayref, relationref=relationref)
}

## extract_ref.way_parsed <- function(wparsed){
##   XMLclone<- lapply(wparsed$elements, xmlClone)
##   XMLclone<- removeKids(XMLclone, "tag")
##   XMLclone<- XMLclone[which(sapply(XMLclone, xmlSize)!=0)]
##   if(length(XMLclone)==0)
##     return(data.frame(id=numeric(), ref=numeric()))
##   ret <- do.call("rbind", lapply(XMLclone, xml2long, "member"))
##   ret$id <- as.numeric(as.character(ret$id))
##   ret$ref <- as.numeric(as.character(ret$ref))
##   ret
## }
extract_ref.way_parsed <- function(wparsed){
  clone <- lapply(wparsed$elements, xmlToList)
  clone <- lapply(clone,
                  function(x) {
                    x[names(x) != "tag"]
                  })
  clone <- clone[sapply(clone, length) > 1]  # .attrs element

  if ( length(clone) == 0 ) {
    return(data.frame(id = numeric(),
                      ref = numeric()))
  }

  ret2 <- lapply(clone, xml2long2)
  ret2 <- do.call(rbind, ret2)

  ret2$id <- as.numeric(ret2$id)
  ret2$ref <- as.numeric(ret2$ref)

  attr(ret2$ref, "names") <- NULL

  ret2
}



## extract_ref.relation_parsed <- function(rparsed){
##   XMLclone<- lapply(rparsed$elements, xmlClone)
##   XMLclone<- removeKids(XMLclone, "tag")
##   XMLclone<- XMLclone[which(sapply(XMLclone, xmlSize)!=0)]
##   if(length(XMLclone)==0)
##     return(data.frame(id=numeric(), ref=numeric(), role=factor(), type = factor()))
##   ret <- do.call("rbind", lapply(XMLclone, xml2long, "member"))
##   ret$id <- as.numeric(as.character(ret$id))
##   ret$ref <- as.numeric(as.character(ret$ref))
##   ret$role <- as.factor(as.character(ret$role))
##   ret$type <- as.factor(as.character(ret$type))
##   ret
## }
extract_ref.relation_parsed <- function(rparsed){
  clone <- lapply(rparsed$elements, xmlToList)
  clone <- lapply(clone,
                  function(x) {
                    x[names(x) != "tag"]
                  })
  clone <- clone[sapply(clone, length) > 1]  # .attrs element

  if ( length(clone) == 0 ) {
    return(data.frame(id = numeric(),
                      ref = numeric(),
                      role = factor(),
                      type = factor()))
  }

  ret2 <- lapply(clone, xml2long2)
  ret2 <- do.call(rbind, ret2)

  ret2$id <- as.numeric(ret2$id)
  ret2$ref <- as.numeric(ret2$ref)
  ret2$type <- as.factor(ret2$type)
  ret2$role <- as.factor(ret2$role)

  attr(ret2$ref, "names") <- NULL
  attr(ret2$role, "names") <- NULL
  attr(ret2$type, "names") <- NULL

  ret2
}


