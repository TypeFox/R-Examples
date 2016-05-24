### Method to import a SmartPLS XML  model specification file: .splsm
# Uses: 'path', 'order', 'innerW', 'initM1', 'block'
read.splsm <-
function(file=character(), order=c("generic", "alphabetical")){
  if(require(XML)==FALSE) stop("Package 'XML' is required.")
  order <- match.arg(order)
  model <- xmlInternalTreeParse(file)

  # latents
  name1 <- unlist(xpathApply(model, "/model/latent", xmlGetAttr, "name"))
  id1  <- as.numeric(unlist(xpathApply(model, "/model/latent", xmlGetAttr, "id")))
  xPos1 <- as.numeric(unlist(xpathApply(model, "/model/latent", xmlGetAttr, "x")))
  yPos1 <- as.numeric(unlist(xpathApply(model, "/model/latent", xmlGetAttr, "y")))
  # manifests
  name2 <- unlist(xpathApply(model, "/model/manifest", xmlGetAttr, "datacolumn"))
  id2  <- as.numeric(unlist(xpathApply(model, "/model/manifest", xmlGetAttr, "id")))
  xPos2 <- as.numeric(unlist(xpathApply(model, "/model/manifest", xmlGetAttr, "x")))
  yPos2 <- as.numeric(unlist(xpathApply(model, "/model/manifest", xmlGetAttr, "y")))

  # all variables with id and postition as built in SmartPLS
  variables <- data.frame(name=c(name1, name2), id=c(id1, id2),
                          xPos=c(xPos1, xPos2), yPos=c(yPos1, yPos2))

  sources <-  unlist(xpathApply(model, "/model/connection", xmlGetAttr, "source"))
  targets <-  unlist(xpathApply(model, "/model/connection", xmlGetAttr, "target"))

  connections <- data.frame(sourceID=as.numeric(sources), targetID=as.numeric(targets))

  # structural model IDs
  smID <- with(connections, connections[sourceID %in% id1 & targetID %in% id1,])

  # pathes with names only for the structural model
  strucmod <-  path(variables, smID)

  # measurement model IDs
  mmID <- with(connections, connections[(sourceID %in% id1 & targetID %in% id2)
                                          |(sourceID %in% id2 & targetID %in% id1),])
  # pathes with names only for the measurement model
  measuremod <-  path(variables, mmID)

  # all the pathes with names instead of IDs
  ret <- path(variables, connections)

  # Adjacency matrix D for the structural model
  D <- innerW(strucmod=strucmod, latent=name1)

  # Ordering of LVs
  if (order=="generic"){
    tmp <- reorder(D)
    latent <- tmp$chain
    strucmod <- tmp$strucmod
  }
  if (order=="alphabetical"){
    latent <- sort(name1)
  }

  # Arranging the rows and columns according to the order of the LVs
  D <- D[latent, latent]

  # build blocks of manifest variables (including 'measurement mode')
  manifest <- sort(name2)
  blocks <- block(latent, manifest, measuremod=measuremod)

  # Ordering of MVs and measuremod
  MVs <- NULL
  mm <- NULL
  for(i in names(blocks)){
    MVs <- append(MVs, blocks[[i]])
    if(attr(blocks[[i]], "mode") == "A"){
      mm <- rbind(mm, (cbind(i, blocks[[i]])))
    }
    if(attr(blocks[[i]], "mode") == "B"){
      mm <- rbind(mm, (cbind(blocks[[i]], i)))
    }
  }
  dimnames(mm) <- dimnames(measuremod)
  measuremod <- mm

  # Result
  result <- list()
  result$connectionIDs <- connections
  result$variables <- variables
  result$latent <- latent
  result$manifest <- MVs
  result$path <- ret
  result$strucmod <- strucmod
  result$measuremod <- measuremod
  result$D <- D
  result$M <- initM1(model=result)
  result$blocks <- blocks
  result$order <- order
  class(result) <- c("splsm", "plsm")
  return(result)
}
