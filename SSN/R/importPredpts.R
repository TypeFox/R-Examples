importPredpts <-
function(target, predpts, obj.type) {

    old_wd <- getwd()
    on.exit(setwd(old_wd))

    if (obj.type == "glm"){
      setwd(target$ssn@path)
      count <- 0

      if (length(target$ssn@predpoints@ID) > 0) {
        for (m in 1:length(target$ssn@predpoints@ID)) {
          if (target$ssn@predpoints@ID[m] == predpts) {
            pred.num <- m
            count <- count + 1}}}

      if (count > 0){
        stop("GLM object already contains predpoints named ", predpts)}}

    if (obj.type == "ssn") {
      setwd(target@path)
      count <- 0

      if (length(target@predpoints@ID) > 0) {
        for (m in 1:length(target@predpoints@ID)) {
          if (target@predpoints@ID[m] == predpts) {
            pred.num <- m
            count <- count + 1}}}

      if (count > 0) {
        stop("SSN already contains predpoints named ", predpts)}}


    if (file.exists(paste(predpts,".shp",sep = ""))==0) {
      stop(paste(predpts,".shp data is missing from ", target@path, " folder",sep = ""))}

    predpoints <- readShapeSpatial(predpts)
    rownames(predpoints@data) <- predpoints@data[,"pid"]
    rownames(predpoints@coords) <- predpoints@data[,"pid"]

    predpoints@data$locID <- as.factor(predpoints@data$locID)

    if (getinfo.shape(predpts)[[2]] != 1){
      stop(paste(predpts,".shp does not have point geometry", sep = ""))}

    network.point.coords <- data.frame(predpoints@data[,"netID"], predpoints@data[,"rid"], predpoints@data[,"upDist"])
    colnames(network.point.coords)<-c("NetworkID", "SegmentID", "DistanceUpstream")
    network.point.coords <- as.data.frame(network.point.coords)
    row.names(network.point.coords) <- row.names(predpoints@data)
	attributes(network.point.coords)$locID <- as.numeric(levels(predpoints@data$locID))[predpoints@data$locID]

    network.point.coords[,1] <- as.factor(network.point.coords[,1])
    network.point.coords[,2] <- as.factor(network.point.coords[,2])

    # Create SSNPoint object for prediction sites
    pp <- new("SSNPoint",
      network.point.coords = network.point.coords,
      point.coords = predpoints@coords,
      point.data = predpoints@data,
      points.bbox = predpoints@bbox)


    if (obj.type == "ssn") {
      pred.num <- length(target@predpoints@SSNPoints)
      target@predpoints@SSNPoints[[pred.num + 1]]<-pp
      target@predpoints@ID[[pred.num + 1]]<-predpts}
    if (obj.type == "glm") {
      pred.num <- length(target$ssn@predpoints@SSNPoints)
      target$ssn.object@predpoints@SSNPoints[[pred.num + 1]]<-pp
      target$ssn.object@predpoints@ID[[pred.num + 1]]<-predpts}

    target
}

