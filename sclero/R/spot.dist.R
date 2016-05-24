#' @title Distance of sampling spots from margin along a measurement axis
#' 
#' @description Scale the location of sampling spots to visible growth bands, project the sampling spots along a measurement axis (\code{main}) and calculate the distance from margin. Useful for locating high-resolution LA-ICP-MS or SIMS sample spots in samples with non-linear growth bands.
#' 
#' @param rawDist a \code{\link[=convert.ijdata]{rawDist}} object for which the alignment should be done.
#' @param sample.name An optional parameter over-riding the sample name. If \code{NULL} (default), the sample name will be passed from the previous steps (\code{\link{read.ijdata}}, \code{\link{convert.ijdata}})
#' @param run.ae a logical indicating whether to run averaging error estimation, if \link[=assign.size]{rawDist contains spot size information}. Defaults to TRUE.
#' @param use.centroids a logical indicating whether to use centroids of spot.owins instead of spots, if \link[=assign.size]{rawDist contains spot size information}. Defaults to TRUE.
#' @details The alignment information with sample spot numbers is stored as a sublist called \code{output} and can be extracted to a data.frame. Otherwise the object behaves like any list in R. Relevant data can be subsetted as needed. Detailed data containing information of the alignment process is stored in a sublist called \code{det.dat}.
#' 
#' The function can either project growth lines on the distance (main) axis or use the crossing points between growth lines and the main axis. These two types of the main axis can be used for different applications. The main axis type is automatically selected by the following criteria:
#'\describe{
#'\item{\strong{along}}{Approperiate for samples with cut-off growth lines such as bivalve margin cross-sections and tree, sediment or ice-cores. This option is selected by placing the measurement axis such that \strong{it does not cross any of the marked growth lines}. The location of each growth line is projected along the measurement axis from the beginning of the growth line (the point where you started marking the growth line in ImageJ).}
#'\item{\strong{cross}}{Approperiate for approximately round cross-sections: samples where the growth lines continue through the entire width of the sample (such as tree, coral or calcareous algae cross-sections and umbo-regions of bivalves). This type is selected by making the \strong{main axis to cross each individual marked growth line}. The location of each growth line along the main axis is considered as a crossing point.}
#'}
#'
#'These criteria are set due to the need of defining a location for each marked growth line along the distance (main) axis. The choice is rigid, to simplify calculations, and to avoid bias in results by allowing two different methods for growth line locations. The easiest way to test which \code{type} suits a particular sample best is to save two sets of ImageJ zip files by moving the measurement axis.
#' @return Returns a list of class \code{spotDist} containing information of the aligned sample spots and the digitized representation of the shell cross-section, which was already included in the \code{\link[=convert.ijdata]{rawDist}} object. 
#' @author Mikko Vihtakari
#' @seealso \code{\link{read.ijdata}} for reading zip files containing ImageJ ROIs.
#' 
#' \code{\link{order.ijdata}} for ordering and subsetting \code{read.ijdata} output.
#' 
#' \code{\link{convert.ijdata}} for converting the coordinate information to \link[spatstat]{spatstat} point patterns. 
#' 
#' \code{\link{plot.spotDist}} for plotting.
#' 
#' \code{\link{print.spotDist}} for printing.
#' 
#' @examples data(shellspots)
#' shell_map <- convert.ijdata(shellspots)
#' x <- spot.dist(shell_map)
#' plot(x) 
#' @import spatstat
#' @export

spot.dist <- function(rawDist, sample.name = NULL, run.ae = TRUE, use.centroids = TRUE){

##These are test parameters. Remove when ready
#rawDist <- convert.ijdata(shellspots); alignment = "closest"; sample.name = NULL; run.ae = TRUE; use.centroids = TRUE

#####################################
### Section I. Helper functions #####
#####################################
# Function 1. Cuts main axis (dist.lines) into gaps and generates marks
seg.fun <- function(z) {
zz <- z[match(dist.lines$line[dist.lines$line %in% z$marks], z$marks),]
marks <- paste(marks(zz)[1:npoints(zz)-1], marks(zz)[2:npoints(zz)], sep = "-")
xx <- coords(zz)$x
yy <- coords(zz)$y
psp(x0 = xx[1:length(xx)-1], x1 = xx[-1], y0 = yy[1:length(yy)-1], y1 = yy[-1], window = window, marks = marks)}

# Function 2. Calculates spot distance along main axis

spot.dist.calculation <- function(spots, gbs, lines, dist.main, l.1st) {
  
## Part 1. Calculate shortest distance (c.dist) from a hole the closest subannual line (close1)

dat <- list()
for(l in 1:length(spots)){
d <- data.frame(spot = spots[[l]]$marks)
for(j in 1:nlines){
temp <- nncross(spots[[l]], gbs[gbs$marks == lines[j]])
colnames(temp) <- c(lines[j], "which")
d <- cbind(d, temp[1])}

tp <- d[!colnames(d) %in% colnames(d)[1]]

d <- cbind(d, close1 = apply(tp, 1, function(x) names(which.min(x))))
d <- cbind(d, c.dist = apply(tp, 1, function(x) (min(x))))

d <- d[!colnames(d) %in% lines]
dat <- c(dat, list(d))}
names(dat) <- names(spots)
d <- dat

## Part 2. Calculate if spots are occuring after or before (from l.1st) the closest subannual line

# Calculate shortest distance from each hole to the 1st growth line (l.1st, when the clam died, margin). Called start here because photographs should be marked against direction of growth.

tmp <- lapply(spots, function(x) nncross(x, l.1st)[,1]) 
d <- mapply(`[<-`, d, 'dist.1st', value = tmp, SIMPLIFY = FALSE)

# Calculate the closest point along the closest growth line ($close1) for each hole. Then calculate the distance from this point to the 1st growth line (l.1st). This distance is called dist.cross. If the logic behind this is hard to grasp, try plotting 'tmp' after project2segment part.

tmp.l <- lapply(seq_along(spots), function(k) {
  tmp <- c()
  for(j in 1:npoints(spots[[k]])){
    tmp <- superimpose(tmp, project2segment(spots[[k]][j], gbs[gbs$marks == d[[k]]$close1[j]])$Xproj)}
    nncross(tmp, l.1st)[,1]})

d <- mapply(`[<-`, d, 'dist.cross', value = tmp.l, SIMPLIFY = FALSE)

# If $dist.1st > $dist.cross, the sample spot occurs after the closest growth line ($close1). If the opposite, then the sample spot occurs before $close1.

tmp <- lapply(d, function(x) factor(ifelse(x$dist.1st >= x$dist.cross, "aft", "bef"))) 
d <- mapply(`[<-`, d, 'location', value = tmp, SIMPLIFY = FALSE)

dat <- d

## Part 3. Calculate the gap between which growth bands each sampling spot is located.

dat <- lapply(seq_along(dat), function(i){
  tmp <- unlist(ifelse(dat[[i]]$location  == "aft", 
    sapply(seq_along(dat[[i]]$location), function(k) {
      grep(paste("^", as.character(dat[[i]]$close1)[k], "-", sep = ""), dist.main$gap, value = TRUE)}),
    sapply(seq_along(dat[[i]]$location), function(k) {
      grep(paste("-", as.character(dat[[i]]$close1)[k], "$", sep = ""), dist.main$gap, value = TRUE)})))
  dat[[i]]$gap <- tmp
  dat[[i]]})

names(dat) <- names(spots)

## Part 4. Shortest distance to the 2nd growth band in the gap the spot is located

# Locate the 2nd closest subannual line
dat <- lapply(dat, function(i){
      i$close2 <- sapply(seq_along(i$close1), function(k) {
      tmp <- unlist(strsplit(i$gap[k], "-"))
      tmp[!tmp %in% i$close1[k]]})
    i})

#Calculate the distance from each hole to the 2nd closest growth band

dat <- lapply(seq_along(dat), function(i){
  dat[[i]]$c2.dist <- sapply(seq_along(dat[[i]]$close2), function(k) 
    nncross(spots[[i]][spots[[i]]$marks==dat[[i]]$spot[k]], gbs[gbs$marks==dat[[i]]$close2[k]])$dist); dat[[i]]})

names(dat) <- names(spots)

## Part 5. Distance ratio to the closest and 2nd closest growth band

dat <- lapply(dat, function(i){
  i$dist.ratio <- ifelse(i$location == "aft", i$c.dist/(i$c.dist + i$c2.dist), i$c2.dist/(i$c.dist + i$c2.dist))
  i})

## Part 6. Compile dat and dist.main datasets 

dat <- lapply(dat, function(i) merge(i, dist.main[-1], by = "gap", all.x = TRUE, sort = FALSE))

## Part 7. Calculate distance of spots along the main axis

tmp <- lapply(dat, function(i) i$dist0 + i$gap.l.main*i$dist.ratio)
dat <- mapply(`[<-`, dat, 'spot.dist.type1', value = tmp, SIMPLIFY = FALSE)

## Part 8. Reorganize dat

colord <- c("spot", "gap", "location", "close1", "c.dist", "close2", "c2.dist", "dist.1st", "dist.cross", "gap.l.main", "dist0", "dist.ratio", "spot.dist.type1")

dat <- lapply(dat, function(k) k[colord])
return(dat)}

###########################################
### Section II. Obligatory calculations ###
###########################################
## Part 1. Define parameters 
## Spots and holes are used synonymously due to historical coding reasons.
x <- rawDist

if(!is.null(sample.name)) x$sample.name <- sample.name
gbs <- x$gbs
main <- x$main
start <- x$start.main
end <- x$end.main
lines <- unique(gbs$marks)
nlines <- length(lines)
spots <- x$spots
window <- x$window
if(exists("spot.area", where = x) & run.ae) spot.area <- x$spot.area

## Part 1.2 Test whether the object is 'along' or 'cross' type (see the tutorial)
test <- crossing.psp(main, superimpose(gbs))
if(test$n == 0){type <- 'along'} else {
  if(test$n == length(unique(gbs$marks))){type <- 'cross'} else {
    stop(paste0("Number of growth line and main axis crossing points is neither 0 or ", length(unique(gbs$marks)), ". Cannot define the aligment type. See ?spot.dist"))
  }
}

## Part 2. Calculate the distance of growth lines on the main axis

if(type == "along"){
tmp <- endpoints.psp(gbs[gbs$marks == lines[1]])
tmp <- superimpose(tmp[1], tmp[npoints(tmp)])
pnt <- tmp[apply(nncross(tmp, main)[1],2, function(x) which.min(x))]
marks(pnt) <- lines[1]

for(j in 2:nlines){
tmp <- endpoints.psp(gbs[gbs$marks == lines[j]])
tmp <- superimpose(tmp[1], tmp[npoints(tmp)])
pnt2 <- tmp[apply(nncross(tmp, main)[1],2, function(x) which.min(x))]
marks(pnt2) <- lines[j]
pnt <- superimpose(pnt, pnt2)}

int <- project2segment(pnt, main)$Xproj}

if(type=="cross"){
int <- ppp()
marks(int) <- as.character(c())
for(j in 1:nlines){
tmp <- crossing.psp(main, gbs[gbs$marks == lines[j]]) 
marks(tmp) <- lines[j]
int <- superimpose(int, tmp)}}

tmp <- t(crossdist(int, start)) 
colnames(tmp) <- int$marks

# The origo for the main axis
p.1st <- int[int$marks == apply(tmp, 1, function(x) names(which.min(x)))]
# First line, i.e. margin, i.e. the line when the animal died (margin)
l.1st <- gbs[gbs$marks == apply(tmp, 1, function(x) names(which.min(x)))] 
# Tip of the margin
p.margin.tip <- endpoints.psp(l.1st)[1]
# Last line, i.e. the line where the animal started growing (holes sequence numbering is against the direction of growth)
l.last <- gbs[gbs$marks == apply(tmp, 1, function(x) names(which.max(x)))]
# The last growth band projection point along the main axis
p.Last <- int[int$marks == apply(tmp, 1, function(x) names(which.max(x)))]

# Compile to a data.frame
dist.lines <- data.frame(line = int$marks, dist0 = crossdist(int, p.1st))
dist.lines <- dist.lines[order(dist.lines$dist0),]

## Part 3. Calculate the gaps growth lines form along the main axis using the function from Section I.

tmp <- seg.fun(int)
dist.main <- data.frame(gap.l.main = lengths.psp(tmp), gap = marks(tmp), line = gsub("\\-.*", "", as.character(marks(tmp))))

# Add distance along the main axis (dist0)

dist.main <- merge(dist.main, dist.lines, by = "line", all = TRUE, sort = FALSE)
dist.main <- dist.main[c("line", "gap", "gap.l.main", "dist0")]

## Part 4. Align the spots along main 

if(exists("spot.area") & use.centroids){
  dat <- spot.dist.calculation(spots = spot.area$centroids, gbs = gbs, lines = lines, dist.main = dist.main, l.1st = l.1st)
} else {
dat <- spot.dist.calculation(spots = spots, gbs = gbs, lines = lines, dist.main = dist.main, l.1st = l.1st)}

#########################################
### Section III Optional calculations ###
#########################################
## Part 1. Averaging error.
if(exists("spot.area")){

## Part 1.1 Create a psp object from owins for distance calculations
SP <- lapply(spot.area$spot.dat, function(k) {
  tmp <- k$spot.owins
  tmp3 <- lapply(seq_along(tmp), function(j) {
    tmp2 <- edges(tmp[[j]], window = x$window)
    marks(tmp2) <- factor(k$spot[j])
    tmp2})
  names(tmp3) <- k$spot
  tmp3})

## Part 1.2 Marks along main axis
mark.list <- unique(marks(x$gbs))

## Part 1.3 Find crossing growth lines
crossing.lines <- lapply(SP, function(k) {
  lapply(k, function(j) {
    unlist(lapply(mark.list, function(i) if(npoints(crossing.psp(j, x$gbs[marks(x$gbs) %in% i])) != 0) i))
    })
  })

crossing.lines.main <- lapply(crossing.lines, function(k) {
    lapply(k, function(j) dist.main[dist.main$line %in% j,])})

## Part 1.4 Find first and last crossing line.

start.crossing.line <- lapply(crossing.lines.main, function(j) {
  lapply(j, function(k) as.character(k[which.min(k$dist0),]$line))})
end.crossing.line <- lapply(crossing.lines.main, function(j) {
  lapply(j, function(k) as.character(k[which.max(k$dist0),]$line))})

## Part 1.5 Find the lines outside spot areas on both sides
start.out.line  <- lapply(seq_along(start.crossing.line), function(k){
    out <- lapply(seq_along(start.crossing.line[[k]]), function(i) {
      if(length(start.crossing.line[[k]][[i]]) == 0) {
        tmp <- names(start.crossing.line[[k]][i])
        tmp2 <- dat[[k]]
        unlist(strsplit(tmp2[tmp2$spot == tmp,]$gap, "-"))[1]} else {
        as.character(dist.main$line[which(dist.main$line %in% start.crossing.line[[k]][[i]])-1])}})
    names(out) <- names(start.crossing.line[[k]])
    out})

names(start.out.line) <- names(start.crossing.line)

end.out.line <- lapply(seq_along(end.crossing.line), function(k){
    out <- lapply(seq_along(end.crossing.line[[k]]), function(i) {
      if(length(end.crossing.line[[k]][[i]]) == 0) {
        tmp <- names(end.crossing.line[[k]][i])
        tmp2 <- dat[[k]]
        unlist(strsplit(tmp2[tmp2$spot == tmp,]$gap, "-"))[2]} else {
        as.character(dist.main$line[which(dist.main$line %in% end.crossing.line[[k]][[i]])+1])}})
    names(out) <- names(end.crossing.line[[k]])
    out})

names(end.out.line) <- names(end.crossing.line)

## Part 1.6 Find start points for each spot

start.point <- lapply(seq_along(SP), function(k) {
    out <- lapply(seq_along(SP[[k]]), function(i) {
    obj <- as.ppp(SP[[k]][[i]])
    obj[which.min(nncross(obj, gbs[marks(gbs) %in% start.out.line[[k]][[i]]])$dist)]})
    names(out) <- names(SP[[k]])
    out})
  
names(start.point) <- names(SP)

start.points <- lapply(start.point, function(k) {
  tmp <- c()
  for(i in seq_along(k)) {
    tmp <- superimpose(tmp, k[[i]])}
  tmp})

## Part 1.7 Find end points for each spot

end.point <- lapply(seq_along(SP), function(k) {
    out <- lapply(seq_along(SP[[k]]), function(i) {
    obj <- as.ppp(SP[[k]][[i]])
    obj[which.min(nncross(obj, gbs[marks(gbs) %in% end.out.line[[k]][[i]]])$dist)]})
    names(out) <- names(SP[[k]])
    out})
  
names(end.point) <- names(SP)

end.points <- lapply(end.point, function(k) {
  tmp <- c()
  for(i in seq_along(k)) {
    tmp <- superimpose(tmp, k[[i]])}
  tmp})

## Part 1.8 Calculate the cross distance for each point

start.point.dat <- spot.dist.calculation(spots = start.points, gbs = gbs, lines = lines, dist.main = dist.main, l.1st = l.1st)
end.point.dat <- spot.dist.calculation(spots = end.points, gbs = gbs, lines = lines, dist.main = dist.main, l.1st = l.1st)

area.dat <- list(centroids = spot.area$centroids, start.points = start.points, end.points = end.points, spot.dat = spot.area$spot.dat)
}

## ADD: traverse type, measurement from tip, measurement from closest point along margin.

##############################
### Section IV Print data  ###
##############################

## The main function is ready now. The rest is cleaning up and setting parameters for the object
# Short output
if(exists("spot.area")){
  output <- lapply(seq_along(dat), function(i) {
    cols <- c("spot", "spot.dist.type1")
    tmp <- dat[[i]][cols]
    colnames(tmp) <- c("spot", "dist")
    tmp2 <- start.point.dat[[i]][cols]
    colnames(tmp2) <- c("spot", "dist.min")
    tmp3 <- end.point.dat[[i]][cols]
    colnames(tmp3) <- c("spot", "dist.max")
    tmp <- merge(tmp, tmp2, by = c("spot"), all = TRUE, sort = FALSE)
    tmp <- merge(tmp, tmp3, by = c("spot"), all = TRUE, sort = FALSE)
    tmp <- tmp[with(tmp, order(spot)),]
    tmp$spot.area <- area.dat$spot.dat[[i]]$spot.area
    tmp$spot.diameter <- area.dat$spot.dat[[i]]$spot.diameter
    return(tmp)})
  names(output) <- names(dat)
  } else {
    output <- lapply(dat, function(k) {
    cols <- c("spot", "spot.dist.type1")
    tmp <- k[cols]
    colnames(tmp) <- c("spot", "dist")
    tmp})  
  }

if(exists("spot.area")){
  lapply(seq_along(output), function(i) {
    test <- which(output[[i]]$dist < output[[i]]$dist.min)
    if(length(test) != 0) warning(paste0(" dist.min > dist in ", names(output[i]), " spot ", output[[i]]$spot[test], ". Check growth lines."))
    test <- which(output[[i]]$dist > output[[i]]$dist.max)
    if(length(test) != 0) warning(paste0(" dist.max < dist in ", names(output[i]), " spot ", output[[i]]$spot[test], ". Check growth lines."))})
}
  
x$main.type <- type
#x$alignment.type <- alignment
x$gb.projections <- int
x$gb.start <- p.1st
x$gb.end <- p.Last
#x$line.cross <- line.cross
#x$cross.segs <- cross.segs
x$gb.projections.dist <- dist.main
x$mid.point.type <- ifelse(exists("spot.area") & use.centroids, "centroids", "spots")
if(exists("spot.area")) x$spot.area <- area.dat
if(exists("spot.area")) {
  x$det.data$start.points <- start.point.dat
  x$det.data$mid.points <- dat
  x$det.data$end.points <- end.point.dat} else {
  x$det.dat$mid.points <- dat
  }
x$output <- output
class(x) <- "spotDist"
  
return(x)}
