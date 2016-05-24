library(testthat)

context("Test rcanvec-package examples that take too long to check")

test_that("No errors on manpage examples", {
  #bbox functions make it easy to manipulate bounding boxes
  wolfville <- searchbbox("wolfville ns") #requires {ggmap} to geocode query
  wolfvillezoomedout <- zoom(wolfville, 0.5)
  
  #easy plotting with canvec.qplot()
  canvec.qplot(bbox=searchbbox("wolfville ns"))
  
  #download canvec or canvec+ data. 250k references use canvec+ (large amounts of data)
  #and 50k references use canvec data (older but distributed in smaller chunks).
  canvec.download(nts('21h1'))
  
  #load data
  buildings <- canvec.load(nts("21h1"), "building")
  lakes <- canvec.load(nts("21h1"), "waterbody")
  rivers <- canvec.load(nts('21h1'), "river")
  roads <- canvec.load(nts('21h1'), "road")
  contours <- canvec.load(nts('21h1'), "contour")
  
  #plot data
  sp::plot(lakes, col="lightblue", border="lightblue")
  sp::plot(rivers, add=TRUE, col="lightblue")
  sp::plot(buildings, add=TRUE, pch=".")
  
  #zoomed in
  sp::plot(lakes, col="lightblue", border="lightblue", 
       xlim=c(-64.4,-64.35), ylim=c(45.05,45.1))
  sp::plot(contours, add=TRUE, col="brown", lwd=0.2)
  sp::plot(rivers, add=TRUE, col="lightblue")
  sp::plot(buildings, add=TRUE, pch=".")
  sp::plot(roads, add=TRUE, lwd=0.5)
  
  
  #equivalent syntax in canvec.qplot()
  canvec.qplot(nts("21h1"), layers=c("waterbody", "contour", "river", "road"))
  canvec.qplot(bbox=makebbox(45.1, -64.35, 45.05, -64.4), 
               layers=c("waterbody", "contour", "river", "building", "road"))
  
  #method returns plot data argument so data does not need to be loaded each time. 
  #this will not work when changing nts sheets.
  plotdata <- canvec.qplot(nts("21h1"), layers=c("waterbody", "contour", "river"))
  plotdata <- canvec.qplot(bbox=makebbox(45.1, -64.35, 45.05, -64.4), 
                           layers=c("waterbody", "contour", "river"),
                           data=plotdata)
  
  #easy exporting with human readable names
  canvec.export(nts("21h01"), "exportdata", layerids=c("road", "river"))
  
  #cleanup
  canvec.cleanup(all=TRUE)
  unlink("exportdata", recursive=TRUE)
  
})