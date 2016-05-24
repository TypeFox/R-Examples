# This is hardly graphics-system-neutral
# e.g., it uses the grid notion of what constitutes a viewport
# But anyway ...

#################
# General Device Stuff
#################

setClass("graphicsDevice",
         representation(name="character",
                        width="numeric",
                        height="numeric"))

setClass("graphicsParams",
         representation(pars="list"))

setGeneric("inchToDevX",
           function(x, device) {
             standardGeneric("inchToDevX")
           })

setGeneric("inchToDevY",
           function(x, device) {
             standardGeneric("inchToDevY")
           })

setGeneric("devArrow",
           function(arrow, gp, device) {
             standardGeneric("devArrow")
           })

setGeneric("devStartElement",
           function(element, gp, device) {
             standardGeneric("devStartElement")
           })

setGeneric("devEndElement",
           function(name, device) {
             standardGeneric("devEndElement")
           })

setGeneric("devTextNode",
           function(text, device) {
             standardGeneric("devTextNode")
           })

setGeneric("devStartClip",
           function(clip, gp, device) {
             standardGeneric("devStartClip")
           })

setGeneric("devStartClipPath",
           function(clippath, gp, device) {
             standardGeneric("devStartClipPath")
           })

setGeneric("devEndClipPath",
           function(clippath, gp, device) {
             standardGeneric("devEndClipPath")
           })

setGeneric("devStartClipPathGroup",
           function(clippath, gp, device) {
             standardGeneric("devStartClipPathGroup")
           })

setGeneric("devStartMask",
           function(mask, gp, device) {
             standardGeneric("devStartMask")
           })

setGeneric("devEndMask",
           function(mask, gp, device) {
             standardGeneric("devEndMask")
           })

setGeneric("devStartMaskGroup",
           function(mask, gp, device) {
             standardGeneric("devStartMaskGroup")
           })

setGeneric("devStartGroup",
           function(group, gp, device) {
             standardGeneric("devStartGroup")
           })

setGeneric("devEndGroup",
           function(name, vp, device) {
             standardGeneric("devEndGroup")
           })

setGeneric("devLines",
           function(lines, gp, device) {
             standardGeneric("devLines")
           })

setGeneric("devPolygon",
           function(polygon, gp, device) {
             standardGeneric("devPolygon")
           })

setGeneric("devPath",
           function(path, gp, device) {
             standardGeneric("devPath")
           })

setGeneric("devRaster",
           function(raster, gp, device) {
             standardGeneric("devRaster")
           })

setGeneric("devRect",
           function(rect, gp, device) {
             standardGeneric("devRect")
           })

setGeneric("devText",
           function(text, gp, device) {
             standardGeneric("devText")
           })

setGeneric("devCircle",
           function(circle, gp, device) {
             standardGeneric("devCircle")
           })

setGeneric("devClose",
           function(device) {
             standardGeneric("devClose")
           })

setGeneric("devStartSymbol",
          function(pch, device) {
            standardGeneric("devStartSymbol")
          })

setGeneric("devPoint",
          function(pch, device) {
            standardGeneric("devPoint")
          })

setGeneric("devEndSymbol",
           function(device) {
             standardGeneric("devEndSymbol")
           })

setGeneric("devUseSymbol",
           function(point, gp, device) {
             standardGeneric("devUseSymbol")
           })

