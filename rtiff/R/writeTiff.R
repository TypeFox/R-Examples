"writeTiff" <-
function(pixmap, fn) {
      if(class(pixmap) == "pixmapRGB"){
        .Call("writeTiff", pixmap@red, pixmap@green, pixmap@blue, fn, PACKAGE="rtiff")
      } else if(class(pixmap) == "matrix") {
        pixmap = newPixmapRGB(pixmap, pixmap, pixmap);
        .Call("writeTiff", pixmap@red, pixmap@green, pixmap@blue, fn, PACKAGE="rtiff")
      } else {
        stop(paste("writeTiff expects a pixmapRGB or matrix, got ", class(pixmap)))
      }
      gc();
      return();
}
