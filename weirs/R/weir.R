"weir.broadcrest" <-
function(..., type=c("TWRI3A5", "BOS", NA)) {
   type <- match.arg(type);
   if(type == "TWRI3A5") {
      return(weir3a5.broadcrest(...));
   } else if(type == "BOS") {
      return(weirbos.broadcrest(...));
   } else {
      stop("unable to dispatch broad-crested weir function");
   }
}




"weir.sharpcrest" <-
function(..., type=c("TWRI3A5", NA)) {
   type <- match.arg(type);
   if(type == "TWRI3A5") {
      return(weir3a5.sharpcrest(...));
   } else {
      stop("unable to dispatch sharp-crested weir function");
   }
}




