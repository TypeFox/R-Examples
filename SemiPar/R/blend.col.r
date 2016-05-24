########## R function: blend.col ##########

# For obtaining a vector of blended colours
# between two given colours.

# Last changed: 19 NOV 2003

blend.col <- function(stt.col,end.col,n=100)
{

   stt.rgb <- as.vector(col2rgb(stt.col))
   end.rgb <-  as.vector(col2rgb(end.col))

   blend.mat <-   cbind(seq(stt.rgb[1],end.rgb[1],length=n),
                        seq(stt.rgb[2],end.rgb[2],length=n),
                        seq(stt.rgb[3],end.rgb[3],length=n))

   hex.vec <- NULL
   for (i in 1:n)
      hex.vec <- c(hex.vec,rgb(blend.mat[i,1],blend.mat[i,2],
                   blend.mat[i,3],maxColorValue=255))

   return(hex.vec)
}

########## End of blend.col ##########
