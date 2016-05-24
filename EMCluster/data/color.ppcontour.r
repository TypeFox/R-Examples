### For ppcontour only.
col.ppcontour <- list(
  point = list(re = rgb(115,   0,   0, maxColorValue = 255),
               gr = rgb(  0,  90,   0, maxColorValue = 255),
               bl = rgb(  0,   0, 195, maxColorValue = 255),
               or = rgb(225,  39,   0, maxColorValue = 255),
               cy = rgb(  0, 145, 145, maxColorValue = 255),
               vi = rgb(118,  23, 206, maxColorValue = 255),
               aq = rgb( 49, 119,  96, maxColorValue = 255)
              ),
  # symbol = c(8, 21, 17, 15, 3, 4, 5),
  # symbol = c(0, 1, 2, 15, 16, 17, 18),
  symbol = c("1", "2", "3", "4", "5", "6", "7"),
  line = list(trre = rgb(175,   0,   0, 50, maxColorValue = 255),
              trgr = rgb(  0, 150,   0, 50, maxColorValue = 255),
              trbl = rgb(  0,   0, 255, 50, maxColorValue = 255),
              tror = rgb(255,  69,   0, 50, maxColorValue = 255),
              trcy = rgb(  0, 175, 175, 50, maxColorValue = 255),
              trvi = rgb(138,  43, 226, 50, maxColorValue = 255),
              traq = rgb( 69, 139, 116, 50, maxColorValue = 255)
              ),
  fill = list(trRed = function(k = 64){
                        s <- seq(0, 200, length.out = k)
                        rgb(175, 0, 0, s, maxColorValue = 255)
                      },
              trGreen = function(k = 64){
                          s <- seq(0, 200, length.out = k)
                          rgb(0, 150, 0, s, maxColorValue = 255)
                        },
              trBlue = function(k = 64){
                          s <- seq(0, 200, length.out = k)
                          rgb(0, 0, 255, s, maxColorValue = 255)
                       },
              trOrange = function(k = 64){
                           s <- seq(0, 200, length.out = k)
                           rgb(255, 69, 0, s, maxColorValue = 255)
                         },
              trCyan = function(k = 64){
                         s <- seq(0, 200, length.out = k)
                         rgb(0, 175, 175, s, maxColorValue = 255)
                        },
              trViolet = function(k = 64){
                           s <- seq(0, 200, length.out = k)
                           rgb(138, 43, 226, s, maxColorValue = 255)
                        },
              trAqua = function(k = 64){
                         s <- seq(0, 200, length.out = k)
                         rgb(69, 139, 116, s, maxColorValue = 255)
                       }
             )
) # End of col.ppcontour.
