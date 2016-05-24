
blockColors.orange <- c("yellow", "gold", "orange", "DarkOrange", 
                        "DarkOrange", "OrangeRed", "red", "DarkRed")

blockColors.violet <- c("goldenrod", "peru", "salmon", "IndianRed", 
                        "maroon", "VioletRed", "DarkMagenta", "DarkViolet")

blockColors.blue <- c("green", "LawnGreen", "GreenYellow", "PaleGreen", 
                      "aquamarine", "turquoise", 
                      "DarkTurquoise", "DeepSkyBlue", 
                      "DodgerBlue","RoyalBlue", "SlateBlue", "MediumOrchid", 
                      "orchid", "violet", "pink", "MistyRose")

blockColorsMatrix <- cbind(blockColors.orange[c(1:8, 1:8)], 
                           blockColors.violet[c(1:8, 1:8)], 
                           blockColors.blue[c(1:16)])

blockColors <- c(t(blockColorsMatrix[,c(1,3)]))

