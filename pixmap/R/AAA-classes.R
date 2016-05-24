setClass("pixmap",
         representation(size="integer",
                        cellres="numeric",
                        bbox="numeric",
                        bbcent="logical"),
         prototype(size=integer(2),
                   cellres=numeric(2),
                   bbox=numeric(4)))

setClass("pixmapChannels",
         representation(channels="character"),
         contains="pixmap")

setClass("pixmapGrey",
         representation(grey="matrix"),
         contains="pixmapChannels",
         prototype=prototype(new("pixmap"), channels="grey"))


setClass("pixmapIndexed",
         representation(index="matrix",
                        col="character"),
         contains="pixmap",
         prototype=prototype(new("pixmap")))

setClass("pixmapRGB",
         representation(red="matrix",
                        green="matrix",
                        blue="matrix"),
         contains="pixmapChannels",
         prototype=prototype(new("pixmap"),
                             channels=c("red", "green", "blue")))




