setGeneric("play",
function(object, player, ...) standardGeneric("play"))

setMethod("play", signature(object = "character", player = "ANY"),
function(object, player, ...){
    if(missing(player)){
        player <- getWavPlayer()
        if(.Platform$OS.type == "windows" && is.null(player)){
            player <- "c:/Program Files/Windows Media Player/wmplayer.exe"
            if(!file.exists(player)){
                player <- "mplay32"
                if(missing(...))
                    player <- paste(player, "/play /close")
            } else {
                player <- shQuote(player)
            }
        }
    }
    if(.Platform$OS.type == "windows"){
        shell(paste('"', paste(player, ..., shQuote(object)), '"', sep=""))
    } else {
        system(paste(player, ..., object))
    }
})

setMethod("play", signature(object = "WaveGeneral", player = "ANY"),
function(object, player, ...){
    filename <- "tuneRtemp.wav"
    wd <- getwd()
    setwd(tempdir())
    on.exit({unlink(filename); setwd(wd)})
    writeWave(object, filename)
    play(file.path(tempdir(), filename), player, ...)
})
