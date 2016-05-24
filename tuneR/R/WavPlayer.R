setWavPlayer <- function(player){
    options(WavPlayer = player)
}

getWavPlayer <- function(){
    options("WavPlayer")$WavPlayer
}
