tts_google <- function(content, destfile=tempfile("Rtts", fileext = ".mp3")){

  .Deprecated("tts_ITRI",
              package = "Rtts",
              msg = "function 'tts_google' is deprecated due to API issue. Please use 'tts_ITIR' instead.\n(Your request was done with 'tts_ITRI' this time)")
  tts_ITRI(content)

}
