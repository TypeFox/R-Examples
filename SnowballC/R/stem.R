wordStem <- function(words, language = "porter")
{
  words <- as.character(words)
  language <- as.character(language[1])

  .Call("R_stemWords", words, language, PACKAGE="SnowballC")
}


getStemLanguages <- function()
{
 .Call("R_getStemLanguages", PACKAGE="SnowballC")
}
