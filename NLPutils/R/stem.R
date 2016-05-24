Snowball_Stem_Annotator <-
function(language = "porter")
{
    f <- function(s) wordStem(s, language)
    description <-
        sprintf("A Snowball word stem annotator for language '%s'",
                language)
    Simple_Stem_Annotator(f, list(description = description))
}
