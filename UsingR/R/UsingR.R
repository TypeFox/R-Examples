UsingR = function(what = c("webpage","errata","changes","exercises")) {

  what = match.arg(what)

  urls = c(
    "webpage"= "http://www.math.csi.cuny.edu/UsingR/",
    "errata"= "http://www.math.csi.cuny.edu/UsingR/Errata/",
    "changes"= "http://www.math.csi.cuny.edu/UsingR/Changes",
    "exercises"= "http://www.math.csi.cuny.edu/UsingR/AnswersToSelectedProblems"
    )
  url = urls[what]

  browseURL(url)

}

## alias
usingr = UsingR
