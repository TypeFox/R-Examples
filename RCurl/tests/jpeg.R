library(RCurl)

url = "http://www.omegahat.net/Rcartogram/demo.jpg"
header = RCurl:::basicTextGatherer()
ans = getBinaryURL(url, headerfunction = header$update)


url = "http://eeyore.ucdavis.edu/stat141/data/Elections/NYExitData.rda"
header = RCurl:::basicTextGatherer()
ans = getBinaryURL(url, headerfunction = header$update)
