library(RCurl)
u = "http://localhost/~duncan"

curl = getCurlHandle()
r = dynCurlReader(curl)

x = getURLContent(u, header = r, curl = curl)
r$value()
rm(r, curl)
gc()


x = getURLContent(u, followlocation = TRUE)

curl = getCurlHandle()
x = getURLContent(u, followlocation = TRUE, curl = curl)


