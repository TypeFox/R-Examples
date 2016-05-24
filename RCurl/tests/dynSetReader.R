#
# This is an example of where we read the header and then set the options
# in the curl handle to read the body ofthe response.
#

library(RCurl)

reader = dynCurlReader()
status = curlPerform(url = "http://ws.audioscrobbler.com/2.0/?method=library%2Egetartists&user=oso&api_key=de933ba987ea45f82a73acd0d82071c3",
                     headerfunction = reader$update,
                     curl = reader$curl())

try(Encoding(reader$value()))    # UTF-8, sometimes

