# issue spotted by Alexandre Sieira
# and his example code.
library(RJSONIO)

fromJSON("[ 1, null, 4 ]", asText = TRUE, simplify = TRUE)
fromJSON('[ "a", null, "d" ]', asText = TRUE, simplify = TRUE)

fromJSON("[ 1, null, 4 ]", asText = TRUE, simplify = TRUE, nullValue = -999L)
fromJSON('[ "a", null, "d" ]', asText = TRUE, simplify = TRUE, nullValue = "999")
