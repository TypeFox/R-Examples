firstDigit <- function(x) {
    trunc(x / 10^(floor(log(abs(x),10))))
}
# lengths (mi) of 141 major North American rivers
table(firstDigit(rivers))              # lengths in miles
table(firstDigit(1.61 * rivers))       # lengths in km
table(firstDigit(5280 * rivers))       # lengths in feet
