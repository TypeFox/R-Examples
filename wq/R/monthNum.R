monthNum <- 
function(y) {
  match(months(y, abbreviate = TRUE), month.abb)
}