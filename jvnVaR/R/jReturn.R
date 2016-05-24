jReturn <-
function(s){
# Asset return
leng <- length(s)
object <- (s[1:leng-1]-s[2:leng]) / s[2:leng]
return(object)
}
