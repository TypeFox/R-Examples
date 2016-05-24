ft2miles <- 
function(x) x/5280

miles2ft <-
function(x) x*5280

ft2km <-
function(x) x/3280.839895

km2ft <-
function(x) x*3280.839895

add.alpha <- 
function (hex.color.list,alpha) sprintf("%s%02X",hex.color.list,floor(alpha*256))
