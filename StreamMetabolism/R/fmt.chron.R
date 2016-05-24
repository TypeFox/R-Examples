`fmt.chron` <-
function(x) {
   chron(sub(" .*", "", x), gsub(".* (.*)", "\\1:00", x))
}

