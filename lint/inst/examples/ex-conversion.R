# conversion examples
library(stringr)
text <- "
hw <- function(){
my.msg <- 
'Hello
World'
cat(my.msg, '\n')
}
"
lines <- readLines(textConnection(text))
p <- parse(text=text, keep.source=TRUE)
(pd <- getParseData(p))
(l <- str_locate('Hello', string=lines))

(f <- locate2find(l))
(r <- find2replace(f))

(s <- lint:::find_string(parse.data=pd))
parse2find(s)
