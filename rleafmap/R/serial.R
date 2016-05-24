# 
# serial <- function(FUN, ...){
#   
#   if(!is.character(FUN)){
#     FUN <- deparse(substitute(FUN))
#   }
#   
#   x <- list(...)
#   x.len.max <- max(sapply(x, length))
#   seq.len.max <- seq(1, x.len.max)
#   xb <- lapply(x, rep, length.out = x.len.max)
#   
#   ds.arg <- deparse(substitute(list(...)))
#   ds.arg <- substr(ds.arg, 6, nchar(ds.arg)-1)
#   ds.arg <- strsplit(ds.arg, ",")
#   n.arg <- length(ds.arg)
#   arg.comp <- matrix(NA, nrow = x.len.max, ncol = n.arg)
#   for(i in n.arg){
#     arg.comp[, i] <- paste(ds.arg[i], "[[", seq.len.max, "]]", sep = "")
#   }
#   arg.comp <- apply(arg.comp, 1, paste, collapse = "", sep = "")
#   
#   res <- paste(FUN, "(", arg.comp, ")", sep = "")
#   return(res)
# }
#   
# 
# 
# serial2 <- function(FUN, ...){
#   fo <- deparse(substitute(list(...)))
#   return(fo)
# }
# FUN <- mean
# L1 <- lapply(1:10, rep, 10)
# L2 <- lapply(5:1, rep, 3)
# x <- list(L1, L2)
# serial("mean", a = L1, b=3)
# serial2(FUN="mean", a=L1, b=3)
# 
# 
# x <- function(y, z) print(deparse(substitute(y)))
# x(y= sum)
# 
# 
# writeLines(c("# hello markdown", "```{r hello-random, echo=TRUE}", "rnorm(5)", "```"), 
#            "test.Rmd")
# knit2html(text=c("# hello markdown", "```{r hello-random, echo=TRUE}", "rnorm(5)", "```"),
#           options = 'fragment_only')
# library(knitr)
# ss <- knit2html(output=NULL, text="# hello markdown\n```{r hello-random, echo=FALSE}\nplot(3)\n```",
#                 options = c('fragment_only', 'base64_images'))
# write(x=ss, "tt.html")
