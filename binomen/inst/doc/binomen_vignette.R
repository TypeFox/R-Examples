## ----echo=FALSE----------------------------------------------------------
library("knitr")
hook_output <- knitr::knit_hooks$get("output")
knitr::knit_hooks$set(output = function(x, options) {
   lines <- options$output.lines
   if (is.null(lines)) {
     return(hook_output(x, options))  # pass to default hook
   }
   x <- unlist(strsplit(x, "\n"))
   more <- "..."
   if (length(lines)==1) {        # first n lines
     if (length(x) > lines) {
       # truncate the output, but add ....
       x <- c(head(x, lines), more)
     }
   } else {
     x <- c(if (abs(lines[1])>1) more else NULL,
            x[lines],
            if (length(x)>lines[abs(length(lines))]) more else NULL
           )
   }
   # paste these lines together
   x <- paste(c(x, ""), collapse = "\n")
   hook_output(x, options)
 })

knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)

## ----eval=FALSE----------------------------------------------------------
#  install.packages("binomen")

## ----eval=FALSE----------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("ropensci/binomen")

## ------------------------------------------------------------------------
library("binomen")

## ------------------------------------------------------------------------
(obj <- make_taxon(genus="Poa", epithet="annua", authority="L.",
  family='Poaceae', clazz='Poales', kingdom='Plantae', variety='annua'))

## ------------------------------------------------------------------------
obj$binomial

## ------------------------------------------------------------------------
obj$binomial$authority

## ------------------------------------------------------------------------
obj$grouping

## ------------------------------------------------------------------------
obj$grouping$family

## ------------------------------------------------------------------------
obj %>% pick(family)
obj %>% pick(family, genus)

## ------------------------------------------------------------------------
obj %>% pop(family)
obj %>% pop(family, genus)

## ------------------------------------------------------------------------
obj %>% span(kingdom, family)

## ------------------------------------------------------------------------
gethier(obj)

## ------------------------------------------------------------------------
df <- data.frame(order = c('Asterales','Asterales','Fagales','Poales','Poales','Poales'),
  family = c('Asteraceae','Asteraceae','Fagaceae','Poaceae','Poaceae','Poaceae'),
  genus = c('Helianthus','Helianthus','Quercus','Poa','Festuca','Holodiscus'),
  stringsAsFactors = FALSE)
(df2 <- taxon_df(df))

## ------------------------------------------------------------------------
df2 %>% pick(order)

## ------------------------------------------------------------------------
df2 %>% pick(order, family, genus)

## ------------------------------------------------------------------------
df2 %>% span(family, genus)

## ----output.lines=1:20---------------------------------------------------
scatter(df2)

## ------------------------------------------------------------------------
out <- scatter(df2)
assemble(out)

