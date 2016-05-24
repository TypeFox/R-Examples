## ----echo=FALSE----------------------------------------------------------
library('knitr')

## ------------------------------------------------------------------------
set.seed(4711) # better reproducibility

library('pubprint') # load library
pp_opts_out$set(pp_init_out("plain")) # better readability in this document

a <- rnorm(40)
b <- a + .3

pprint(t.test(a, b))
pprint(cor.test(a, b))

## ------------------------------------------------------------------------
args(pubprint:::pprint)

## ------------------------------------------------------------------------
pprint(t.test(a, b))
pprint(t.test(a, b),
       concat = FALSE)
pprint(t.test(a, b),
       separator = NULL)
pprint(t.test(a, b),
       concat = FALSE,
       separator = NULL)
pprint(t.test(a, b),
       mmode = FALSE)
pprint(pprint(t.test(a, b),
              concat = FALSE,
              separator = NULL)[c(1, 2)])

## ------------------------------------------------------------------------
pprint(t.test(a, b),
       format = "object")
pprint(t.test(a, b),
       format = "t.test")
pprint(t.test(a, b),
       format = "chisquared")

## ------------------------------------------------------------------------
pprint(list(t.test(a, b), 0.2828363))

## ------------------------------------------------------------------------
pprint(t.test(a, b), print.estimate = FALSE)
pprint(t.test(a, b), estimate.names = c("control", "treatment"))

## ------------------------------------------------------------------------
ppo <- pubprint()

## ------------------------------------------------------------------------
push(ppo) <- t.test(a, b)
pull(ppo)

## ------------------------------------------------------------------------
args(pubprint:::`push<-.pubprint`)
args(pubprint:::pull.pubprint)

## ----error=TRUE----------------------------------------------------------
# save items in pipe and named memory
push(ppo) <- t.test(a)
push(ppo) <- t.test(a, b)
push(ppo, item = "i1") <- t.test(a, b + .2)

# retrieve items from pipe
pull(ppo) # item is removed from pipe
pull(ppo) # here as well
pull(ppo) # error because there are no more items in pipe

# retrieve items from named memory
pull(ppo, item = "i1") # item is not removed
pull(ppo, item = "i1", remove = TRUE) # item is removed
pull(ppo, item = "i1") # error, item does not more exist

## ------------------------------------------------------------------------
push(ppo) <- t.test(a)
push(ppo) <- t.test(a, b)
# add to last pipe item (n = 1 is default)
push(ppo, add = TRUE) <- 0.2828363

pull(ppo) # retrieve one way t-test
pull(ppo) # retrieve two way t-test with Cohen's d

## ------------------------------------------------------------------------
pp_opts_out$set(pp_init_out("html"))
pprint(t.test(a, b))
pp_opts_out$set(pp_init_out("latex"))
pprint(t.test(a, b))

## ------------------------------------------------------------------------
pp_opts$set(mmode = FALSE)

## ------------------------------------------------------------------------
myttest <- function(...) return("Hello World!")
pp_opts_style$set("t.test" = myttest)
pprint(t.test(a, b))

