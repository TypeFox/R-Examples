# copied from Rdpack on 2014-03-30 and modified
#
#    TODO: remove the original from Rdpack after publishing "gbutils" on CRAN
#
                                                                       # 2014-04-06 cleaned up
parse_text <- function(text, ...,  keep = TRUE){            # see comments in parse_usage_text
   ks <- getOption("keep.source")
   if(!identical(keep, ks)){
       on.exit(options(keep.source = ks)) # restore previous value on exit
       options(keep.source = keep)
   }

   parse(text=text, ...)
}
