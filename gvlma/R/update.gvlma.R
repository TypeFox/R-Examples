"update.gvlma" <-
function(object, formula, ...)
{
   gvlmaobj <- object
   lmobj <- gvlmaobj[-grep("GlobalTest", names(gvlmaobj))]
   class(lmobj) <- "lm"
   extras <- match.call(expand.dots = FALSE)$...
   if (!missing(formula)) extras$formula <- formula
   args <- list(as.name("update"))
   args$object <- lmobj
   if (length(extras) > 0)
     {
       gvlmaspecial <- !is.na(match(names(extras),
                                    c("alphalevel","timeseq", "warn")))
       args <- append(args, extras[!gvlmaspecial])              
     }
   argscall <- as.call(args)
   newlm <- eval(argscall, parent.frame())
######################################################################
   indxalpha <- match("alphalevel", names(extras), nomatch = 0)
   if (indxalpha>0) alpha <- extras[[indxalpha]]
   else alpha <- gvlmaobj$GlobalTest$LevelOfSignificance
   indxtimeseq <- match("timeseq", names(extras), nomatch = 0)
   if (indxtimeseq > 0) timevals <- extras[[indxtimeseq]]
   else timevals <- gvlmaobj$GlobalTest$timeseq
   indxwarn <- match("warn", names(extras), nomatch = 0)
   if (indxwarn > 0) warn <- extras[[indxwarn]]
   else warn <- TRUE
   if (length(timevals) != length(fitted(newlm)))
     {
       if (warn)
         warning("\n** Timeseq of incorrect length, using 1:nobs.**\n")
       timevals <- 1:length(fitted(newlm))
     }
######################################################################   
   gvlma(newlm, alphalevel = alpha, timeseq = timevals)
}

