crr.Lexis <-
function( obj, mod, quiet=FALSE, ... )
{
# Model formula to be transmitted
md <- mod
# Outcome variable
fc <- as.character( mod[2] )
# Censoring levels
cn <- unique( as.character(obj$lex.Xst)[obj$lex.Xst==obj$lex.Cst] )
# The competing causes (the remaining levels)
cc <- setdiff( levels(obj$lex.Xst), union(fc,cn) )
# No l.h.s. side of formula when deriving model matrix
mod[2] <- NULL
# Remember no intercept term
cv <- model.matrix( mod, data=obj )[,-1]
# Then do it
M <- crr( ftime = obj$lex.dur,
         fstatus = obj$lex.Xst,
         failcode = fc,
         cencode = cn,
         cov1 = cv,
         ... )
# A table of the no of transitions
N <- with( obj, table( Relevel( lex.Xst, c(fc,cc,cn) ) ) )
names( N ) <- paste( rep( c("Event:"," comp:"," cens:"),
                          c(length(fc),length(cc),length(cn)) ),
                      names(N) )
# add model an table to the resulting object
M <- c( M, list( model.Lexis = md,
                 transitions = cbind(N) ) )
# remember the class attribute (lost by doing "c")
class( M )  <- "crr"
# print overview if desired
if( !quiet )
cat(   "crr analysis of event", paste('"',fc,'"',sep=''),
     "\n               versus", paste('"',cc,'"',sep=''),
     "\n                 with", paste('"',cn,'"',sep=''),
             "as censoring.\n" )
     M
     }
