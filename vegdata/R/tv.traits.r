tv.eco <- function (...) stop('This function is deprecated, use tv.traits instead.')
meanTraits <- function (...) stop('This function is deprecated, use ics with method "mean" instead.')

tv.traits <- function (db, trait.db = 'ecodbase.dbf', refl, ...) {
    tv_home <- tv.home()
    if(missing(refl))  refl <- if(missing(db)) tv.refl() else tv.refl(db = db)
    ecodb <- read.dbf(file.path(tv_home, 'Species', refl, trait.db), as.is = TRUE)
    empty <- function(x) all(is.na(x) | x == 0)
    na <- apply(ecodb, 2, empty)
    if(any(na)) {
#       if(!quiet) {
#         cat("\n The following columns contain no data and are omitted: \n")
#         cat(names(ecodb)[na])}
        ecodb <- ecodb[, !na]
                }
#    if(!quiet) message("Changing character fields into logical, integer or numericals if appropriate.")
# ecoDB <- apply(ecodb, 2, function(x) type.convert(as.character(x))) # doesnt work 
    ecoDB <- ecodb
#    for(i in 1:ncol(ecodb)) if(is.factor(ecodb[,i])) {
#       ecoDB[,i] <- iconv(as.character(ecodb[,i], "CP437", ""))
#       ecoDB[,i] <- type.convert(ecoDB[,i]) }
#     for(i in 1:ncol(ecoDB))  if(class(ecodb[,i]) != class(ecoDB[,i]))
#       ecoDB$ABBREVIAT <- as.character(ecoDB$ABBREVIAT)
#       ecoDB$LETTERCODE <- as.character(ecoDB$LETTERCODE)
#       if(!quiet) message('Data type of ', names(ecoDB)[i], ' changed to ', class(ecoDB[,i]))    
    return(ecoDB)
}


# meanTraits <- function(trait, veg, refl, trait.db = 'ecodbase.dbf', join = 'LETTERCODE', zero.is.NA = TRUE, ...) {
#   cat('Maximum performance value:', max(veg, na.rm=TRUE),'\n')
#   if(missing(refl)) refl <- attr(veg, 'taxreflist')
#   if(is.null(refl)) refl <- tv.refl()
#   if(missing(trait.db)) {
#     trait.db <- 'ecodbase.dbf'
#     cat('Using trait database:', trait.db, '\n') 
#     }
#   eco <- tv.traits(trait.db = trait.db, refl=refl, quiet=TRUE)
#   if(!trait %in% names(eco)) stop(paste('This trait name does not occur in column names of', trait.db))
#   IV <- as.numeric(eco[,trait][match(names(veg), eco[,join])])
#   names(IV) <- names(veg)
#   if(zero.is.NA) IV[IV == 0] <- NA
#   ind <- !is.na(IV)
#   veg <- veg[,ind]; IV <- IV[ind]
#   out <- rowSums(t((t(veg) * IV)) / rowSums(veg))
#   return(out)
# }
# 
# 
# 
# 
