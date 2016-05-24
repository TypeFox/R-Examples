# 2013-12-08 removed.
#     The name of the function (with this implementation) is misleading
#     since other "locate" functions look recurisvely.
#
# Rdo_locate_predefined_section <- function(rdo, sec){                       # 2013-03-29
#     Rdo_which_tag_eq(rdo, sec)
# }

# ## > tools:::RdTags
# toolsdotdotdotRdTags <-
# function (Rd)
# {
#     res <- sapply(Rd, attr, "Rd_tag")
#     if (!length(res))
#         res <- character()
#     res
# }
                                                                         # removed: 2013-12-01
# .aux_Rdo_locate_leaves2 <- function(x, i){                                  #new: 2013-11-28
#
#     if(is.logical(x)){
#         if(isTRUE(x))
#             i
#         else
#             FALSE
#     }else if(is.list(x) && length(x)==0){        # :TODO: 2013-11-28 VERY dangewroug change!
#                                        # arguments section: has list() with Rd_Tag = "\dots"
#         i     # consider the empty list a leave!
#
#     }else{     # here x is a list of integer vectors
#         lapply(x, function(x){ c(i,x)})
#     }
# }
#
# Rdo_locate_leaves2 <- function(object, f = function(x) TRUE){
#     fxx <- function(x){
#         if(is.character(x)){
#             return(f(x))         # evaluate `f' for this leaf
#         }else if(is.list(x)){
#             if(length(x)==0)   # list()
#                 return(x)
#             wrk <- lapply(x, fxx)
#             for(i in seq_along(wrk))
#                 wrk[[i]] <- .aux_Rdo_locate_leaves2(wrk[[i]], i)
#             indx <- sapply(wrk, function(y) !identical(y,FALSE) )
#             wrk <- wrk[indx]
#             for(i in seq_along(wrk))  # 2013-11-28 :TODO: krapka!!!
#                 if(is.list(wrk[[i]]) && length(wrk[[i]])==0){
#                     va <- attr(wrk[[i]], "Rd_tag")
#                     if(!is.null(va))
#                         wrk[[i]] <- attr(wrk[[i]], "Rd_tag")
#                 }
#
#             wrk <- do.call("c",wrk)
#             return(wrk)
#         }else{
#             return(FALSE)     # todo: replace this with a function and additional argument
#         }
#     }
#
#     fxx(object)
# }

# removed: 2013-12-01
# .locate_item_label <- function(rdo, pos){#find the position of the label of the specified item
#     # 2013-11-28 "\\dots" in "\\ arguments":  Rdo_locate_leaves(rdo[[pos]])[[1]]
#     # 2013-12-01 was: Rdo_locate_leaves2(rdo[[pos]])[[1]]
#     #            this function became redundant after starting to use .ascharRd()
#     1
# }

### 2013-04-09 macham fromgbRd.R - izglezhda nisto ne se izpolzva veche.
                                                                 ####### begin file fromgbRd.R
# strRd <- function(rdo, indent = "\t", verbose_blanks = FALSE, omit_blanks = FALSE){
#     n <- length(rdo)
#     if(n==0)
#         return("Empty list")
#
#     tags <- tools:::RdTags(rdo)
#
#     res <- character(0)
#     for(i in 1:n){
#         if(is.list(rdo[[i]])){
#             wrk <- Recall(rdo[[i]], indent = indent
#                           , verbose_blanks = verbose_blanks
#                           , omit_blanks = omit_blanks
#                           )
#             res <- c(res, tags[i], paste(indent, wrk, sep="") )
#                                                             # if(is.null(tags[[i]])) browser()
#         }else{
#             wrk <- if(is.null(tags[[i]])){      "_ No Rd_tag"
#                    }else if(length(rdo[[i]])==1 && rdo[[i]] == "\n"){    # note: new lines are
#                        if(omit_blanks)          NULL                     # not only in "TEXT"
#                        else if(verbose_blanks)  "NEWLINE:"
#                        else                     rdo[[i]]
#                    }else if( length(grep("^[[:space:]]+$", rdo[[i]])) == 1 ){
#                        if(omit_blanks)          NULL
#                        else if(verbose_blanks)  paste("BLANK:", rdo[[i]])
#                        else                     rdo[[i]]
#                    }else                        paste("", rdo[[i]])
#
#             if(!is.null(wrk))
#                 res <- c(res, paste(tags[i], ":", wrk, sep=""))
#         }
#     }
#     res
# }
#
# compare_attributes <- function(x,y, nams){                                  # todo: not tested
#     xattr <- attributes(x)
#     yattr <- attributes(y)
#     if(identical(xattr,yattr))   # is this comparison meaningful?
#         return(TRUE)
#
#     if(missing(nams))
#         nams <- unique(c(names(xattr), names(yattr)))
#
#     tbl <- matrix(NA, ncol = 3, nrow = length(nams))
#     rownames(tbl) <- nams
#     colnames(tbl) <- c("x","y", "compare")
#     tbl[,1] <- nams %in% names(xattr)
#     tbl[,2] <- nams %in% names(yattr)
#
#     tbl[,3] <- sapply(nams, function(x) if(tbl[x,1] && tbl[x,2])
#                                              identical(xattr[[x]], yattr[[x]])
#                                         else FALSE)
#     structure(FALSE, details = tbl)
# }
#
# .compareRdo_elem <- function(xrdo, yrdo){
#     list( list(xrdo, yrdo) )
# }
#
# .compareRdo_core <- function(xrdo, yrdo){
#     if(identical(xrdo,yrdo))
#         res <- return( structure(TRUE, Rd_tag = attr(xrdo,"Rd_tag")) )
#
#     if( is.character(xrdo) && is.character(yrdo)){
#         if( xrdo == yrdo ) {             # todo: ignoring attributes for now
#             res <- TRUE
#             attr(res,"Rd_tag") <- attr(xrdo,"Rd_tag")
#         }else
#             res <- .compareRdo_elem(xrdo, yrdo)  # todo: refine! if both are non-lists
#         return(res)                          #       the values may be the same,
#     }                                        #       and only the attributes different.
#
#     if( !all(c(is.list(xrdo),is.list(yrdo))) ){
#         res <- .compareRdo_elem(xrdo, yrdo)  # todo: refine! if both are non-lists
#         return(res)                          #       the values may be the same,
#     }                                        #       and only the attributes different.
#     # xrdo and yrdo are lists below
#
#     xattr <- attributes(xrdo)
#     yattr <- attributes(yrdo)
#     same_attr <- identical( xattr, yattr)
#
#     # if(!same_attr) browser()
#
#     xn <- length(xrdo)
#     yn <- length(yrdo)
#
#     # todo: zero xn or yn
#
#     res <- vector(max(xn,yn), mode="list")  # todo: better handling!
#     attr(res,"Rd_tag") <- attr(xrdo,"Rd_tag")
#
#     attributes(res) <- attributes(xrdo)
#     if( identical(class(res),"Rd") )
#         class(res) <- "cmpRd"     # "list"
#
#     attr(res, ".compare_attr") <- same_attr   # todo: needs better handling.
#
#     if(xn >= yn) xc <- yc <- 1:yn
#     else         xc <- yc <- 1:xn
#
#     if(xn == yn) xindmore <- yindmore <- integer(0)
#     else         xindmore <- yindmore <- (min(xn,yn)+1) : max(xn,yn)
#
#     for( i in seq_along(xc) ){
#         res[[i]] <- Recall(xrdo[[ xc[i] ]], yrdo[[ yc[i] ]])
#     }
#
#     for( i in seq_along(xindmore) ){
#         if(xindmore[i] > xn)
#             res[[i]] <- .compareRdo_elem( NA, yindmore[i])
#         else if(yindmore[i] > yn)
#             res[[i]] <- .compareRdo_elem( xindmore[i], NA)
#         #else{ # this should not occur
#         #    res[[i]] <- Recall(xrdo[[ xc[i] ]], yrdo[[ yc[i] ]])
#         #}
#     }
#
#     res
# }
#
# # todo:   obrabotka na blanks i new lines i pri sravnenieto (?)
# #         obrabotka na comments?
# #     !!! option da vrasta samo "leaves" koito ne sa ednakvi.
# #     !!! da pechata i indeksite (optional?)
# #     !!! obrabotka na attributes ( v momenta prosto gi sravnyava s identical)
# #
# strcmpRd <- function(rdo, indent = "\t", verbose_blanks = FALSE, omit_blanks = FALSE){
#     n <- length(rdo)
#     if(n==0)
#         return("Empty list")
#
#     tags <- tools:::RdTags(rdo)              # note: tags is usually a character vector
#     if(is.list(tags))                        # but it can be a list if some elements are NULL.
#         tags <- sapply(tags,
#                        function(x){
#                            if(length(x)>1) stop("Tags must be of length 1.")
#                            if(!is.null(x) && length(x)==0)
#                                            stop("Tags cannot be of length 0 (unless NULL.")
#                            if(is.character(x)) return(x)
#                            if(is.null(x)  ||
#                               is.list(x) && length(x)==1 && is.null(x[[1]]))
#                                               return("NULLTAG")
#                            x[[1]]})
#
#     if(!is.character(tags)){
#         print("tags is not a character vector!")
#         browser()
#     }
#
#     # tags[is.null(tags)] <- "None"
#     # tags[tags==""] <- "Ouch"
#
#     if(n==1  &&  is.logical(rdo[[1]])){
#         res <- rdo[[1]]
#         return(res)
#     }
#
#     res <- character(0)
#     for(i in 1:n){
#         if(is.list(rdo[[i]])){
#             wrk <- Recall(rdo[[i]], indent = indent
#                           , verbose_blanks = verbose_blanks
#                           , omit_blanks = omit_blanks
#                           )
#             if(is.logical(wrk))
#                 res <- c(res, paste(if(is.null(tags[[i]])) "None1" else tags[[i]],
#                                     ":", wrk, sep="") )
#             else
#                 res <- c(res, tags[[i]], paste(indent, wrk, sep="") )
#
#         }else if(is.logical(rdo[[i]])){
#             wrktag <- if(is.null(tags[[i]])) "None2" else tags[[i]]
#             res <- c(res, paste( wrktag,      ":", rdo[[i]], sep="") )
#
#         }else if(is.character(rdo[[i]])){
#             res <- c(res, tags[[i]], paste(indent, rdo[[i]], sep="") )
#
#         }else{
#             res <- c(res, tags[[i]], paste(indent, rdo[[i]], sep="") )
#         }
#     }
#     res
# }
#
# print.cmpRd <- function(x,...){
#     if(isTRUE(x))
#         print(x)
#     else
#         cat(unlist(strcmpRd(x,...)), sep="\n")
#     invisible(x)
# }
                                                                 #######   end file fromgbRd.R

#                                                    # urdo - usage from Rdo file/object;
# compare_sig1 <- function(urdo, ucur){              # ucur -       generated from actual object
#
#     # browser()
#     # 2011-11-12 commenting out as both should be quoted
#     # urdo$defaults <- gsub("\"","", urdo$defaults) # krapka; todo: elements may have embedded
#     # ucur$defaults <- gsub("\"","", ucur$defaults) # quotes; need to investigate where to put
#     status <- identical(urdo, ucur)               # the blame happens
#     status
#
#     # Note: commented out the rest since currently only the truth value is used.
#     #
#     # obj_removed <- is.null(ucur) || is.na(ucur)
#     # obj_added   <- is.null(urdo) || is.na(urdo)
#     #
#     # if(obj_removed || obj_added)
#     #     return( structure( status, details = list( obj_removed = obj_removed
#     #                                              , obj_added   = obj_added
#     #                                              , rdo_usage               = urdo
#     #                                              , cur_usage               = ucur
#     #                                )) )
#     #
#     # identical_names <- urdo$name == ucur$name
#     #
#     # identical_argnames <- identical(urdo$argnames, ucur$argnames)
#     # identical_defaults <- identical(urdo$defaults, ucur$defaults)
#     # identical_formals <- identical_argnames & identical_defaults
#     #
#     # added_argnames <- ucur$argnames[ !(ucur$argnames %in% urdo$argnames) ]
#     # removed_argnames <- urdo$argnames[ !(urdo$argnames %in% ucur$argnames) ]
#     #
#     #                                     # note: !!! intersect() is binary operation
#     # s <- intersect( intersect(names(urdo$argnames), names(ucur$argnames)),
#     #                 intersect(names(urdo$defaults), names(ucur$defaults)) )
#     #
#     # unchanged_defaults <- urdo$defaults[ ucur$defaults[s] == urdo$defaults[s] ]
#     #
#     # names_unchanged_defaults <- names(unchanged_defaults)[unchanged_defaults]
#     #
#     # # todo: more details for the case when !identical, e.g. equal up to reordering,
#     # #       added/removed defaults
#     #
#     # structure( status, details = list( identical_names          = identical_names
#     #                                  , obj_removed              = obj_removed
#     #                                  , obj_added                = obj_added
#     #                                  , identical_argnames       = identical_argnames
#     #                                  , identical_defaults       = identical_defaults
#     #                                  , identical_formals        = identical_formals
#     #                                  , added_argnames           = added_argnames
#     #                                  , removed_argnames         = removed_argnames
#     #                                  , names_unchanged_defaults = names_unchanged_defaults
#     #                                  , rdo_usage                = urdo
#     #                                  , cur_usage                = ucur
#     #                    ))
# }

# # note: todo: !!! is.pairlist(NULL) gives TRUE, so NULL indeed is a pairlist with zero elem.
# #             but inherits(NULL, "pairlist") gives FALSE
# #
# # this function is probably redundant and is called only by get_usage (and with one element,
# #  so the return in the first `if' below is invoked (todo: check and remove if so!)
# pairlist2f_usage <- function(x, nams, S3class = "", S4sig = "", infix = FALSE, fu  = TRUE,
#                              verbose = TRUE){
#                                 # 2012-09-27 !!! smenyam inherits(x, "pairlist") s is.pairlist
#     if(is.pairlist(x)){   # TRUE also when x is NULL
#         return(pairlist2f_usage1(x, name = if(missing(nams)) "fun_1" else nams,
#                                  S3class = S3class, S4sig=S4sig, infix = infix, fu = fu))
#     }
#     if(missing(nams))
#        nams <- names(x)
#     if(is.null(nams)  &&  length(x) > 0){
#         if(verbose)
#             cat("Argument 'x' is not named, supplying dummy function name(s).\n")
#         nams <- paste("fun_", seq_along(x), sep="")
#     }
#
#     if(is.null(names(x)))          # note: S3class below needs to have the same length as nams
#         names(x) <- nams
#     mapply("pairlist2f_usage1", x, nams, S3class, S4sig, infix, fu)
# }

