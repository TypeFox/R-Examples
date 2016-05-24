set_power <-
function(x)
{
    if (length(x) < 1L) return(set(set()))
    set_union(set(set()),
              as.set(unlist(lapply(seq_along(x),
                                   function(i)
                                   apply(combn(.as.list(x), i), 2L, as.set)
                                   ),
                            recursive = FALSE)
                     )
              )
}

gset_power <-
function(x)
{
    if (length(x) < 1L) return(set(set()))
    .make_gset_from_list(c(list(gset()),
                           unlist(lapply(seq_along(x),
                                         function(i) gset_combn(x, i)
                                         ),
                                  recursive = FALSE)
                           )
                         )
}

cset_power <-
function(x)
{
    if (length(x) < 1L) return(set(set()))
    .make_gset_from_list(c(list(gset()),
                           unlist(lapply(seq_along(x),
                                         function(i) cset_combn(x, i)
                                         ),
                                  recursive = FALSE)
                           )
                         )
}

