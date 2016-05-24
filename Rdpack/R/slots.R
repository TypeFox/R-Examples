inspect_slots <- function(rdo, final=TRUE){
    fullname <- .get.name_content(rdo)$name     # name of the class with '-class' suffix
    cur <- .capture_promptAny(fullname, final=final)
    curnames <- .get_top_labels(cur, "Slots")   # current slots

    rdonames <- .get_top_labels(rdo, "Slots")   # slots in rdo

    icmp <- .asym_compare(rdonames, curnames)   # compare; get fields $i_new,  $i_rem, $i_com

    if(length(icmp$i_new)>0){

        ## 2014-06-18 the labels contain Rd markup, as in \code{mo.col}, but this seems
        ##            harmless.

        ## 2014-06-21 new slots were not handled correctly, see the changes below;
        ##            see also inspect_clmethods - there is a lot of common ground.
        ##            todo: consolidate?

        newnames <- names(icmp$i_new)  # 2014-06-21 new

        cat("Undocumented slots:", newnames  # 2014-06-21 was: names(icmp$i_new)
                                             # 2014-06-18 was: icmp$i_new
            , "\n")
        cat("\tAdding items for them.\n")
                                                           # todo: insert in particular order?
        cnt_newslots <- .get_top_items(cur, "Slots", newnames) # 2014-06-21 was: icmp$i_new
        cnt_newslots <- .nl_and_indent(cnt_newslots) # 2014-06-21 new

                                        # this ensures that the closing brace for "describe"
                                        # is on new line. todo: this needs to be indented.
        cnt_newslots <- c(cnt_newslots, list(Rdo_newline())) # 2014-06-21 new

        dindx <- .locate_top_tag(rdo, "Slots")

        rdo <- append_to_Rd_list(rdo, cnt_newslots, dindx)
    }

    if(length(icmp$i_removed)>0){             # todo: maybe put this note in a section in rdo?
        cat("Slots:", icmp$i_removed, "\n")
        cat("are no longer present. Please remove their descriptions manually.\n")
    }
    rdo
}
