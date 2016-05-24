"pickFrom" <-  function (vec, nsets = 1, return.indices = FALSE,
                         setlabels = NULL,
                         edit.setlabels = TRUE,
                         subset = TRUE,
                         warningText = "one or more selections empty",
                         title = "Subset picker",
                         items.label = "Pick from",
                         labels.prompt = "Your label for this set",
                         list.height = 20,
                         items.scrollbar = TRUE,
                         preserve.order = TRUE,
                         graphics = TRUE,
                         listFont = "Courier 12",
                         labelFont = "Helvetica 11",
                         windowPos = "+150+30")
{
    if (!interactive())
        stop("Attempt to use interactive selection function when R is not ",
             "running interactively")

    if (!is.vector(vec))
        stop("argument `vec' muct be a vector")
    vec.is.numeric <- if (is.numeric(vec)) TRUE else FALSE
    vec.as.char <- as.character(vec)
    vec.to.pickfrom <- vec.as.char[subset]
    ni <- length(vec.to.pickfrom)

    if (is.character(subset)) subset <- match(subset(names(vec)))
    if (is.logical(subset)) subset <- seq(along = vec)[subset]

    setlabels <- if (!is.null(setlabels))
        as.list(setlabels)
    else as.list(rep("", nsets))

    items.label <- paste(items.label, ":", sep = "")

    if (graphics & capabilities("tcltk")) {
        requireNamespace("tcltk", quietly = TRUE)
        ppp <- NULL ## only to avoid a NOTE at package check time

        string.to.vector <- function(string.of.indices) {
            as.numeric(strsplit(string.of.indices, split = " ")[[1]])
        }
        base <- tcltk::tktoplevel(takefocus = 1)
        tcltk::tkwm.title(base, title)
        tcltk::tkwm.geometry(base, windowPos)
        tcltk::tkwm.resizable(base, 0, 0)

        right.frm <- tcltk::tkframe(base)
        left.frm <- tcltk::tkframe(base)

        items.list <- as.character(tcltk::tclVar(paste("{", paste(vec.to.pickfrom,
                                                           collapse = "} {"),
                                                "}", sep = "")))
        items.frm <- tcltk::tkframe(left.frm)

        items.label <- tcltk::tklabel(items.frm,
                               text = items.label,
                               anchor = "w",
                               justify = "left")
        tcltk::tkgrid(items.label, row = 0, columnspan = 2, sticky = "w")
        items.height <- min(list.height, ni)
        items.width <- max(8, max(nchar(vec.to.pickfrom)))
        items <- tcltk::tklistbox(items.frm,
                           listvar = items.list,
                           bg = "grey50",
                           selectmode = "extended",
                           fg = "white",
                           font = listFont,
                           width = items.width,
                           height = items.height)
        tcltk::tkgrid(items, row = 1, column = 0)
        preserve.order <- tcltk::tclVar(as.numeric(preserve.order))
        buttons.frm <- tcltk::tkframe(left.frm)
        buttonA <- tcltk::tkradiobutton(buttons.frm,
                                 text = "Sort sets in\nthe above order\nupon \"Add\"",
                                 justify = "left",
                                 variable = preserve.order,
                                 value = "1",
                                 command = function(){NULL}
                                 )

        buttonB <- tcltk::tkradiobutton(buttons.frm,
                                 text = "Place\nnewly added\nitems last",
                                 justify = "left",
                                 variable = preserve.order,
                                 value = "0",
                                 command = function(){NULL}
                                 )

        if (items.scrollbar && (length(vec) > items.height)) {
            items.scrollbar <- tcltk::tkscrollbar(items.frm,
                                           orient = "vertical",
                                           repeatinterval = 1,
                                           command = function(...) {
                                               tcltk::tkyview(items, ...)
                                           })
            tcltk::tkconfigure(items, yscrollcommand = function(...) {
                tcltk::tkset(items.scrollbar, ...)
                xy <- string.to.vector(tcltk::tclvalue(tcltk::tkget(items.scrollbar)))
                tcltk::tkyview.moveto(items, xy[1])
            })
            tcltk::tkgrid(items.scrollbar, row = 1, column = 1, sticky = "ns")
        }
        tcltk::tkpack(buttonA, buttonB, pady = 1, padx = 5, side = "top",
               anchor = "nw")
        tcltk::tkpack(items.frm, buttons.frm,
               pady = 1, padx = 5, side = "top")

        tcltk::tkpack(left.frm, side = "top", expand = "true", anchor = "n")

        sets.frm <- tcltk::tkframe(right.frm)
        setframe <- list()
        label <- list()
        setlabeltext <- list()
        labelentry <- list()
        TCLlabel <- list()
        listbox <- list()
        add.but <- list()
        labelbox <- list()
        listvarname <- list()
        remove.but <- list()
        tkset <- list()
        set <- list()
        Rtkset <- list()
        subset.height <- min(list.height - 5, ni)
        for (i in 1:nsets) {
            tkset[[i]] <- tcltk::tclVar("")
            TCLlabel[[i]] <- tcltk::tclVar(setlabels[[i]])
            setframe[[i]] <- tcltk::tkframe(sets.frm,
                                     width = 250,
                                     relief = "groove",
                                     borderwidth = 2)
            label[[i]] <- tcltk::tklabel(setframe[[i]], text = setlabels[[i]])
            listvarname[[i]] <- as.character(tkset[[i]])
            listbox[[i]] <- tcltk::tklistbox(setframe[[i]],
                                      listvar = listvarname[[i]],
                                      bg = "white",
                                      height = subset.height,
                                      font = listFont,
                                      width = items.width,
                                      selectmode = "extended")
            labelbox[[i]] <- tcltk::tkframe(setframe[[i]], width = 250)
            setlabeltext[[i]] <- tcltk::tklabel(labelbox[[i]], text =
                                         paste(labels.prompt, ":", sep = ""))
        }
        add.cmd <- deparse(function() {
            set[[ppp]] <- match(Tcl.to.R(tcltk::tclvalue(tkset[[ppp]])),
                                vec.to.pickfrom)
            set[[ppp]] <- union(set[[ppp]], 1 +
                            string.to.vector(tcltk::tclvalue(tcltk::tkcurselection(items))))
            if (as.logical(tcltk::tclObj(preserve.order)))
                set[[ppp]] <- sort(set[[ppp]])
            tcltk::tclvalue(tkset[[ppp]]) <- R.to.Tcl(vec.to.pickfrom[set[[ppp]]])
            tcltk::tkconfigure(add.but[[ppp]], state = "disabled")
        })
        remove.cmd <- deparse(function() {
            Rtkset[[ppp]] <- Tcl.to.R(tcltk::tclvalue(tkset[[ppp]]))
            out <- 1 +
                string.to.vector(tcltk::tclvalue(tcltk::tkcurselection(listbox[[ppp]])))
            if (length(Rtkset[[ppp]]) == length(out))
                tcltk::tclvalue(tkset[[ppp]]) <- ""
            else tcltk::tclvalue(tkset[[ppp]]) <- R.to.Tcl(Rtkset[[ppp]][-out])
            tcltk::tkconfigure(remove.but[[ppp]], state = "disabled")
            tcltk::tkselection.clear(listbox[[ppp]], "0", "end")
        })
        for (i in 1:nsets) {
            add.but[[i]] <- tcltk::tkbutton(setframe[[i]], text = "Add",
                                     fg = "darkgreen",
                                     disabledforeground = "darkgrey",
                                     width = 10,
                                     state = "disabled",
                                     command = eval(parse(text = gsub("ppp",
                                                          as.character(i),
                                                          add.cmd))))
            remove.but[[i]] <- tcltk::tkbutton(setframe[[i]], text = "Remove",
                                        fg = "darkred",
                                        disabledforeground = "darkgrey",
                                        width = 10,
                                        state = "disabled",
                                        command = eval(parse(text = gsub("ppp",
                                                             as.character(i),
                                                             remove.cmd))))
            labelentry[[i]] <- tcltk::tkentry(labelbox[[i]],
                                       textvariable = as.character(TCLlabel[[i]]),
                                       font = labelFont,
                                       bg = "white")
            if (edit.setlabels) {
                tcltk::tkpack(setlabeltext[[i]], labelentry[[i]], side = "top",
                       anchor = "w")
            }
            tcltk::tkpack(label[[i]], add.but[[i]], remove.but[[i]], listbox[[i]],
                   labelbox[[i]], side = "top", padx = 5, pady = 5)
            tcltk::tkpack(setframe[[i]], side = "left", padx = 3, pady = 10)
        }
        fun1 <- deparse(function() {
            if (tcltk::tclvalue(tcltk::tkcurselection(listbox[[ppp]])) != "") {
                for (j in 1:nsets) {
                    tcltk::tkconfigure(add.but[[j]], state = "disabled")
                }
                tcltk::tkconfigure(remove.but[[ppp]], state = "normal")
            }
            for (j in (1:nsets)[-ppp]) {
                tcltk::tkconfigure(remove.but[[j]], state = "disabled")
            }
            tcltk::tkfocus(listbox[[ppp]])
        })
        for (i in 1:nsets) {
            tcltk::tkbind(listbox[[i]], "<<ListboxSelect>>",
                   eval(parse(text = gsub("ppp", as.character(i), fun1))))
        }
        tcltk::tkbind(items, "<<ListboxSelect>>", function() {
            items.selected <- vec.to.pickfrom[1 +
                         string.to.vector(tcltk::tclvalue(tcltk::tkcurselection(items)))]
            for (i in 1:nsets) {
                set[[i]] <- Tcl.to.R(tcltk::tclvalue(tkset[[i]]))
                if (setequal(items.selected, intersect(items.selected,
                                                       set[[i]]))) {
                    tcltk::tkconfigure(add.but[[i]], state = "disabled")
                }
                else tcltk::tkconfigure(add.but[[i]], state = "normal")
                tcltk::tkconfigure(remove.but[[i]], state = "disabled")
            }
        })
        tcltk::tkbind(items, "<Button-1>", function() tcltk::tkfocus(items))
        buttons.frame <- tcltk::tkframe(right.frm)
        OK <- tcltk::tclVar(0)
        ok.but <- tcltk::tkbutton(buttons.frame, text = "OK", width = 10,
                           command = function() {
                               tcltk::tkdestroy(base)
                               tcltk::tclvalue(OK) <- 1
                           })
        tcltk::tkconfigure(ok.but, state = "normal")
        cancel.but <- tcltk::tkbutton(buttons.frame, text = "Cancel", width = 10,
                               command = function() {
                                   tcltk::tkdestroy(base)
                               })
        tcltk::tkpack(ok.but, cancel.but, side = "left", padx = 20, pady = 20)
        tcltk::tkpack(sets.frm, buttons.frame, side = "top")
        tcltk::tkpack(left.frm, side = "left", anchor = "nw", padx = 1)
        tcltk::tkpack(right.frm, anchor = "ne")
        tcltk::tkwait.window(base)
        tcltk::.Tcl("update idletasks")
        if (tcltk::tclvalue(OK) == "1") {
            sets <- lapply(tkset, function(set) {
                match(Tcl.to.R(tcltk::tclvalue(set)), vec.to.pickfrom)
            })
            if (any(sapply(sets, length) == 0)) {
                warning(warningText)
            }
            labels <- lapply(TCLlabel, tcltk::tclvalue)
            names(sets) <- labels
            result <- sets
        } else return(NULL)
    }
    else {
        result <- list()
        cat("**", title, "**\n\n")
        cat(items.label, "\n")
        op <- paste(format(seq_len(ni)), ": ", vec.to.pickfrom, sep = "")
        if (ni > 10) {
            fop <- format(op)
            nw <- nchar(fop[1], "w") + 2
            ncol <- getOption("width")%/%nw
            if (ncol > 1)
                op <- paste(fop, c(rep("  ", ncol - 1), "\n"),
                            sep = "", collapse = "")
            cat("", op, sep = "\n")
        }
        else cat("", op, "", sep = "\n")
        cat("Enter sequences of numbers separated by commas:\n")
        for (i in 1:nsets) {
            ind <- readline(paste(ifelse(nchar(setlabels[[i]]),
                                         setlabels[[i]],
                                         paste("Set", i)), ": ", sep = ""))
            ind <- eval(parse(text = paste("c(", ind, ")")))
            if (edit.setlabels){
                tmp <- readline(paste(labels.prompt,
                                      ifelse(nchar(setlabels[i]),
                                             paste(" [", setlabels[i], "]",
                                                   sep = ""), ""),
                                      ": ", sep = ""))
                if (nchar(tmp)) setlabels[i] <- tmp
            }
            if (all(invalid <- !ind %in% seq(ni)))
                result[[i]] <- numeric(0)
            else if (any(invalid)){
                warning("Ignored invalid selection(s): ",
                        paste(ind[invalid], sep = ", "),
                        ".\n", immediate. = TRUE)
                result[[i]] <- ind[!invalid]
            }
            else
                result[[i]] <- ind
            if (!preserve.order) result[[i]] <- sort(result[[i]])
        }
        if (!all(sapply(result, length)))
            warning(warningText)
        names(result) <- setlabels
    }
    return(
           if (return.indices)
           lapply(result, function(set) subset[set])
           else lapply(result, function(set) (vec[subset])[set])
           )
}
