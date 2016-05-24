inspectCorpus <- function() {
    setBusyCursor()
    on.exit(setIdleCursor())

    objects <- .getCorpusWindow()
    window <- objects$window
    txt <- objects$txt
    listbox <- objects$listbox

    tkwm.title(window, .gettext("Current Corpus"))

    mark <- 0

    tktag.configure(txt, "heading", font="sans 13 bold")
    tktag.configure(txt, "articlehead", font="sans 12 bold")
    tktag.configure(txt, "details", font="sans 10 italic")
    tktag.configure(txt, "small", font="sans 5")
    tktag.configure(txt, "fixed", font="courier 11")

    tkinsert(txt, "end", paste(sprintf(.gettext("Current corpus contains %i documents and %i terms."),
                                       nrow(dtm), ncol(dtm)), "\n\n", sep=""), "body")

    for(i in seq_along(corpus)) {
        id <- names(corpus)[i]
        tkinsert(txt, "end", paste(id, "\n", sep=""),
                 "articlehead")
        tkmark.set(txt, paste("mark", mark, sep=""), tkindex(txt, "insert-1c"))
        mark <- mark + 1
        tkinsert(listbox, "end", id)

        origin <- meta(corpus[[id]], "Origin")
        date <- meta(corpus[[id]], "DateTimeStamp")
        if(length(origin) > 0 && length(date) > 0)
            tkinsert(txt, "end", paste(origin, " - ", date, "\n", sep=""), "details")
        else if(length(origin) > 0)
            tkinsert(txt, "end", paste(origin, "\n", sep=""), "details")
        else if(length(origin) > 0)
            tkinsert(txt, "end", paste(date, "\n", sep=""), "details")

         if(length(origin) > 0 || length(date) > 0)
            tkinsert(txt, "end", "\n", "small")

        tkinsert(txt, "end", paste(paste(corpus[[id]], collapse="\n"), "\n\n"), "body")
    }

    # Only raise the window when we're done, as filling it may take some time
    tkraise(window)
}

