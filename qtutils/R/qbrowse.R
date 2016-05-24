
## "object browser" type tools


qbrowser <- function(namespaces = FALSE)
    ## if (namespaces), then loaded namespaces are shown (as well as unexported objects in them?)
{
    pkgs <- .packages()
    nsps <- loadedNamespaces()
    if (namespaces)
    {
        entries <- union(pkgs, nsps)
        hasns <- entries %in% nsps
        ans <- vector(mode = "list", length = length(entries) + 1)
        names(ans) <- c("Global Workspace",
                        paste(ifelse(hasns, "Namespace", "Package"), entries, sep = ":"))
        ans[[1]] <- .GlobalEnv
        for (i in seq_along(ans)[-1])
        {
            ans[[i]] <-
                if (hasns[i-1]) getNamespace(entries[i-1])
                else as.environment(paste("package", entries[i-1], sep = ":"))
        }
    }
    else
    {
        entries <- pkgs
        ans <- vector(mode = "list", length = length(entries) + 1)
        names(ans) <- c("Global Workspace", paste("Package", entries, sep = ":"))
        ans[[1]] <- .GlobalEnv
        for (i in seq_along(ans)[-1])
        {
            ans[[i]] <-
                as.environment(paste("package", entries[i-1], sep = ":"))
        }
    }
    w <- qstr(ans)
    w$resize(600, 400)
    w
}


qrecover <- function()
{
    if(.isMethodsDispatchOn()) {
        ## turn off tracing
        tState <- tracingState(FALSE)
        on.exit(tracingState(tState))
    }
    ## find an interesting environment to start from
    calls <- sys.calls()
    from <- 0L
    n <- length(calls)
    if (identical(sys.function(n), qrecover))
        ## options(error=qrecover) produces a call to this function as an object
        n <- n - 1L
    ## look for a call inserted by trace() (and don't show frames below)
    ## this level.
    for(i in rev(seq_len(n))) {
        calli <- calls[[i]]
        fname <- calli[[1L]]
        ## deparse can use more than one line
        if(!is.na(match(deparse(fname)[1L],
                        c("methods::.doTrace", ".doTrace")))) {
            from <- i-1L
            break
        }
    }
    ## if no trace, look for the first frame from the bottom that is not
    ## stop or recover
    if (from == 0L)
        for (i in rev(seq_len(n)))
        {
            calli <- calls[[i]]
            fname <- calli[[1L]]
            if(!is.name(fname) || is.na(match(as.character(fname), c("qrecover", "stop", "Stop"))))
            {
                from <- i
                break
            }
        }
    if (from > 0L && interactive())
    {
        if (identical(getOption("show.error.messages"), FALSE)) # from try(silent=TRUE)?
            return (NULL)
        calls <- limitedLabels(calls[1L:from])
        env.collection <- new.env(parent = emptyenv())
        for (i in seq_along(calls))
        {
            env.collection[[ sprintf("%03d: %s", i, calls[i]) ]] <- sys.frame(i)
        }
        qstr(env.collection)$show()
    }
    else
        cat(gettext("No suitable frames for qrecover()\n"))
}


