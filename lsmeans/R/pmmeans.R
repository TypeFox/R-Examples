### Support for overriding "ls" with "pm" in names

### general-purpose wrapper for creating pmxxxxx functions
.pmwrap = function(lsfcn, ...) {
    result = lsfcn(...)

        if (inherits(result, "ref.grid"))
        result = .sub.ls.pm(result)
    else if(inherits(result, "lsm.list")) {
        for (i in seq_along(result))
            result[[i]] = .sub.ls.pm(result[[i]])
        names(result) = gsub("^ls", "pm", names(result))
    }
    result
}

# returns an updated ref.grid or lsmobj with setName "ls..." replaced by "pm..."
.sub.ls.pm = function(object) {
    nm = object@misc$estName
    update(object, estName = gsub("^ls", "pm", nm))
}

### Exported implementations

pmmeans = function(...)
    .pmwrap(lsmeans, ...)

# uh, maybe not...   pmms = pmmeans

pmtrends = function(...)
    .pmwrap(lstrends, ...)


pmmip = function(...)
    lsmip(...)

pmm = function(...)
    lsm(...)

pmmobj = function(...)
    .pmwrap(lsmobj, ...)

pmm.options = function(...)
    lsm.options(...)

get.pmm.option = function(...)
    get.lsm.option(...)
