### Functions for handling coded data

# parse a coding formula.
# Normally expects the form coded ~ (orig - center) / divisor,
# But any linear expression in one variable is OK.
.parse.coding = function(form) {
#print(form)    
    if (!inherits(form, "formula")) 
        stop("Coding formulas must be of class \"formula\"")
    if (length(form) < 3) 
        stop("Formula lacks a left-hand-side")
    nm = all.vars(form)
    if (length(nm) < 2)
        stop(paste("Error in coding formula:", .form2str(form), 
                   "\nCoded and uncoded names must differ"))
    names(nm) = c("coded", "orig")
    rhs = as.character(form)[3]
    a = eval(parse(text = sub(nm[2], "0", rhs)))
    b = eval(parse(text = sub(nm[2], "1", rhs)))
    d = 1 / (b - a)
    list(names = nm, const=c(center = signif(-a * d, 4), divisor = signif(d, 4)))
}

### figure out the "rsdes" attribute for given data
# arguments data and primary are REQUIRED
# block, if non-missing, is PROPOSED name(s) for blocks. If pmatched (case insensitively)
# to variables in data, those columns are designated blocking factors
.rsdesattr = function(data, primary, block, call) {
    rsd = list(primary=primary)
    if (!missing(block)) {
        bidx = pmatch(tolower(block), tolower(names(data)))
        bidx = bidx[!is.na(bidx)]
        if (length(bidx) > 0) {
            rsd$block = names(data)[bidx]
            blk = data[[bidx[1]]]
            if (length(bidx) > 1) 
                for (i in bidx[-1]) blk = paste(blk, data[[i]], sep=".")
            bidxs = split(1:nrow(data), blk)
        }
    }
    
    if (!missing(call)) rsd$call = call
    rsd
}

### Create a string that looks the same as a printed formula
.form2str = function(formula) {
    if (inherits(formula, "formula")) {
        formula = as.character(formula)
        if (length(formula) == 3)
            formula = paste(formula[c(2,1,3)], collapse=" ")
        else
            formula = paste(formula, collapse=" ")
    }
    formula
}


# Code the data in a data.frame; may specify as arguments or in a list
coded.data = function(data, ..., formulas=list(...), block="block") {
    CALL = match.call()
    nm = names(data)
    if (length(formulas) == 0) {
#        stop("must provide coding formulas")
# auto-generated codings ...       
        codables = nm[sapply(data, function(x) length(unique(x)) < 6)]
        if (any(!is.na(exc <- pmatch(tolower(block), tolower(codables)))))
            codables = codables[-exc]
        if (length(codables) == 0)
            stop("No codings supplied and no variables with 5 or fewer distinct values")
        for (i in 1:length(codables)) {
            rng = range(as.numeric(data[[codables[i]]]))
            ctr = round(mean(rng), 3)
            div = round(rng[2] - ctr, 3)
            formulas[[i]] = as.formula(paste(
                "x", i, "~(", codables[i], "-", ctr, ")/", div, sep=""))
            
        }
        warning("Automatic codings created -- may not be what you want")
    }

    codings = list()
    for (f in formulas) {
        attr(f, ".Environment") = .GlobalEnv # keeps it from showing env when printed
        info = .parse.coding(f)
        cod = info$names[["coded"]]
        org = info$names[["orig"]]
        codings[[cod]] = f
        if (!is.null(data[[org]])) {
            if (is.factor(data[[org]])) 
                data[[org]] = as.numeric(data[[org]])
            data[[org]] = (data[[org]] - info$const[["center"]]) / info$const[["divisor"]]
            nm[nm==org] = cod
        }
    }
    names(data) = nm
    attr(data, "design.info") = attr(data, "desnum") = attr(data, "run.order") = NULL  # will no longer obey "design" class
    attr(data, "codings") = codings
    attr(data, "rsdes") = .rsdesattr(data, primary=names(codings), block=block, CALL)
    if (!is.coded.data(data))
        class(data) = c("coded.data", "data.frame")
    data
}

# Add coding attributes to already-coded data
as.coded.data = function(data, ..., formulas=list(...), block="block") {
    if (!is.data.frame(data))
        stop("'data' must inherit from \"data.frame\"")
    CALL = match.call()
    if (length(formulas) == 0) {
        if (is.coded.data(data))
            formulas = codings(data)
        else {
            codable = sapply(data, function(x) zapsmall(c(mean(x), max(x)))[1] == 0)
            if (sum(codable) == 0)
                stop("No codings supplied and no variables look like coded variables")
            formulas = sapply(names(data)[codable], function(nm) 
                as.formula(paste(nm, "~", nm, ".as.is", sep="")))
            warning("Default codings created -- may not be what you want")
        }
    }
    codings = list()
    for (f in formulas) {
        attr(f, ".Environment") = .GlobalEnv
        info = .parse.coding(f)
        cod = info$names[["coded"]]
        codings[[cod]] = f
    }
    mismatch = is.na(match(names(codings), names(data)))
    if (any(mismatch))
        stop("mismatch between coded names and data names")
    attr(data, "design.info") = attr(data, "desnum") = attr(data, "run.order") = NULL  # will no longer obey "design" class
    attr(data, "codings") = codings
    attr(data, "rsdes") = .rsdesattr(data, primary=names(codings), block=block, CALL)
    if (!is.coded.data(data))
        class(data) = c("coded.data", "data.frame")
    data
}

is.coded.data = function(x)
    inherits(x, "coded.data")

print.coded.data = function(x, ..., decode = TRUE) {
    if (!decode) {
        print.data.frame (x, ...)
        cat ("\nVariable codings ...\n")
    }
    else {
        print.data.frame(decode.data(x), ...)
        cat ("\nData are stored in coded form using these coding formulas ...\n")
    }
    sapply (attr(x, "codings"), print, showEnv=FALSE)
    invisible (x)
}

# Transform coded data back to original scale
decode.data = function(data) {
    nm = names(data)
    codings = attr(data, "codings")
    if (!is.null(codings)) for (f in codings) {
        info = .parse.coding(f)
        cod = info$names[["coded"]]
        org = info$names[["orig"]]
        if (!is.null(data[[cod]])) {
            data[[cod]] = info$const[["divisor"]] * data[[cod]] + info$const[["center"]]
            nm[nm==cod] = org
        }
    }
    names(data) = nm
    attr(data, "codings") = NULL
    attr(data, "rsdes") = NULL
    cls = class(data)
    class(data) = cls[cls != "coded.data"]
    data
}

### Recode a set of coded data, e.g. to a new center
recode.data = function(data, ..., formulas = list(...)) {
    rsd = attr(data, "rsdes")
    ddc = decode.data(data)
    data = coded.data(ddc, formulas=formulas)
    rsd$primary = attr(data, "rsdes")$primary
    attr(data, "rsdes") = rsd
    data
}

# Convert values in X to original based on codings
# Returns an object of the same type (data.frame, matrix, or vector)
# names (or column names) of X must match those used in codings
code2val = function(X, codings) {
    if (is.matrix(X))
        Z = as.matrix (code2val(as.data.frame(X), codings))
    else if (is.vector(X)) {
        nm = names(X)
        X = as.data.frame(as.list(X))
        names(X) = nm
        X = code2val (X, codings)
        Z = as.numeric (as.matrix (X))
        names(Z) = names(X)
    }
    else if (is.data.frame(X)) {
        attr(X, "codings") = codings
        Z = decode.data(X)
    }
    else stop ("Can't convert this object")
    Z
}

# Convert values in X to original based on codings
# Returns an object of the same type (data.frame, matrix, or vector)
# names (or column names) of X must match those used in codings
val2code = function(X, codings) {
    if (is.matrix(X))
        Z = as.matrix (val2code (as.data.frame(X), codings))
    else if (is.vector(X)) {
        nm = names(X)
        X = as.data.frame(as.list(X))
        names(X) = nm
        X = val2code (X, codings)
        Z = as.numeric (as.matrix (X))
        names(Z) = names(X)
    }
    else if (is.data.frame(X)) {
        Z = coded.data(X, formulas=codings)
        attr(Z, "codings") = NULL
        cls = class(Z)
        class(Z) = cls[cls != "coded.data"]
    }
    else stop ("Can't convert this object")
    Z
}

# Primitive accessor methods
codings = function(object)
    UseMethod("codings")

# S3 method for codings in coded.data
codings.coded.data = function(object)
    attr(object, "coding")

"codings<-" = function(object, value) {
    as.coded.data(object, formulas = value)
}

# Needed because we lose some attributes when subsetting
# Also we remove codings of variables that are lost
"[.coded.data" = function(x, ...) {
    cod = attr(x, "codings")
    rsd = attr(x, "rsdes")
    cls = class(x)
    x = get("[.data.frame")(x, ...)
    if (!("coded.data" %in% cls)) 
        return (x)
    lost = which(is.na(match(nm <- names(cod), names(x))))
    for (i in lost) cod[[nm[i]]] <- NULL
    if (length(cod) > 0) {
        attr(x, "codings") = cod
        attr(x, "rsdes") = rsd
    }
    else { # no longer a coded dataset
        class(x) = cls[-1]
        attr(x, "codings") = attr(x, "rsdes") = NULL
    }
    x
}

# When renaming coded data, change the formulas accordingly
"names<-.coded.data" = function(x, value) {
    if (!is.coded.data(x)) stop("not a coded.data object")
    oldnm = attr(x, "names")
    cod = codings(x)
    for (i in 1:length(oldnm)) {
        if (value[i] != oldnm[i]) {
            idx = match(oldnm[i], names(cod))
            if (!is.na(idx)) {
                names(cod)[idx] = value[i]
                cod[[idx]][[2]] = as.name(value[i])
            }
        }
    }
    attr(x, "names") = value
    attr(x, "codings") = cod
    x
}


# This is the flip side of names<-.coded.data: rename the original variables
"truenames" = function(x) {
    UseMethod("truenames")
}

"truenames.coded.data" = function(x) {
    nm = names(x)
    if(is.coded.data(x)) {
        for (f in codings(x)) {
            vn = all.vars(f)
            if (!is.na(idx <- grep(vn[1], nm))) 
                nm[idx] = vn[2]
        }
    }
    nm
}

"truenames<-" = function(x, value) {
    UseMethod("truenames<-")
}
    
"truenames<-.coded.data" = function(x, value) {
    if (!is.coded.data(x)) stop("not a coded.data object")
    oldnm = newnm = attr(x, "names")
    oldtrue = truenames(x)
    cod = codings(x)
    whichcoded = match(names(cod), oldnm)
    for (i in 1:length(oldnm)) {
        if (i %in% whichcoded) { # replace in coding formulas
            fstr = paste(as.character(cod[[oldnm[i]]])[c(2,1,3)], collapse=" ")
            cod[[oldnm[i]]] = as.formula(gsub(oldtrue[i], value[i], fstr))
        }
        else
            newnm[i] = value[i]
    }
        
    attr(x, "names") = newnm
    attr(x, "codings") = cod
    x
}
