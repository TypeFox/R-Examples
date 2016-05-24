# with contributions by Adam Kramer and Michael Friendly (originally by J. Fox)
# last modified 2015-06-09 by J. Fox

globalVariables("dot")

Greek <- read.table(system.file("etc/GreekLetters.txt", package="sem"), as.is=TRUE)

math <- function(text, html.only=FALSE, hat=FALSE){
    if (length(text) > 1) {
        result <- sapply(text, math, html=html.only, hat=hat)
        names(result) <- names(text)
        return(result)
    }
#     subscripts <- c("&#8320;", "&#8321;", "&#8322;", "&#8323;", "&#8324;", "&#8325;", "&#8326;",
#         "&#8327;", "&#8328;", "&#8329;")
    subscripts <- c("&#x2080;", "&#x2081;", "&#x2082;", "&#x2083;", "&#x2084;", "&#x2085;", "&#x2086;",
        "&#x2087;", "&#x2088;", "&#x2089;")
    superscripts <- c("&#x2070;", "&sup1;", "&sup2;", "&sup3;", "&#x2074;", "&#x2075;", 
        "&#x2076;", "&#x2077;", "&#x2078;", "&#x2079;")
    names(subscripts) <- names(superscripts) <- 0:9
    hat <- if (hat) "&#770;" else ""
    text <- gsub(" ", "", text)
    symbol <- regexpr("^[a-zA-Z]+", text)
    if (symbol != 1) stop("text does not start with an alphabetic symbol name")
    symbol <- if (html.only) {
        paste0("&", substring(text, 1, attr(symbol, "match.length")), ";")
    }
    else{
        s <- substring(text, 1, attr(symbol, "match.length"))
        s <- Greek[s, "decimal"]
        if (is.na(s)) stop(s, " is not a Greek letter")
        s
    }
    subscript <- regexpr("_\\{", text)
    subscript <- if (subscript >= 1){
        subscript <- substring(text, subscript + 2)
        endbrace <- regexpr("\\}", subscript)
        if (endbrace < 1) stop("unmatched closing brace in ", text)
        substring(subscript, 1, endbrace - 1)
    }
    else ""
    if (subscript != ""){
        subscript <- unlist(strsplit(subscript, split=""))
        subscript <- subscripts[subscript]
        if (any(is.na(subscript))) stop ("invalid non-numeral subscript")
        subscript <- paste(subscript, collapse="")
    }
    superscript <- regexpr("\\^\\{", text)
    superscript <- if (superscript >= 1){
        superscript <- substring(text, superscript + 2)
        endbrace <- regexpr("\\}", superscript)
        if (endbrace < 1) stop("unmatched closing brace in ", text)
        substring(superscript, 1, endbrace - 1)
    }
    else ""
    if (superscript != ""){
        superscript <- unlist(strsplit(superscript, split=""))
        superscript <- superscripts[superscript]
        if (any(is.na(superscript))) stop ("invalid non-numeral superscript")
        superscript <- paste(superscript, collapse="")
    }
    paste0(symbol, hat, subscript, superscript)
}

path.diagram <- function(...) {
    .Deprecated("pathDiagram", package = "sem")
    pathDiagram(...)
}

pathDiagram <- function (model, ...)
{
    UseMethod("pathDiagram")
}

pathDiagram.semmod <- function(model, obs.variables, ...) {
    parse.path <-
        function(path) {
            path.1 <- gsub("-", "", gsub(" ","", path))
            direction <- if (regexpr("<>", path.1) > 0)
                2
            else if (regexpr("<", path.1) > 0)
                - 1
            else if (regexpr(">", path.1) > 0)
                1
            else
                stop(paste("ill-formed path:", path))
            path.1 <- strsplit(path.1, "[<>]")[[1]]
            list(first = path.1[1], second = path.1[length(path.1)], direction =
                    direction)
        }
    if ((!is.matrix(model)) |
            ncol(model) != 3)
        stop ("model argument must be a 3-column matrix")
    startvalues <- as.numeric(model[,3])
    par.names <- model[,2]
    n.paths <- length(par.names)
    heads <- from <- to <- rep(0, n.paths)
    for (p in 1:n.paths) {
        path <- parse.path(model[p,1])
        heads[p] <- abs(path$direction)
        to[p] <- path$second
        from[p] <- path$first
        if (path$direction == -1) {
            to[p] <- path$first
            from[p] <- path$second
        }
    }
    ram <- matrix(0, p, 5)
    all.vars <- unique(c(to, from))
    latent.vars <- setdiff(all.vars, obs.variables)
    vars <- c(obs.variables, latent.vars)
    pars <- na.omit(unique(par.names))
    ram[,1] <- heads
    ram[,2] <- apply(outer(vars, to, "=="), 2, which)
    ram[,3] <- apply(outer(vars, from, "=="), 2, which)
    par.nos <- apply(outer(pars, par.names, "=="), 2, which)
    if (length(par.nos) > 0)
        ram[,4] <-
        unlist(lapply(par.nos, function(x)
            if (length(x) == 0)
                0
            else
                x))
    ram[,5] <- startvalues
    colnames(ram) <-
        c("heads", "to", "from", "parameter", "start value")
    pars <- unique(na.omit(par.names))
    coeff <- rep(0, length(pars))
    names(coeff) <- pars
    fake.sem <-
        list(
            ram = ram, n = length(obs.variables), var.names = vars, coeff = coeff,
            semmod = model
        )
    class(fake.sem) <- "sem"
    pathDiagram(fake.sem, ...)
}


pathDiagram.sem <-
    function (model, file = "pathDiagram", style = c("ram", "traditional"),
        output.type = c("html", "graphics", "dot"), graphics.fmt = "pdf", dot.options = NULL,
        size = c(8, 8), node.font = c("Helvetica", 14),
        edge.font = c("Helvetica", 10), digits = 2, rank.direction = c("LR", "TB"),
        min.rank = NULL, max.rank = NULL, same.rank = NULL,
        variables = model$var.names, var.labels, parameters, par.labels,
        ignore.double = TRUE, ignore.self = FALSE, error.nodes = TRUE,
        edge.labels = c("names", "values", "both"),  edge.colors = c("black", "black"),
        edge.weight = c("fixed", "proportional"),
        node.colors = c("transparent", "transparent", "transparent"),
        standardize = FALSE, ...) {
            Dot <- function(..., semicolon = TRUE, newline = TRUE) {
                cat(file = handle, paste0(..., if (semicolon)
                    ";"
                    else
                        "",
                    if (newline)
                        "\n"
                    else
                        ""))
            }
            style <- match.arg(style)
            output.type <- match.arg(output.type)
            edge.labels <- match.arg(edge.labels)
            edge.weight <- match.arg(edge.weight)
            rank.direction <- match.arg(rank.direction)
            if (output.type == "html") {
                handle <- textConnection("dot", "w")
            }
            else {
                dot.file <- paste0(file, ".dot")
                handle <- file(dot.file, "w")
                if (output.type == "graphics")
                    graph.file <- paste0(file, ".", graphics.fmt)
            }
            on.exit(close(handle))
            Dot("digraph \"", deparse(substitute(model)), "\" {", semicolon = FALSE)
            Dot("  rankdir=", rank.direction)
            Dot("  size=\"", size[1], ",", size[2], "\"")
            Dot(
                "  node [fontname=\"", node.font[1],
                "\" fontsize=", node.font[2], " fillcolor=\"", node.colors[1],
                "\" shape=box style=filled]"
            )
            Dot("  edge [fontname=\"", edge.font[1],
                "\" fontsize=", edge.font[2], "]")
            Dot("  center=1")
            if (!is.null(min.rank)) {
                min.rank <- paste0("\"", min.rank, "\"")
                min.rank <- gsub(",", "\" \"", gsub(" ", "", min.rank))
                Dot("  {rank=min ", min.rank, "}", semicolon = FALSE)
            }
            if (!is.null(max.rank)) {
                max.rank <- paste0("\"", max.rank, "\"")
                max.rank <- gsub(",", "\" \"", gsub(" ", "", max.rank))
                Dot("  {rank=max ", max.rank, "}", semicolon = FALSE)
            }
            if (!is.null(same.rank)) {
                for (s in 1:length(same.rank)) {
                    same <- paste0("\"", same.rank[s], "\"")
                    same <- gsub(",", "\" \"", gsub(" ", "", same))
                    Dot("  {rank=same ", same, "}", semicolon = FALSE)
                }
            }
            latent <- variables[-(1:model$n)]
            for (lat in latent) {
                Dot("  \"", lat, "\" [shape=ellipse]", semicolon = FALSE)
            }
            endogenous <- classifyVariables(model$semmod)$endogenous
            endogenous <-
                variables[apply(outer(endogenous, model$var.names, "=="), 1, which)]
            if (style == "traditional") {
                variables <- c(variables, paste0(endogenous, ".error"))
                error.color <-
                    if (length(node.colors) < 3)
                        node.colors[1]
                else
                    node.colors[3]
            }
            for (endog in endogenous) {
                Dot("  \"", endog, "\" [fillcolor=\"", node.colors[2], "\"]", semicolon =
                        FALSE)
                if (style == "traditional") {
                    if (error.nodes) Dot(
                        "  \"", endog, ".error\" [shape=ellipse] [fillcolor=\"", error.color, "\"]",
                        semicolon = FALSE
                    )
                    else Dot(
                        "  \"", endog, 
                        ".error\" [shape=ellipse width=0 height=0 fixedsize=true label=\"\"] [fillcolor=\"", 
                        error.color, "\"]",
                        semicolon = FALSE
                    )
                }
            }
            ram <- model$ram
            if (missing(parameters)) {
                par.names <- names(coef(model))
                rownames(ram)[ram[, "parameter"] != 0] <-
                    par.names[ram[, "parameter"]]
                rownames(ram)[ram[, "parameter"] == 0] <-
                    ram[ram[, "parameter"] == 0, "start value"]
                parameters <- rownames(ram)
            }
            if (standardize)
                ram[, 5] <- stdCoef(model)[, 2]
            else
                ram[names(model$coeff), 5] <- model$coeff
            coefs <- ram[, 5]  # handle equality constraints, if any
            na.coefs <- is.na(coefs)
            if (any(na.coefs)) {
                for (coef in which(na.coefs)) {
                    ram[coef, 5] <-
                        (ram[ram[coef, 4] == ram[, 4], 5])[1] # paste in the estimate
                }
            }
            values <- round(ram[, 5], digits)
            heads <- ram[, 1]
            to <- ram[, 2]
            from <- ram[, 3]
            if (!missing(par.labels)) {
                check <- names(par.labels) %in% parameters
                if (any(!check)) {
                    msg <- if (sum(!check) > 1)
                        paste(
                            "The following parameters do not appear in the model:", paste(names(par.labels)[!check], collapse =
                                    ", ")
                        )
                    else
                        paste("The following parameter does not appear in the model:", names(par.labels)[!check])
                    warning(msg)
                    par.labels <- par.labels[check]
                }
                names(parameters) <- parameters
                parameters[names(par.labels)] <- par.labels
            }
            labels <- if (edge.labels == "names")
                parameters
            else if (edge.labels == "values")
                values
            else
                paste(parameters, values, sep = "=")
            colors <- ifelse(values > 0, edge.colors[1], edge.colors[2])
            direction <- ifelse((heads == 2), " dir=both", "")
            lineweight <- rep(1, nrow(ram))
            if (edge.weight == "proportional") {
                lineweight <- abs(values) / mean(values)
                if (!standardize)
                    warning("proportional edge weights requested for an unstandardized model")
            }
            if (style == "ram") {
                for (par in 1:nrow(ram)) {
                    if ((!ignore.double) || (heads[par] == 1)) {
                        if (ignore.self && to[par] == from[par]) next
                        Dot(
                            "  \"", variables[from[par]],
                            "\" -> \"", variables[to[par]], "\" [label=\"",
                            labels[par], "\"", direction[par], " color=", colors[par],
                            " penwidth=", round(lineweight[par] + 0.001, 3), "]"
                        )
                    }
                }
            }
            else
                for (par in 1:nrow(ram)) {
                    # style == "traditional"
                    if (heads[par] == 1) {
                        Dot(
                            "  \"", variables[from[par]],
                            "\" -> \"", variables[to[par]], "\" [label=\"",
                            labels[par], "\"", direction[par], " color=", colors[par],
                            " penwidth=", round(lineweight[par] + 0.001, 3), "]"
                        )
                    }
                    else if (variables[to[par]] %in% endogenous) {
                        if (variables[to[par]] == variables[from[par]]) {
                            # convert self-arrow to residual path
                            lab <- labels[par]
                            val <-
                                round(sqrt(values[par]), digits = digits)
                            lab <-
                                if (edge.labels == "names")
                                    paste0("sqrt(", lab, ")")
                            else if (edge.labels == "values")
                                val
                            else
                                paste0("sqrt(", parameters[par], ")=", val)
                            Dot(
                                "  \"", variables[to[par]], ".error\" -> \"",
                                variables[to[par]], "\" [color=", edge.colors[1], " label=\"", lab,
                                "\" penwidth=", round(sqrt(lineweight[par]) + 0.001, 3)," ]"
                            )
                        }
                        else{
                            # convert endogenous covariance to error covariance
                            Dot(
                                "  \"", variables[to[par]], ".error\" -> \"",
                                variables[from[par]], ".error\" [dir=both label=\"", labels[par],
                                "\" color=", colors[par],
                                " penwidth=", round(lineweight[par] + 0.001, 3), "]"
                            )
                        }
                    }
                    else if (!ignore.double &&
                            (variables[to[par]] != variables[from[par]])) {
                                Dot(
                                    "  \"", variables[from[par]],
                                    "\" -> \"", variables[to[par]], "\" [label=\"",
                                    labels[par], "\"", direction[par], " color=", colors[par],
                                    " penwidth=", round(lineweight[par] + 0.001, 3), "]"
                                )
                            }
                }
            if (!missing(var.labels)) {
                check <- names(var.labels) %in% variables
                if (any(!check)) {
                    msg <- if (sum(!check) > 1)
                        paste(
                            "The following variables do not appear in the model:", paste(names(var.labels)[!check], collapse =
                                    ", ")
                        )
                    else
                        paste("The following variable does not appear in the model:", names(var.labels)[!check])
                    warning(msg)
                    var.labels <- var.labels[check]
                }
                Dot("  // variable labels: ", semicolon = FALSE)
                lines <-
                    paste0('    "', names(var.labels), '" [label="', var.labels, '"];\n')
                Dot(paste(lines, collapse = ""), semicolon = FALSE, newline = FALSE)
            }
            Dot("}", semicolon = FALSE)
            if (output.type == "graphics") {
                cmd <-
                    paste0("dot -T", graphics.fmt, " -o ", graph.file,  " -Gcharset=latin1 ",
                        dot.options, " ", dot.file)
                cat("Running ", cmd, "\n")
                result <- try(system(cmd))
            }
            if (output.type == "html" && requireNamespace("DiagrammeR")) {
                print(DiagrammeR::DiagrammeR(textConnection(dot), type = "grViz"))
            }
            result <-
                if (output.type == "html")
                    dot
            else
                readLines(dot.file)
            invisible(result)
        }

