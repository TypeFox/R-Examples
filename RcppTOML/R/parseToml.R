
parseTOML <- function(filename, verbose=FALSE) {
    fullfile <- path.expand(filename)
    toml <- tomlparseImpl(fullfile, verbose)
    class(toml) <- c("toml", "list")
    attr(toml, "file") <- filename
    toml
}

## alias for now, to be renamed
tomlparse <- function(...) parseTOML(...)

## alias for now, to be renamed
parseToml <- function(...) parseTOML(...)


print.toml <- function(x, ...) {
    print(utils::str(x, give.attr=FALSE))           # convenient shortcut
    #klass <- oldClass(x)
    #oldClass(x) <- klass[klass != "toml"]
    #NextMethod("print")
    invisible(x) 
}

summary.toml <- function(object, ...) {
    cat("toml object with top-level slots:\n")
    cat("  ", paste(names(object), collapse=", "), "\n")
    cat("read from", attr(object, "file"), "\n")
    invisible(NULL)
}
