# Display debug message before evaluating a code block
ok_group <- function (message, tests = NULL) {
    cat(paste0("# ", unlist(strsplit(message, "[\r\n]+")), "\n", collapse=""), sep = "")
    tests
    invisible(NULL)
}
