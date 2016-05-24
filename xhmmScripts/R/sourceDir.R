sourceDir <- function(path, trace = FALSE, recursive = FALSE, ...) {
        for (nm in list.files(path, pattern = "\\.[RrSsQq]$", recursive = recursive)) {
                if(trace) cat(nm,":")           
                fp = file.path(path, nm)
                source(fp, ...);
                if(trace) cat("\n")
        }
}
