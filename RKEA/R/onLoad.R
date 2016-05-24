.onLoad <-
function(libname, pkgname)
{
    ## Create and populate RKEA_work_dir to avoid the warnings resulting
    ## from KEA's hard-coding of stopword file paths.
    RKEA_work_dir <<- tempfile()
    dir.create(file.path(RKEA_work_dir, "data", "stopwords"),
               recursive = TRUE)
    file.copy(Sys.glob(file.path(system.file("stopwords",
                                             package = "RKEAjars"),
                                 "stopwords_*.txt")),
              file.path(RKEA_work_dir, "data", "stopwords"))
}
