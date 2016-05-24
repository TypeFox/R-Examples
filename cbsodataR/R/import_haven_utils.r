import_haven_utils <- function(){
  file <- "R/util.R"
  cat("# Borrow functions from haven, automatically created by 'import_haven_utils.R'\n", file=file)
  
  dump( list=c( "labelled"
                , "as_factor"
                , "as_factor.labelled"
                , "as_factor.character"
                , "[.labelled"
  )
  , file   = file
  , envir  = getNamespace('haven')
  , append = TRUE
  )
}
