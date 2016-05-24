sourceDir <-
function(absPath, trace = TRUE, ...)
{
    for (funcName in list.files(absPath, pattern = "\\.[R]$"))
      {
       source(file.path(absPath, funcName), ...)
       if(trace){cat("    Loaded:", funcName, "\n")}
      }
}
