## simple script for converting C++ Qt code to R Qt code

qtranslate <- function(f = "tmp.cpp") {
  lines <- readLines(f)
  table <- c("new " = "Qt$", "->" = "$", ".* \\*" = "", ";" = "", "=" = "<-",
             "::" = "$", "tr\\((.*?)\\)" = "\\1", "true" = "TRUE",
             "false" = "FALSE")
  for(i in seq_along(table))
    lines <- gsub(names(table)[i], table[i], lines)
  rf <- sub("cpp$", "R", f)
  writeLines(lines, rf)
}
