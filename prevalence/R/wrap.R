wrap <-
function(me){
  # 1e-06 => 1.0E-06
  me <- gsub("e+", "E+", me, fixed = TRUE)
  me <- gsub("e-", "E-", me, fixed = TRUE)

  for (line in seq_along(me)){
    s <- me[line]
    split <- strsplit(s, "[0-9]|[0-9]\\.[0-9]|\\.[0-9]")[[1]]

    for (i in seq_along(split))
      if (split[i] != "E+" & split[i] != "E-" & split[i] != "")
        s <- paste(strsplit(s, split[i], fixed = TRUE)[[1]], collapse = "@")

    split <- unlist(strsplit(s, "@"))

    for (i in seq_along(split))
      if (regexpr("E+", split[i])[1] != -1){
        dig <- nchar(strsplit(split[i],"E")[[1]][1])
        me[line] <- sub(split[i],
          formatC(as.numeric(split[i]), format = "E", digits = dig),
          me[line], fixed = TRUE)
      }
  }

  write(me, "modelTempFile.txt")
}