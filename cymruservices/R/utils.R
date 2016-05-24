# remove leading/trailing blanks around all strings in data.frame cols
trim_df <- function(df, stringsAsFactors=FALSE) {
  data.frame(lapply(df, function (v) {
    if (is.character(v)) {
      trimws(v)
    } else {
      v
    }
  }), stringsAsFactors=stringsAsFactors)
}
