rtconnect <-
function(daily.dir="~/data/daily") {
  daily.files <- list.files(daily.dir)
  daily.files <- paste(daily.dir, "/", daily.files, sep="")

  data <- NULL
  for (f in daily.files) {
    df <- read.table(f, sep="\t", header=T)
    # 2013/10頃以前のレポートについてはCategoryが存在しないため
    # Category列を削除して処理する
    df$Category <- NULL
    data <- rbind(data, df)
  }

  data$Version <- as.character(data$Version)
  data$Product.Type.Identifier <- as.character(data$Product.Type.Identifier)
  data$Begin.Date <- as.Date(data$Begin.Date, "%m/%d/%Y")
  data$End.Date <- as.Date(data$End.Date, "%m/%d/%Y")

  class(data) <- c("rtconnect", "data.frame")

  return(invisible(data))
}
