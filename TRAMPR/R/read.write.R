## read.write.R (part of the TRAMPR package)

## TODO: This is something of a fudge for now, but will do.  I might
## come up with nicer I/O later.
write.TRAMPknowns <- function(x, file.pat=x$file.pat, warn=TRUE) {
  valid.TRAMPknowns(x)
  if ( !is.null(file.pat) ) {
    write.csv(x$info, sprintf("%s_info.csv", file.pat),
              row.names=FALSE)
    write.csv(x$data, sprintf("%s_data.csv", file.pat),
              row.names=FALSE)
  } else if ( warn ) {
    warning("No file pattern specified and x$file.pat is NULL")
  }
}

read.TRAMPknowns <- function(file.pat, auto.save=TRUE, overwrite=FALSE) {

  f.info <- sprintf("%s_info.csv", file.pat)
  f.data <- sprintf("%s_data.csv", file.pat)
  info <- read.csv.safe(f.info)
  data <- read.csv.safe(f.data)
  if ( auto.save )
    file.copy(c(f.info, f.data),
              sprintf("%s_%s_%s.csv", file.pat, c("info", "data"),
                      format(Sys.Date(), "%Y%m%d")), overwrite)
  else
    file.pat <- NULL
  TRAMPknowns(data, info, file.pat=file.pat)
}

write.TRAMPsamples <- function(x, file.pat) {
  valid.TRAMPsamples(x)
  write.csv(x$info, sprintf("%s_info.csv", file.pat),
            row.names=FALSE)
  write.csv(x$data, sprintf("%s_data.csv", file.pat),
            row.names=FALSE)
}

read.TRAMPsamples <- function(file.pat) {
  f.info <- sprintf("%s_info.csv", file.pat)
  f.data <- sprintf("%s_data.csv", file.pat)
  info <- read.csv.safe(f.info)
  data <- read.csv.safe(f.data)
  TRAMPsamples(data, info)
}
