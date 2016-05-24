wqData <-
function(data, locus, wqdata, site.order, time.format = "%Y-%m-%d",
  type = c("long", "wide")) {    

  # Validate args
  if (length(locus) != 3)
      stop("locus must be of length 3")
  cnames <- colnames(data)
  if (is(wqdata, "character"))
      wqdata <- match(wqdata, cnames, nomatch=0)
  if (any(identical(wqdata, 0)) || max(wqdata) > ncol(data))
      stop("wqdata not in data")
  type <- match.arg(type)

  # Reshape data
  if (identical(type, "long")) {
      data <- data.frame(data[, locus], data[, wqdata])
      names(data) <- c("time", "site", "depth", "variable", "value")
  } else {
      if (identical(length(wqdata), 1L)) {		
          data <- data.frame(data[, locus], variable =
            rep(cnames[wqdata], nrow(data)), value = data[, wqdata])
          names(data)[1:3] <- c("time", "site", "depth")
      } else {
          # Avoid possible duplicate names
          wqd <- data[, wqdata]
          ind <- match(c("time", "site", "depth"), names(wqd), nomatch=0)
          names(wqd)[ind] <- paste(names(wqd)[ind], 1, sep="")
          # Assemble and reshape data
          data <- data.frame(data[, locus], wqd)
          names(data)[1:3] <- c("time", "site", "depth")
          data <- melt(data, id.vars = 1:3)
      }
  }
  
  # Change time to correct format and class if needed
  if (grepl('H', time.format)) {
    data <- within(data, time <- as.POSIXct(time, format =
      time.format))
  } else {
      data <- within(data, time <- as.Date(time, format =
        time.format))
  }
          
  # Remove NAs
  data <- data[!is.na(data$value), ]
  rownames(data) <- 1:nrow(data)
  
  # Remove unneeded factor levels
  data <- within(data, site <- factor(site, ordered = site.order))
  levels(data$site) <- gsub('X','s', make.names(levels(data$site),
    unique = TRUE))
  
  # Make sure variable is a factor
  data <- within(data, variable <- as.factor(variable))

  #
  new(Class="WqData", data)

}
