"getExt" <-
  function(x)
{
  hasPeriod <-  grep("\\.", x)
  if(length(x) == 1 & length(hasPeriod) == 0)
    {
      z <- NA
    } else {
      if(length(hasPeriod) < length(x)) x[-hasPeriod] <- NA
      y <- strsplit(x, ".", fixed = TRUE)
      z <- lapply(y, function(data) data[length(data)])
    }
  unlist(z)
}

"getExt.new" <-
  function(x)
{
  hasPeriod <-  grep("\\.", x)
  if(length(hasPeriod) == 0) {
    z <- rep(NA, length(x))
  } else {
    if(length(hasPeriod) < length(x))
      x[-hasPeriod] <- NA
    y <- strsplit(x, ".", fixed = TRUE)
    z <- lapply(y, function(data)
                { if (length(data) > 1)
                    {
                      data[length(data)]
                    } else NA
                }
                )
  }
  unlist(z)
}


## Test it - last two examples work only for the new version
## getExt.new(c("a.txt", "b.R"))
## getExt.new("a.txt")
## getExt.new("a")
## getExt.new(c("a.txt", "b"))
## getExt.new(c("a", "b"))
## getExt.new(c("a.txt", "b."))




