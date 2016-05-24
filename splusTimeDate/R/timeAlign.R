"timeAlign" <- 
function(x, by = "days", k.by = 1, direction = 1, week.align = NULL, holidays
	 = timeDate())
{
  ## direction is +1 for after, -1 for before
  k.by <- as(k.by, "integer")
  if(k.by < 1)
    stop("k.by must be >= 1")
  ## first take care of a couple of special cases
  switch(by,
         weekdays = k.by <- 1,
         bizdays = k.by <- 1,
         weeks = {
           k.by <- 1
           if(is.null(week.align))
             by <- "days"
         }
         )
  k.by <- k.by * direction
  rtto <- timeRelative(by = by, k.by = k.by, align.by = TRUE, week.day = 
                       week.align, holidays. = holidays)
  rtaway <- timeRelative(by = by, k.by =  - k.by, align.by = TRUE, week.day
                         = week.align, holidays. = holidays)
  ## go away from and then back to correct direction
  rtto@Data <- paste(rtaway@Data, rtto@Data)
  x + rtto
}

