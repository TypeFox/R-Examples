## trimws was added in R 2.3.0.  If we're using a previous version of
## R we need to define it.
if(!exists('trimws', mode='function'))
   trimws <-  function(s)
                {
                  s <- sub(pattern="^[[:blank:]]+", replacement="", x=s)
                  s <- sub(pattern="[[:blank:]]+$", replacement="", x=s)
                  s
                }

   
