# We used to call the internal seek() function directly, because it had
# less overhead (almost twice as fast) compared to seek(con=con, 
# where=offset, origin="start", rw="read").  However, since R v2.15.0
# .Internal() calls are no longer allowed. /HB 2012-04-12
##  origin <- pmatch("start", c("start", "current", "end"));
##  rw <- pmatch("read", c("read", "write"), 0);
##  .Internal(seek(con, offset, origin, rw));

.seekCon <- function(con, where, rw="read", ...) {
  base::seek.connection(con=con, where=where, origin="start", rw=rw);
} # .seekCon()


############################################################################
# HISTORY:
# 2012-04-15
# o Added .seekCon() to avoid using .Internal(seek(...)).
# o Created.
############################################################################
