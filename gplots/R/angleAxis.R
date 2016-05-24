angleAxis <- function(side, labels, at=1:length(labels), srt=45, adj, xpd=TRUE, ...)
{
  usr <- par("usr")
  emH <- strheight("M")
  emW <- strwidth("M")
  
  if(missing(adj))
    adj <- switch(side, 1, 1, 0, 0)
  
  switch(side,
         #1 - below
         text(x=at,           y=usr[3]-emH/2, labels=labels, srt=srt, adj=adj, xpd=xpd, ...), 
         #2 - left
         text(x=usr[1]-emW/2, y=at,           labels=labels, srt=srt, adj=adj, xpd=xpd, ...), 
         #3 - above
         text(x=at,           y=usr[4]+emH/2, labels=labels, srt=srt, adj=adj, xpd=xpd, ...), 
         #4 - right
         text(x=usr[2]+emW/2, y=at,           labels=labels, srt=srt, adj=adj, xpd=xpd, ...)
  )       
  
  invisible(NULL)
}