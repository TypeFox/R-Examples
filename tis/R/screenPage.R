screenPage <- function(head = NULL,
                       sub = NULL,
                       foot = NULL,
                       date = FALSE,
                       dateFormat = "%x",
                       time = FALSE,
                       topLeft = character(0),
                       topRight = character(0),
                       headFont = par("font.main"),
                       subFont = par("font.sub"),
                       footFont = par("font"),
                       cex = 1.001,
                       headCex = 1.5,
                       subCex = 0.85,
                       footCex = 0.75,
                       topLeftCex = 0.85,
                       topRightCex = 0.85,
                       footAlign = 0,
                       leftMargin = 0,
                       rightMargin = leftMargin,
                       topMargin = 0,
                       bottomMargin = topMargin){
  ## Calculate outer margin space needed for outer margin text items: titles,
  ## subtitles, footnotes, and various labels.  Set par(oma) appropriately,
  ## then draw the outer margin items.
  
  value <- list(head = head, sub = sub, foot = foot,
                date = date, dateFormat = dateFormat, time = time,
                topLeft = topLeft, topRight = topRight,
                headFont = headFont, subFont = subFont, footFont = footFont,
                cex = cex, headCex = headCex, subCex = subCex,
                footCex = footCex,
                topLeftCex = topLeftCex, topRightCex = topRightCex,
                footAlign = footAlign, 
                leftMargin = leftMargin, rightMargin = rightMargin,
                topMargin = topMargin, bottomMargin = bottomMargin)
  
  if(length(dev.list())== 0) ## no graphics device is active
    return(value)

  plot.new()
  ## unsplit screen
  if(length(ss <- split.screen()) != 1 || ss) close.screen(all.screens = TRUE)

  par(mex = 1, cex = 1, fig = c(0, 1, 0, 1),
      omi = c(bottomMargin, leftMargin, topMargin, rightMargin))
  
  if(!is.null(head))    head <- evalOrEcho(head)
  if(!is.null( sub))    sub  <- evalOrEcho(sub)
  if(!is.null(foot))    foot <- evalOrEcho(foot)

  footAlign <- rep(footAlign, length.out = length(foot))
  
  ## All cex's are relative to cex
  topLeftCex   <- cex * topLeftCex
  topRightCex  <- cex * topRightCex
  headCex <- cex * headCex
  subCex  <- cex * subCex
  footCex <- cex * footCex
  ## topRight (Upper Right Corner text) will hold date and time
  topRight <- character(0)
  if(is.character(date))
    topRight <- c(topRight, date)
  else
    if(date) topRight <- c(topRight, format(Sys.time(), dateFormat))
  if(time) topRight <- c(topRight, format(Sys.time(), "%X"))
  ## Figure top margin space
  nHead <- length(head)
  nSub <- length(sub)
  nTopLeft  <- length(topLeft)
  nTopRight <- length(topRight)
  topLeftSpace  <- topLeftCex * (nTopLeft + 1)
  topRightSpace <- topRightCex * (nTopRight + 1)
  subSpace <- subCex * nSub
  headSpace <- headCex * nHead
  topSpace <- max(max(topLeftSpace, topRightSpace) + headSpace + subSpace + 0.5)
  ## Figure additional bottom margin space
  nFoot <- length(foot)
  bottomSpace <- footCex * (nFoot + 1)

  ## Set par("oma")
  oma <- par("oma")
  par(oma = c(bottomSpace + oma[1], oma[2], topSpace + oma[3], oma[4]))
  
  ## place topLeft
  if(nTopLeft)
    for(i in 1:nTopLeft)
      mtext(text = topLeft[i], side = 3, adj = 0, outer = TRUE, cex = topLeftCex,
            line = topSpace - i * topLeftCex)
  
  ## place topRight
  if(nTopRight)
    for(i in 1:nTopRight)
      mtext(text = topRight[i], side = 3, adj = 1, outer = TRUE, cex = topRightCex,
            line = topSpace - i * topRightCex)
  ## place head
  if(nHead)
    for(i in 1:nHead)
      mtext(head[i], side = 3, adj = 0.5, outer = TRUE, cex = headCex,
            line = headSpace + subSpace - i * headCex + 0.5,
            font = headFont)
  ## place sub
  if(nSub)
    for(i in 1:nSub)
      mtext(sub[i], side = 3, adj = 0.5, outer = TRUE, cex = subCex,
            line = subSpace - i * subCex + 0.5, font = subFont)
  ## place foot
  if(nFoot)
    for(i in 1:nFoot)
      mtext(text = foot[i], side = 1, adj = footAlign[i],
            outer = TRUE, cex = footCex, line = (i - 1) * footCex)

  omd <- par("omd")
  ## omd is the boundaries of the
  ## c(left, right, bottom, top) outer margins, as fractions (in [0,1])
  ## of the device surface.
  omd[3] <- max(omd[3], 0.002)
  pagehead <- c(omd[1], omd[2], omd[4], 1)
  pagefoot <- c(omd[1], omd[2], 0,      omd[3])
  pagebody <- c(omd[1], omd[2], omd[3], omd[4])
  ## reset outer margins to zero
  ## par(oma = c(0, 0, 0, 0))
  ## figs:c(left, right, bottom,  top )
  if(pagefoot[4] > 0.001)
    figs <- rbind(pagehead, pagefoot, pagebody)
  else
    figs <- rbind(pagehead, pagebody)
  split.screen(figs, erase = FALSE)
  screens <- split.screen()
  screen(screens[length(screens)], new = TRUE)
  ## collecting garbage seems to solve problem of mysterious crashes
  gc()
  invisible(value)
}

