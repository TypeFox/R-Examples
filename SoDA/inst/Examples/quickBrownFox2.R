txt <- "The quick brown fox"

require(grid)

{if(interactive())
  grid.newpage()
else
  pdf("Examples/quickBrownFox2.pdf", width = 4, height=1.5)
}

pushViewport(viewport(gp = gpar(fontsize=18)))
grid.rect(width=unit(1, "strwidth", txt),
                      height=unit(1,   "lines"), gp = gpar(lwd=2))
grid.text(txt)
