basicStyles <- getStyleDefs()

basicStyles$lowerBorder$backgroundColor <-  "#FFFFFF"

basicStyles$lowerBorderGray <-  basicStyles$lowerBorder
basicStyles$lowerBorderGray$backgroundColor <-  "#E8E8E8"

basicStyles$ArialCenteredItalics <- basicStyles$ArialCenteredBold
basicStyles$ArialCenteredItalics$fontSize <- "10pt"
basicStyles$ArialCenteredItalics$fontType <- "italic"

basicStyles$ArialHighlight$fontColor <- "#000080"
basicStyles$ArialHighlight$fontType <- ""

basicStyles$highlight <- basicStyles$noBorder
basicStyles$highlight$backgroundColor <- "#D4EBF3"

basicStyles$RTable2 <- basicStyles$RTable1
basicStyles$RTable2$marginLeft <- "1.00in"
basicStyles$RTable2$marginRight <- "1.00in"
basicStyles$RTable2$marginTop <- "0.05in"
basicStyles$RTable2$marginBottom <- "2in"

basicStyles$frameWithBorders <- basicStyles$basicFigFrame
basicStyles$frameWithBorders$bottomBorder <- "0.01in solid #000000"
basicStyles$frameWithBorders$topBorder <- "0.01in solid #000000"
basicStyles$frameWithBorders$leftBorder <- "0.01in solid #000000"
basicStyles$frameWithBorders$rightBorder <- "0.01in solid #000000"

basicStyles$paddingEx <- basicStyles$frameWithBorders
basicStyles$paddingEx$padding <- "0.5in"

basicStyles$anchor2 <- basicStyles$basicFigFrame
basicStyles$anchor2$imageAnchor <- "as-char"

basicStyles$anchor1 <- basicStyles$basicFigFrame
basicStyles$anchor1$frameAnchor <- "as-char"

basicStyles$wrapping <- basicStyles$basicFigFrame
basicStyles$wrapping$wrap <- "right"

basicStyles$wideBullet <- basicStyles$Rbullet
basicStyles$wideBullet$spaceBefore <- "0.1in"
basicStyles$wideBullet$minLabelWidth <- "0.5in"
basicStyles$wideBullet$paraStyle <- "ttBlue"
basicStyles$wideBullet$bulletChar <- "\342\234\224"
Encoding(basicStyles$wideBullet$bulletChar) <- "UTF-8"

setStyleDefs(basicStyles)

plotInfo <- getImageDefs()
plotInfo$dispHeight <- 4  
plotInfo$dispWidth <- 4 
setImageDefs(plotInfo)

options(width = 80)

odfWeave("formatting.odt", "formattingOut.odt")
