require("TSmisc")
require("xts")
require("tfplot")

pit <- TSconnect("zip", dbname="http://pitrading.com/free_eod_data")

indu <- TSget("INDU", con=pit, TSrepresentation=xts, quote="Close")

tmp <- tempdir()

files <- paste(tmp, c("indusmall.png", "vixsmall.png",
    "tyxsmall.png", "vsUSDsmall.png", "adcdsmall.png" ),sep="/")

png(file=files[1],width=480, height=240, pointsize=12, bg = "white")
# mv indusmall.png indusmall.png.orig ; pngcrush indusmall.png.orig indusmall.png
#png(file="indu.png",    width=960, height=480, pointsize=12, bg = "white")
  tfOnePlot(indu, start=as.Date("2011-01-01"), 
    Title = "Running commentary, blah, blah, blah", 
    subtitle="Dow Jones Industrial Average Index Close",
    ylab= "index",
    xlab= "2011",
    source="Source: pitrading (indu)",
    footnoteRight = paste("Extracted:", date()),
    lastObs = TRUE )
dev.off()


vix <- TSget("VIX", con=pit, TSrepresentation=xts, quote="Close")

png(file=files[2],width=480, height=240, pointsize=12, bg = "white")
# mv vixsmall.png vixsmall.png.orig ; pngcrush vixsmall.png.orig vixsmall.png
#png(file="vix.png",    width=960, height=480, pointsize=12, bg = "white")
  tfOnePlot(vix, start=as.Date("2011-01-01"), 
    Title = "Running commentary, blah, blah, blah", 
    subtitle="CBOE Volatility Index - Close",
    ylab= "index",
    xlab= "2011",
    source="Source: pitrading (vix)",
    footnoteRight = paste("Extracted:", date()),
    lastObs = TRUE )
dev.off()


tyx <- TSget("TYX", con=pit, TSrepresentation=xts, quote="Close")

# not sure what this is, the scale seems off for tyx (compare TShistQuote)
png(file=files[3],width=480, height=240, pointsize=12, bg = "white")
# mv tyxsmall.png tyxsmall.png.orig ; pngcrush tyxsmall.png.orig tyxsmall.png
#png(file="tyx.png",    width=960, height=480, pointsize=12, bg = "white")
  tfOnePlot(tyx, start=as.Date("2011-01-01"), 
    Title = "Running commentary, blah, blah, blah", 
    subtitle="30-Year Treasury Bond",
    ylab= "index",
    xlab= "2011",
    source="Source: pitrading (tyx)",
    footnoteRight = paste("Extracted:", date()),
    lastObs = TRUE )
dev.off()


vsUSD <- TSget(c("EURUSD", "GBPUSD"), con=pit, quote="Close")
png(file=files[4],width=480, height=240, pointsize=12, bg = "white")
# mv vsUSDsmall.png vsUSDsmall.png.orig ; pngcrush vsUSDsmall.png.orig vsUSDsmall.png
#png(file="vsUSD.png",    width=960, height=480, pointsize=12, bg = "white")
  tfOnePlot(vsUSD, start=as.Date("2011-01-01"), 
    Title = "Running commentary, blah, blah, blah", 
    subtitle="Euro /US$ (black) British Pound / US$ (red)  - Close",
    ylab= "per US$",
    xlab= "2011",
    source="Source: pitrading (EURUSD, GBPUSD)",
    footnoteRight = paste("Extracted:", date()),
    lastObs = TRUE )
dev.off()

adcd <- TSget(c("AD", "CD"), con=pit, quote="Close")

png(file=files[5],width=480, height=240, pointsize=12, bg = "white")
# mv adcdsmall.png adcdsmall.png.orig ; pngcrush adcdsmall.png.orig adcdsmall.png
#png(file="adcd.png",    width=960, height=480, pointsize=12, bg = "white")
  tfOnePlot(adcd, start=as.Date("2011-01-01"), 
    Title = "Australian and Canadian Dollar Futures", 
    subtitle="Australian (black) and Canadian (red) - Close",
    ylab= "",
    xlab= "2011",
    source="Source: pitrading (AD, CD)",
    footnoteRight = paste("Extracted:", date()),
    lastObs = TRUE )
dev.off()

unlink(tmp, recursive = TRUE)

