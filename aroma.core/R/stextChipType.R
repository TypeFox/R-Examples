setMethodS3("stextChipType", "character", function(chipType, side=4, fmtstr="%s", pos=1, cex=0.7, col="darkgray", ...) {
  stext(side=side, text=sprintf(fmtstr, chipType), pos=pos, cex=cex, col=col, ...);
}, private=TRUE)


############################################################################
# HISTORY:
# 2008-05-17
# o Created from ditto for AffymetrixCdfFile:s.
############################################################################
