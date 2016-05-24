# Chemical reactor data Table 7.6 of Myers et al. (2009), 
# Response Surface Methodology (3rd ed.), Wiley.

#--- CCD in 2 blocks
ChemReact <-
structure(list(Time = c(80, 80, 90, 90, 85, 85, 85, 85, 85, 85, 
92.07, 77.93, 85, 85), Temp = c(170, 180, 170, 180, 175, 175, 
175, 175, 175, 175, 175, 175, 182.07, 167.93), Block = structure(c(1L, 
1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L), .Label = c("B1", 
"B2"), class = "factor"), Yield = c(80.5, 81.5, 82, 83.5, 83.9, 
84.3, 84, 79.7, 79.8, 79.5, 78.4, 75.6, 78.5, 77)), .Names = c("Time", 
"Temp", "Block", "Yield"), class = "data.frame", row.names = c("1", 
"2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", 
"14"))

#--- Just the first block
ChemReact1 <-
structure(list(Time = c(80, 80, 90, 90, 85, 85, 85), Temp = c(170, 
180, 170, 180, 175, 175, 175), Yield = c(80.5, 81.5, 82, 83.5, 
83.9, 84.3, 84)), .Names = c("Time", "Temp", "Yield"), row.names = c("1", 
"2", "3", "4", "5", "6", "7"), class = "data.frame")

#--- Just the 2nd block
ChemReact2 <-
structure(list(Time = c(85, 85, 85, 92.07, 77.93, 85, 85), Temp = c(175, 
175, 175, 175, 175, 182.07, 167.93), Yield = c(79.7, 79.8, 79.5, 
78.4, 75.6, 78.5, 77)), .Names = c("Time", "Temp", "Yield"), row.names = c("8", 
"9", "10", "11", "12", "13", "14"), class = "data.frame")


# heli data -- CCD with 4 factors in 2 blocks
# Table 12.5 of Box, GEP, Hunter, JS, and Hunter, WG (2005), 
# Statistics for Experimenters (2nd ed.), Wiley.
# NOTE -- subsequent statements translate it to coded form
heli <-
structure(list(block = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L), .Label = c("1", "2"), class = "factor"), 
    x1 = c(-1L, 1L, -1L, 1L, -1L, 1L, -1L, 1L, -1L, 1L, -1L, 
    1L, -1L, 1L, -1L, 1L, 0L, 0L, -2L, 2L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L), x2 = c(-1L, -1L, 1L, 1L, -1L, -1L, 1L, 
    1L, -1L, -1L, 1L, 1L, -1L, -1L, 1L, 1L, 0L, 0L, 0L, 0L, -2L, 
    2L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L), x3 = c(-1L, -1L, -1L, 
    -1L, 1L, 1L, 1L, 1L, -1L, -1L, -1L, -1L, 1L, 1L, 1L, 1L, 
    0L, 0L, 0L, 0L, 0L, 0L, -2L, 2L, 0L, 0L, 0L, 0L, 0L, 0L), 
    x4 = c(-1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, -2L, 
    2L, 0L, 0L, 0L, 0L), ave = c(367L, 369L, 374L, 370L, 372L, 
    355L, 397L, 377L, 350L, 373L, 358L, 363L, 344L, 355L, 370L, 
    362L, 377L, 375L, 361L, 364L, 355L, 373L, 361L, 360L, 380L, 
    360L, 370L, 368L, 369L, 366L), logSD = c(72L, 72L, 74L, 79L, 
    72L, 81L, 72L, 99L, 90L, 86L, 92L, 112L, 76L, 69L, 91L, 71L, 
    51L, 74L, 111L, 93L, 100L, 80L, 71L, 98L, 69L, 74L, 86L, 
    74L, 89L, 76L)), .Names = c("block", "x1", "x2", "x3", "x4", 
  "ave", "logSD"), row.names = c(NA, -30L), class = c("data.frame"))
# TRANSLATE TO CODED FORM
heli = as.coded.data(heli, x1 ~ (A - 12.4)/0.6, x2 ~ (R - 2.52)/0.26, 
            x3 ~ (W - 1.25)/0.25, x4 = x4 ~ (L - 2)/0.5)

# CO emissions data - BH^2 2nd ed., Table 10.17, p.396.
# See also sloping ridge discussion, pp 466-7
codata <-
structure(list(x1 = c(-1L, -1L, 0L, 0L, 1L, 1L, -1L, -1L, 0L, 0L, 1L, 1L, -1L, -1L, 0L, 0L, 1L, 1L), 
x2 = c(-1L, -1L, -1L, -1L, -1L, -1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L), 
y = c(61.9, 65.6, 80.9, 78, 89.7, 93.8, 72.1, 67.3, 80.1, 81.4, 77.8, 74.8, 66.4, 68.2, 68.9, 66, 60.2, 57.9)), 
.Names = c("x1", "x2", "y"), class = "data.frame", row.names = c(NA, -18L))

# coded.data version - created via
# CO <- as.coded.data(data = codata, x1 ~ (Ethanol - 0.2)/0.1, x2 ~ A.F.ratio - 15)
# names(CO)[3] = "CO.conc"
# CO <-
# structure(list(x1 = c(-1L, -1L, 0L, 0L, 1L, 1L, -1L, -1L, 0L, 
# 0L, 1L, 1L, -1L, -1L, 0L, 0L, 1L, 1L), x2 = c(-1L, -1L, -1L, 
# -1L, -1L, -1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L
# ), CO.conc = c(61.9, 65.6, 80.9, 78, 89.7, 93.8, 72.1, 67.3, 
# 80.1, 81.4, 77.8, 74.8, 66.4, 68.2, 68.9, 66, 60.2, 57.9)), .Names = c("x1", 
# "x2", "CO.conc"), class = c("coded.data", "data.frame"), row.names = c(NA, 
# -18L), codings = structure(list(x1 = quote(x1 ~ (Ethanol - 0.2)/0.1), 
# x2 = quote(x2 ~ A.F.ratio - 15)), .Names = c("x1", "x2")), rsdes = structure(list(
# primary = c("x1", "x2"), n0 = c(2L, NA), non0 = c(16L, NA
# ), alpha = "NA", blk.info = list(structure(list(n0 = c(2L, 
# NA), non0 = c(16L, NA), alpha = NA), .Names = c("n0", "non0", 
# "alpha"))), call = quote(as.coded.data(data = co, x1 ~ (Ethanol - 
# 0.2)/0.1, x2 ~ A.F.ratio - 15))), .Names = c("primary", 
# "n0", "non0", "alpha", "blk.info", "call")))


