library(colorspace)
blueShade <- function(inrgb) {
    rgb <- col2rgb(inrgb)
    RGB <- RGB(t(rgb)/255) 
    # Special case "black"
    if (all(coords(RGB) == 0))
        RGB <- RGB(0, 0, .1)
    LCH <- as(RGB, "polarLUV")
    lch <- coords(LCH)
    # Scale the chroma so greys become blues
    hcl(240, 20 + .8*lch[2], lch[1])
}
