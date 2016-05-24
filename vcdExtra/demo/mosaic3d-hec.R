## mosaic3d-hec	2D and 3D visualizations of HairEyeColor data

library(vcdExtra)

# two-way mosaic displays
HairEye <- margin.table(HairEyeColor, c(1,2))
HairEye

mosaic(HairEye, shade=TRUE)

# three-way mosaic displays
structable(HairEyeColor)

# mutual independence model
mosaic(HairEyeColor, shade=TRUE)
# joint independence of Hair*Eye with Sex
mosaic(HairEyeColor, expected =~(Hair*Eye)+Sex)

# observed frequencies, mutual independence
mosaic3d(HairEyeColor)

# expected frequencies under mutual independence
mosaic3d(HairEyeColor, type="expected")

# expected frequencies under joint independence
mosaic3d(HairEyeColor, type="expected", expected =~(Hair*Eye)+Sex)


