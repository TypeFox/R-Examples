library(hexbin)

if(FALSE) { ## the following is still quite a bit from working/useful :

## what should that do? set a palette?
rgb <- matrix(c(
        15,15,15,

        0, 0, 0,
        1, 9,15,
        9,15, 9,
       15, 9, 9,

        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,

        9, 9, 9,
        0, 2, 7,
        0, 7, 1,
        8, 1, 1,

       15, 2, 2,
       11, 1, 1,
        8, 1, 1,
        5, 1, 1,
        5, 1, 1,
       15,15,15), ncol = 3, byrow = TRUE)

##ps.options(rasters=600,color=rgb/15,background=2)
##ps.options(color=rgb/15,background=2)
postscript("large.ps",width = 10,height = 7.5)

plot.hexbin(ans.25mil, style = "nest", lcex = .9)

}## FALSE, i.e. nothing done
