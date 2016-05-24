# example spcosa package: projection

# check if required packages are available
if (suppressWarnings(!(require(grid) & require(sp) & require(ggplot2)))) {
    stop("This demo requires packages 'grid', 'sp', and 'ggplot2'.\nMake sure these packages are installed", call. = FALSE)
}    

# construct two identical grids
grd1 <- expand.grid(
    longitude = 1:50,
    latitude  = 90 - 1:50
)
coordinates(grd1) <- ~ longitude * latitude
gridded(grd1) <- TRUE
grd2 <- grd1

# add projection attributes to grd2
proj4string(grd2) <- CRS("+proj=longlat +ellps=WGS84")

# stratify both grids into 50 strata
strat1 <- stratify(grd1, nStrata = 50)
strat2 <- stratify(grd2, nStrata = 50)

# plot resulting strata
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(
    plot(strat1) + ggtitle("without projection attributes\n"),
    vp = viewport(layout.pos.row = 1, layout.pos.col = 1)
)
print(
    plot(strat2) + ggtitle("with projection attributes\n"),
    vp = viewport(layout.pos.row = 1, layout.pos.col = 2)
)
