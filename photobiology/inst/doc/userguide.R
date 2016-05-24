## ----setup, include=FALSE, cache=FALSE-----------------------------------
library(knitr)
opts_chunk$set(fig.path = 'figure/pos-', fig.align = 'center', fig.show = 'hold',
               fig.width = 7, fig.height = 6, size = "footnotesize")
# options(replace.assign=TRUE,width=60)
options(warnPartialMatchArgs = FALSE)

## ----example-0-hiden, eval=TRUE, include=FALSE---------------------------
library(photobiology)
library(lubridate)

## ----own-set-up, echo=FALSE, include=FALSE-------------------------------
my_version <- packageVersion("photobiology")

## ----set-up-printing-----------------------------------------------------
options(dplyr.print_max = 4)
options(dplyr.print_min = 4)

## ----eval=FALSE----------------------------------------------------------
#  print(sun.spct, n = 3)

## ----eval=FALSE----------------------------------------------------------
#  summary(sun.spct)

## ----example-1, eval=FALSE-----------------------------------------------
#  # not run
#  my.spct <- source_spct(w.length = wavelength/10, s.e.irrad = irrad/1000)

## ------------------------------------------------------------------------
is.any_spct(sun.spct)
is.source_spct(sun.spct)

## ------------------------------------------------------------------------
class_spct(sun.spct)
class(sun.spct)

## ------------------------------------------------------------------------
my.df <- data.frame(w.length = 400:410, s.e.irrad = 1)
my.spct <- as.source_spct(my.df)
class(my.spct)
class(my.df)
my.spct

## ------------------------------------------------------------------------
my.g.spct <- as.generic_spct(my.spct)
class(my.g.spct)

## ------------------------------------------------------------------------
source_spct(w.length = 300:305, s.e.irrad = 1)

## ------------------------------------------------------------------------
z <- 300:305
y <- 2
source_spct(w.length = z, s.e.irrad = y)

## ------------------------------------------------------------------------
w.length <- 300:305
s.e.irrad <- 1
source_spct(w.length, s.e.irrad)

## ------------------------------------------------------------------------
my.d.spct <- as.source_spct(my.df, time.unit = "day")

## ------------------------------------------------------------------------
source_spct(w.length = 300:305, s.e.irrad = -1)
source_spct(w.length = 300:305, s.e.irrad = -1, strict.range = NULL)

## ------------------------------------------------------------------------
my.cm.spct <- source_spct(w.length = 300:305, s.e.irrad = 1,
                          comment = "This is a comment")
comment(my.cm.spct)

## ------------------------------------------------------------------------
my.spct <- sun.spct
setWhenMeasured(my.spct, NULL)
getWhenMeasured(my.spct)
setWhenMeasured(my.spct, ymd_hms("2015-10-31 22:55:00", tz = "EET"))
getWhenMeasured(my.spct)

## ------------------------------------------------------------------------
setWhereMeasured(my.spct, NULL)
getWhereMeasured(my.spct)
setWhereMeasured(my.spct, lat = 60, lon = -10)
getWhereMeasured(my.spct)
getWhereMeasured(my.spct)$lon
my.spct

## ------------------------------------------------------------------------
is_effective(sun.spct)
is_effective(sun.spct * waveband(c(400, 700)))

## ------------------------------------------------------------------------
ten.minutes.spct <-
  convertTimeUnit(sun.spct, time.unit = duration(10, "minutes"))
ten.minutes.spct
getTimeUnit(ten.minutes.spct)

## ------------------------------------------------------------------------
two_suns.mspct <- source_mspct(list(sun1 = sun.spct, sun2 = sun.spct))
two_suns.mspct

## ------------------------------------------------------------------------
mixed.mspct <- generic_mspct(list(filter = clear.spct, source = sun.spct))

## ------------------------------------------------------------------------
two_gen.mscpt <- as.generic_mspct(two_suns.mspct)
class(two_gen.mscpt)
lapply(two_gen.mscpt, class_spct)

## ------------------------------------------------------------------------
two_suns.spct <- rbindspct(list(a = sun.spct, b = sun.spct / 2))
subset2mspct(two_suns.spct)

## ------------------------------------------------------------------------
test1.df <- data.frame(w.length = rep(200:210, 2),
                       s.e.irrad = rep(c(1, 2), c(11, 11)),
                       spectrum = factor(rep(c("A", "B"), c(11,11))))
subset2mspct(test1.df, member.class = "source_spct", idx.var = "spectrum")

## ------------------------------------------------------------------------
test2.df <- data.frame(w.length = 200:210, A = 1, B = 2, z = "A")
split2source_mspct(test2.df)
split2source_mspct(test2.df, spct.data.var = "s.q.irrad")

## ------------------------------------------------------------------------
is.source_mspct(two_suns.mspct)
class(two_suns.mspct)

## ------------------------------------------------------------------------
is.filter_mspct(mixed.mspct)
is.any_mspct(mixed.mspct)
class(mixed.mspct)
lapply(mixed.mspct, class_spct)
lapply(mixed.mspct, class)

## ------------------------------------------------------------------------
two_suns.mspct[1]
two_suns.mspct[2:1]

## ------------------------------------------------------------------------
two_suns.mspct[1:2] <- two_suns.mspct[2:1]

## ------------------------------------------------------------------------
two_suns.mspct[[1]]
two_suns.mspct$sun1
two_suns.mspct[["sun1"]]

## ------------------------------------------------------------------------
two_suns.mspct[["sun1"]] <- sun.spct * 2
two_suns.mspct[["sun2"]] <- NULL
two_suns.mspct

## ------------------------------------------------------------------------
c(two_suns.mspct, mixed.mspct)

## ------------------------------------------------------------------------
two.mspct <- source_mspct(list(A = sun.spct * 1, B = sun.spct * 2))
msmsply(two.mspct, `+`, 0.1)

## ------------------------------------------------------------------------
msmsply(two.mspct, trim_wl, range = c(281, 500), fill = NA)

## ------------------------------------------------------------------------
msdply(two.mspct, max)

## ------------------------------------------------------------------------
ranges.df <- msdply(two.mspct, range)
ranges.df
cat(comment(ranges.df))

## ------------------------------------------------------------------------
msdply(two.mspct, range, na.rm = TRUE)

## ------------------------------------------------------------------------
str(mslply(two.mspct, names))

## ------------------------------------------------------------------------
str(msaply(two.mspct, max))

## ------------------------------------------------------------------------
str(msaply(two.mspct, range))

## ------------------------------------------------------------------------
convolve_each(two.mspct, sun.spct)

## ------------------------------------------------------------------------
convolve_each(sun.spct, two.mspct)

## ------------------------------------------------------------------------
another_two.mspct <- two.mspct
names(another_two.mspct) <- c("a", "b")
convolve_each(another_two.mspct, two.mspct)

## ------------------------------------------------------------------------
convolve_each(two.mspct, sun.spct, oper = `+`)

## ------------------------------------------------------------------------
getWhenMeasured(two.mspct)
setWhenMeasured(two.mspct, ymd("2015-10-31", tz = "EET"))
getWhenMeasured(two.mspct)
setWhenMeasured(two.mspct,
                list(ymd_hm("2015-10-31 10:00", tz = "EET"),
                     ymd_hm("2015-10-31 11:00", tz = "EET")))
getWhenMeasured(two.mspct)
two.mspct

## ------------------------------------------------------------------------
PAR <- waveband(c(400, 700), wb.name = "PAR")
UVA <- waveband(c(315, 400), wb.name = "UVA")
UVB <- waveband(c(280, 315), wb.name = "UVB")
UVC <- waveband(c(100, 280), wb.name = "UVC")
UV  <- waveband(c(100, 400), wb.name =  "UV")
UV_bands <- list(UVC, UVB, UVA)

## ------------------------------------------------------------------------
CIE_e_fun <-
function(w.length){
    CIE.energy <- numeric(length(w.length))
    CIE.energy[w.length <= 298] <- 1
    CIE.energy[(w.length > 298) & (w.length <= 328)] <-
      10^(0.094*(298-w.length[(w.length > 298) & (w.length <= 328)]))
    CIE.energy[(w.length > 328) & (w.length <= 400)] <-
      10^(0.015*(139-w.length[(w.length > 328) & (w.length <= 400)]))
    CIE.energy[w.length > 400] <- 0
    return(CIE.energy)
}

## ------------------------------------------------------------------------
CIE <- waveband(c(250, 400), weight = "SWF",
                SWF.e.fun = CIE_e_fun, SWF.norm = 298)

## ------------------------------------------------------------------------
waveband(sun.spct)

## ------------------------------------------------------------------------
is.waveband(PAR)

## ------------------------------------------------------------------------
is_effective(waveband(c(400,500)))

## ------------------------------------------------------------------------
wavebands <- list(waveband(c(300,400)), waveband(c(400,500)))
wavebands

## ------------------------------------------------------------------------
split_bands(c(200, 225, 300))
split_bands(c(200, 225, 300), length.out = 2)

## ------------------------------------------------------------------------
split_bands(sun.spct, length.out = 2)
split_bands(PAR, length.out = 2)
split_bands(c(200, 800), length.out = 3)

## ------------------------------------------------------------------------
split_bands(list(A = c(200, 300), B = c(400, 500), C = c(250, 350)))
split_bands(list(c(100, 150, 200), c(800, 825)))

## ------------------------------------------------------------------------
split_bands(UV_bands, length.out  =  2)
split_bands(list(c(100, 150, 200), c(800, 825)), length.out = 1)

## ------------------------------------------------------------------------
sun.spct * sun.spct

## ------------------------------------------------------------------------
sun.spct * polyester.spct

## ------------------------------------------------------------------------
sun.spct * polyester.spct * polyester.spct
sun.spct * polyester.spct^2

## ------------------------------------------------------------------------
sun.spct * 2
2 * sun.spct
sun.spct * c(0,1)

## ------------------------------------------------------------------------
sun.spct * UVB

## ----eval=FALSE----------------------------------------------------------
#  sun.spct * CIE

## ------------------------------------------------------------------------
-sun.spct
sqrt(sun.spct)

## ------------------------------------------------------------------------
options(photobiology.radiation.unit = "photon")
sun.spct * UVB
options(photobiology.radiation.unit = "energy")
sun.spct * UVB

## ------------------------------------------------------------------------
# STOPGAP
shade.spct <- sun.spct

## ------------------------------------------------------------------------
rbindspct(list(sun.spct, shade.spct))
rbindspct(list(A = sun.spct, B = shade.spct), idfactor = "site")

## ------------------------------------------------------------------------
sun.spct[1:10, ]
sun.spct[1:10, 1]
sun.spct[1:10, 1, drop = TRUE]
sun.spct[1:10, "w.length", drop = TRUE]

## ------------------------------------------------------------------------
subset(sun.spct, s.e.irrad > 0.2)
subset(sun.spct, w.length > 600)
subset(sun.spct, c(TRUE, rep(FALSE, 99)))

## ------------------------------------------------------------------------
e2q(sun.spct, "add")
e2q(sun.spct, "replace")

## ------------------------------------------------------------------------
normalize(sun.spct)

## ------------------------------------------------------------------------
normalize(sun.spct, range = PAR, norm = "max")

## ------------------------------------------------------------------------
normalize(sun.spct, norm = 600.3)

## ------------------------------------------------------------------------
fscale(sun.spct)
fscale(sun.spct, f = "total")
fscale(sun.spct, range = PAR, f = irrad)

## ------------------------------------------------------------------------
fshift(sun.spct, range = UVB, f = "mean")
fshift(sun.spct, range = c(280,290), f = "min")

## ------------------------------------------------------------------------
clean(sun.spct - 0.01, range = c(280.5, 282))

## ------------------------------------------------------------------------
clean(polyester.spct - 0.053)

## ------------------------------------------------------------------------
interpolate_wl(sun.spct, seq(400, 500, by = 0.1))

## ------------------------------------------------------------------------
clip_wl(sun.spct, range = c(400, 402))
clip_wl(sun.spct, range = c(400, NA))

## ------------------------------------------------------------------------
clip_wl(sun.spct, range = UVA)

## ------------------------------------------------------------------------
clip_wl(sun.spct, range = c(100, 200))

## ------------------------------------------------------------------------
trim_wl(sun.spct, c(282.5, NA))
clip_wl(sun.spct, c(282.5, NA))

## ------------------------------------------------------------------------
trim_wl(sun.spct, PAR)

## ------------------------------------------------------------------------
trim_wl(sun.spct, c(281.5, NA), fill = NA)

## ------------------------------------------------------------------------
trim_wl(sun.spct, c(275, NA), fill = 0)

## ------------------------------------------------------------------------
trim_wl(sun.spct, c(281.5, NA), fill = NA)
trim_wl(sun.spct, c(281.5, NA), fill = NA, use.hinges = FALSE)

## ------------------------------------------------------------------------
sun.spct * CIE

## ----weighted-spectra-01, tidy=FALSE-------------------------------------
weighted.s.e.irrad <-
  with(sun.spct,
       s.e.irrad * calc_multipliers(w.length, CIE)
  )

## ------------------------------------------------------------------------
tag(sun.spct, PAR, byref = FALSE)
tag(sun.spct, UV_bands, byref = FALSE)

## ------------------------------------------------------------------------
tg.sun.spct <- tag(sun.spct, PAR, byref = FALSE)
attr(tg.sun.spct, "spct.tags")

## ------------------------------------------------------------------------
wb2tagged_spct(UV_bands)
wb2rect_spct(UV_bands)

## ------------------------------------------------------------------------
tg.sun.spct
is_tagged(tg.sun.spct)
untag(tg.sun.spct)
is_tagged(tg.sun.spct)

## ------------------------------------------------------------------------
summary(sun.spct)

## ------------------------------------------------------------------------
range(sun.spct)
min(sun.spct)
max(sun.spct)
midpoint(sun.spct)
spread(sun.spct)
stepsize(sun.spct)

## ------------------------------------------------------------------------
filters.mspct <- filter_mspct(list(none = clear.spct,
                                   pet = polyester.spct,
                                   yellow = yellow_gel.spct))
range(filters.mspct)

## ------------------------------------------------------------------------
peaks(sun.spct, span = 51)
valleys(sun.spct, span = 51)

## ------------------------------------------------------------------------
peaks(sun.spct, span = 51, unit.out = "photon")

## ------------------------------------------------------------------------
peaks(sun.spct, span = 21)

## ------------------------------------------------------------------------
msmsply(filters.mspct, peaks, span = 11)

## ------------------------------------------------------------------------
irrad(sun.spct)

## ------------------------------------------------------------------------
irrad(sun.spct, PAR)

## ------------------------------------------------------------------------
e_irrad(sun.spct, PAR) # W m-2
q_irrad(sun.spct, PAR) * 1e6 # umol s-1 m-2

## ------------------------------------------------------------------------
irrad(sun.spct, PAR, time.unit = "hour")
irrad(sun.spct, PAR, time.unit = duration(8, "hours"))

## ------------------------------------------------------------------------
irrad(sun.daily.spct, PAR, time.unit = "second")

## ------------------------------------------------------------------------
e_irrad(sun.spct, UV_bands) # W m-2

## ------------------------------------------------------------------------
irrad(sun.spct, UV_bands, quantity = "total")
irrad(sun.spct, UV_bands, quantity = "contribution")
irrad(sun.spct, UV_bands, quantity = "relative")
irrad(sun.spct, UV_bands, quantity = "average")

## ------------------------------------------------------------------------
names(filters.mspct)

## ------------------------------------------------------------------------
filtered_sun <- convolve_each(filters.mspct, sun.spct)
irrad(filtered_sun, list(UVA, PAR))

## ------------------------------------------------------------------------
irrad(convolve_each(filters.mspct, sun.spct), list(UVA, PAR))

## ------------------------------------------------------------------------
filtered_sun <- msmsply(filters.mspct, `*`, sun.spct)
irrad(filtered_sun, list(UVA, PAR))

## ------------------------------------------------------------------------
with(sun.spct, photon_irradiance(w.length, s.e.irrad, PAR))

## ------------------------------------------------------------------------
fluence(sun.spct, exposure.time = duration(1, "hours"))
fluence(sun.spct, exposure.time = 3600) # seconds

## ------------------------------------------------------------------------
e_fluence(sun.spct, PAR, exposure.time = duration(25, "minutes"))

## ------------------------------------------------------------------------
q_ratio(sun.spct, UVB, PAR)
q_ratio(sun.spct,
        list(UVC, UVB, UVA,
        UV))
q_ratio(sun.spct,
        UVB,
        list(UV, PAR))

## ------------------------------------------------------------------------
qe_ratio(sun.spct, list(UVB, PAR))

## ------------------------------------------------------------------------
q_ratio(filtered_sun, list(UVB, UVA, PAR))

## ----example-ratios-01---------------------------------------------------
with(sun.data,
     photon_ratio(w.length, s.e.irrad, UVB, PAR))

## ------------------------------------------------------------------------
normalized_diff_ind(sun.spct,
                    waveband(c(400, 700)), waveband(c(700, 1100)),
                    irrad)

## ------------------------------------------------------------------------
transmittance(polyester.spct, list(UVB, UVA, PAR))

## ------------------------------------------------------------------------
irrad(sun.spct * polyester.spct, list(UVB, UVA, PAR, wb.trim = TRUE)) /
  irrad(sun.spct, list(UVB, UVA, PAR, wb.trim = TRUE))

## ------------------------------------------------------------------------
transmittance(filters.mspct, list(UVA, PAR))

## ------------------------------------------------------------------------
response(photodiode.spct)

## ------------------------------------------------------------------------
e_response(photodiode.spct, list(UVB, UVA))

## ------------------------------------------------------------------------
sensors <- response_mspct(list(GaAsP = photodiode.spct,
                               CCD = ccd.spct))
response(sensors, list(UVB, UVA, PAR), quantity = "contribution")

## ------------------------------------------------------------------------
integrate_spct(sun.spct)

## ------------------------------------------------------------------------
average_spct(sun.spct)

## ------------------------------------------------------------------------
sun_angles(now(), lat = 34, lon = 0)
sun_angles(ymd_hms("2014-01-01 0:0:0", tz = "UTC"))

## ------------------------------------------------------------------------
sun_angles(getWhenMeasured(sun.spct), geocode = getWhereMeasured(sun.spct))

## ------------------------------------------------------------------------
dates <- seq(from = ymd("2015-03-01"), to = ymd("2015-07-1"), length.out = 3)

## ------------------------------------------------------------------------
noon_time(dates, tz = "UTC", lat =  60)
noon_time(dates, tz = "CET", lat =  60)

## ------------------------------------------------------------------------
day_night(dates, lat =  60)

## ------------------------------------------------------------------------
sunrise_time(lat = 60)

## ------------------------------------------------------------------------
sunrise_time(today(tz = "UTC"), tz = "UTC", lat = 60, lon = 0)
sunrise_time(today(tz = "EET"), tz = "EET", lat = 60, lon = 25)

## ------------------------------------------------------------------------
sunrise_time(dates, lat =  60)
sunrise_time(dates, lat = -60)

## ------------------------------------------------------------------------
sunrise_time(today(tz = "EET"), tz = "EET", lat = 60, lon = 25,
             twilight = "civil")
sunrise_time(today(tz = "EET"), tz = "EET", lat = 60, lon = 25,
             twilight = -10)
sunrise_time(today(tz = "EET"), tz = "EET", lat = 60, lon = 25,
             twilight = +12)

## ------------------------------------------------------------------------
sunrise_time(today(tz = "EET"), tz = "EET", lat = 60, lon = 25,
             unit.out = "hour")

## ------------------------------------------------------------------------
day_length(dates, lat = 60)
night_length(dates, lat = 60)

## ------------------------------------------------------------------------
day_night(dates, lat = 60)
day_night(dates, lat = 60, unit.out = "hour")

## ------------------------------------------------------------------------
w_length2rgb(550) # green
w_length2rgb(630) # red
w_length2rgb(c(550, 630, 380, 750)) # vectorized

## ------------------------------------------------------------------------
w_length_range2rgb(c(400,700))

## ------------------------------------------------------------------------
with(sun.spct, s_e_irrad2rgb(w.length, s.e.irrad))
with(sun.spct, s_e_irrad2rgb(w.length, s.e.irrad, sens = ciexyzCMF2.spct))

## ------------------------------------------------------------------------
rgb_spct(sun.spct)
rgb_spct(sun.spct, sens = ciexyzCMF2.spct)

## ------------------------------------------------------------------------
color(sun.spct)
color(sun.spct * yellow_gel.spct)

