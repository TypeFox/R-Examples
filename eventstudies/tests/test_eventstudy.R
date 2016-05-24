library(eventstudies)

# An example dataset, with 3 firms --
p <- structure(c(33.16, 34.0967, 35.3683, 34.46, 34.17, 35.89, 36.19,
                 37.1317, 36.7033, 37.7933, 37.8533, 285.325, 292.6,
                 290.025, 286.2, 290.075, 295.05, 289.325, 285.625,
                 293.7, 298.5, 289.05, 704.5438, 708.35, 735.8375,
                 710.625, 711.65, 731.0125, 727.575, 715.0187, 724.2,
                 713.1875, 695.1812), .Dim = c(11L, 3L), .Dimnames =
                 list( NULL, c("ITC", "Reliance", "Infosys")), index =
                 structure(c(12418, 12419, 12422, 12423, 12424, 12425,
                 12426, 12429, 12430, 12431, 12432), class = "Date"),
                 class = "zoo")
# An example events list
eventslist <- data.frame(unit=c("ITC","Reliance","Infosys",
                           "ITC","Reliance","Junk"),
                         when=as.Date(c(
                           "2004-01-02", "2004-01-08", "2004-01-14",
                           "2005-01-15", "2004-01-01", "2005-01-01")))
eventslist$unit <- as.character(eventslist$unit)

# What we expect if we don't worry about width --
rawres <- structure(list(z.e = structure(c(NA, NA, NA, NA, NA, NA,
  NA, NA, 33.16, 34.0967, 35.3683, 34.46, 34.17, 35.89, 36.19,
  37.1317, 36.7033, 37.7933, 37.8533, NA, NA, NA, NA, 285.325, 292.6,
  290.025, 286.2, 290.075, 295.05, 289.325, 285.625, 293.7, 298.5,
  289.05, NA, NA, NA, NA, 704.5438, 708.35, 735.8375, 710.625, 711.65,
  731.0125, 727.575, 715.0187, 724.2, 713.1875, 695.1812, NA, NA, NA,
  NA, NA, NA, NA, NA), .Dim = c(19L, 3L), .Dimnames = list( NULL,
  c("1", "2", "3")), index = -9:9, class = "zoo"), outcomes =
  structure(c(1L, 1L, 1L, 3L, 3L, 2L), .Label = c("success",
  "unitmissing", "wrongspan" ), class = "factor")), .Names = c("z.e",
  "outcomes"))

# Check without the width handling --
a <- phys2eventtime(p, eventslist,width=0)
all.equal(a, rawres)
# Check with width of 1 --
a <- phys2eventtime(p, eventslist,width=1)
all.equal(a, rawres)

# But when we go to width=2, column 1 and 3 drop off because they have
# only 1 obs before & after the event date respectively.
a <- phys2eventtime(p, eventslist,width=2)
all.equal(a, structure(list(z.e = structure(c(NA, NA, NA, NA, 285.325,
                              292.6, 290.025, 286.2, 290.075, 295.05,
                              289.325, 285.625, 293.7, 298.5, 289.05,
                              NA, NA, NA, NA), index = -9:9, class =
                              "zoo"), outcomes = structure(c(3L, 1L,
                              3L, 4L, 4L, 2L), .Label = c("success",
                              "unitmissing", "wdatamissing",
                              "wrongspan"), class = "factor")), .Names
                              = c("z.e", "outcomes" )))
