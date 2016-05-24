
 library(aster)

 aster:::setfam(fam.default())

 .C("aster_check_model",
     nind = as.integer(42),
     nnode = as.integer(4),
     pred = as.integer(c(0, 1, 2, 3)),
     fam = as.integer(c(1, 1, 2, 3)))

 options(error=dump.frames)

 .C("aster_check_model",
     nind = as.integer(0),
     nnode = as.integer(4),
     pred = as.integer(c(0, 1, 2, 3)),
     fam = as.integer(c(1, 1, 2, 3)))

 .C("aster_check_model",
     nind = as.integer(42),
     nnode = as.integer(-3),
     pred = as.integer(c(0, 1, 2, 3)),
     fam = as.integer(c(1, 1, 2, 3)))

 .C("aster_check_model",
     nind = as.integer(42),
     nnode = as.integer(4),
     pred = as.integer(c(0, 2, 2, 3)),
     fam = as.integer(c(1, 1, 2, 3)))

 .C("aster_check_model",
     nind = as.integer(42),
     nnode = as.integer(4),
     pred = as.integer(c(0, 1, 2, 3)),
     fam = as.integer(c(1, 1, 2, 25)))

