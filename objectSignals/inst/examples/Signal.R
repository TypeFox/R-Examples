Signal(x, y)
signal <- Signal(x, y, z = NA)
signal$connect(function(n, x, option = "none") message("x:", x),
               namedArgs = TRUE)
signal$connect(function(z, ...) message("z:", z, " x:", list(...)$x),
               namedArgs = TRUE)
signal$emit(0, 1)
##'
signal$connect(function(x, y, option = "none")
               message("y:", y, " op:", option), TRUE)
signal$connect(function(x, y, option = "none")
               message("op:", option), option = "test")
signal$connect(function(x, y, option = "none")
               message("op:", option), FALSE, "test")
id <- signal$connect(function(x, y, option = "none")
                     message("op:", option), TRUE, "test")
##'
signal$emit(0, 1)
##'
signal$disconnect(id)
signal$emit(0, 2)
##'
signal <- Signal(x)
signal$connect(function(i) print(i))
##'
signal$block()
signal$emit(0)
signal$unblock()
signal$emit(0)
##'
signal$buffer()
signal$emit(0); signal$emit(1); signal$emit(3)
signal$flush()
##'
signal$accumulator(function(prev, cur) {
  prev$x <- c(prev$x, cur$x)
  prev
})
signal$buffer()
signal$emit(0); signal$emit(1); signal$emit(3)
signal$flush()
