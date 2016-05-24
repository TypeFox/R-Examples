lung <- data.frame(
volume = c(
  3.3, 3.1, 4.0,
  2.5, 2.6, 2.8,
  3.1, 3.5, 4.1,
  3.0, 3.7, 3.5,
  2.8, 3.6, 3.9,
  2.9, 2.8, 2.9),
method  = gl(3,1,18, labels = LETTERS[1:3]),
subject = gl(6,3,18))
