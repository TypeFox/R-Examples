autolayout <- function (
		N, 
		byrow = TRUE, 
		...) {
  if (N < 1 | N > 33) stop("I am not able to determine layout for more than 33 plots automatically.")
lay <- switch(N,
  matrix(1, ncol = 1),
  matrix(1:2, ncol = 2, byrow = byrow),
  matrix(c(0, 1, 1, 0, 2, 2, 3, 3), nrow = 2, byrow = byrow), 
  matrix(1:4, nrow = 2, byrow = byrow),
  matrix(c(0, 1, 1, 2, 2, 0, 3, 3, 4, 4, 5, 5), nrow = 2, byrow = byrow),
  matrix(1:6, nrow = 2, byrow = byrow),
  matrix(c(0, 1,1,2,2,0,rep(3:5, each=2),0, 6,6,7,7,0), nrow = 3, byrow = byrow),
  matrix(1:8, nrow = 2, byrow = byrow),
  matrix(1:9, nrow = 3, byrow = byrow),
  matrix(1:10, nrow = 2, byrow = byrow),
  matrix(c(rep(1:4,each=2), 0, rep(5:7,each=2), 0, rep(8:11,each=2)), nrow = 3, byrow = byrow),
  matrix(1:12, nrow = 3, byrow = byrow),
  matrix(c(0, rep(c(1:3,0,4:6,0,7:9),each=2),0,rep(10:13,each=2)), nrow = 4, byrow = byrow),
  matrix(c(rep(1:4,each=2), 0, rep(c(5:7,0,8:10),each=2), 0, rep(11:14, each=2)), nrow = 4, byrow = byrow),
  matrix(c(0,1,1,2,2,3,3,0,rep(4:15,each=2)), nrow = 4, byrow = byrow),
  matrix(1:16, nrow = 4, byrow = byrow),
  matrix(c(0, rep(c(1:4,0,5:8,0,9:12),each=2),0, rep(13:17,each=2)), nrow = 4, byrow = byrow),
  matrix(c(rep(1:5,each=2), 0, rep(c(6:9,0,10:13),each=2), 0, rep(14:18,each=2)), nrow = 4, byrow = byrow),
  matrix(c(0, 1,1,2,2,3,3,4,4,0,rep(5:19,each=2)), nrow = 4, byrow = byrow),
  matrix(1:20, nrow = 4, byrow = byrow),
  matrix(c(rep(0:4,each=2), 0, 0, 0, rep(5:9,each=2), 0, rep(10:21,each=2)), nrow = 4, byrow = byrow),
  matrix(c(rep(1:6,each=2), 0, rep(7:11,each=2),0, rep(12:17,each=2),0,rep(18:22,each=2),0), nrow = 4, byrow = byrow),
  matrix(c(0, rep(1:5,each=2), 0, rep(6:23,each=2)), nrow = 4, byrow = byrow),
  matrix(1:24, nrow = 4, byrow = byrow),
  matrix(c(rep(0:5,each=2), 0, 0, 0, rep(6:11,each=2), 0, rep(12:25, each=2)), nrow = 4, byrow = byrow),
  matrix(c(0, rep(1:6,each=2), 0, rep(7:20,each=2), 0, rep(21:26,each=2), 0), nrow = 4, byrow = byrow),
  matrix(c(0, rep(1:6,each=2), 0, rep(7:27,each=2)), nrow = 4, byrow = byrow),
  matrix(1:28, nrow = 4, byrow = byrow),
  matrix(c(rep(0:6,each=2), 0, 0, 0, rep(7:13,each=2), 0, rep(14:29, each=2)), nrow = 4, byrow = byrow),
  matrix(1:30, nrow = 5, byrow = byrow),
  matrix(c(0,rep(1:7,each=2), 0, rep(8:31, each=2)), nrow = 4, byrow = byrow),
  matrix(1:32, nrow = 4, byrow = byrow),
  matrix(1:33,ncol=11, byrow = byrow)
  )
  invisible(layout(lay, ...))
}

