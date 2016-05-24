


test.rand_connComp3D <- function() {
  for (i in 1:10) {
    mat = array(sample(c(0, 1), size = 10^3, replace = TRUE), dim = rep(10, 3))
    dim(mat)
    # cc = connComp3D(mat)
    mask = mat > 0
    cc = connComp3D(mask)
  }
  
  TRUE
  
}

test.twoclust_connComp3D <- function() {
  arr <- array(0,rep(50,3))
  spc <- BrainSpace(rep(50,3), spacing=c(1,1,1))
  vol <- BrainVolume(arr, spc)
  
  sphere1 <- RegionSphere(vol, c(25,25,25), 2, fill=1)
  vol[coords(sphere1)] <- 1
  
  sphere2 <- RegionSphere(vol, c(8,8,8), 1, fill=1)
  vol[coords(sphere2)] <- 1
  
  cc <- connComp3D(as.logical(vol))
  checkEquals(sum(cc$index==1), nrow(coords(sphere1)))
  checkEquals(sum(cc$index == 2), nrow(coords(sphere2)))
  checkEquals(max(cc$index), 2)
  checkEquals(max(cc$size), max(c(nrow(coords(sphere1)), nrow(coords(sphere2)))))
  
}
