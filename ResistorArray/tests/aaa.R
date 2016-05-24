library(ResistorArray)

#first some Platonic results:
stopifnot(all.equal(resistance(cube(),1,7),5/6))
stopifnot(all.equal(resistance(cube(),1,2),7/12))

stopifnot(all.equal(resistance(octahedron(),1,6),1/2))
stopifnot(all.equal(resistance(octahedron(),1,5),5/12))
stopifnot(all.equal(resistance(dodecahedron(),1,5),19/30))

stopifnot(all.equal(Wu(cube())[1,2], 7/12))

#now reproduce opposite corners of skeleton cube:
a <- circuit(cube(),c(0,rep(NA,5),1,NA))
stopifnot(all.equal((a$potentials[7] - a$potentials[1])/a$currents[7] , 5/6))

# series:
stopifnot(all.equal(resistance(series(rep(1,10)),1,11),10))

# Jacob's ladder:
phi <- (sqrt(5)-1)/2
stopifnot(all.equal(resistance(ladder(30),1,2), phi))
stopifnot(all.equal( Wu(ladder(30))[1,2],phi))
