blisters <- data.frame(
	time = c( 5,6,6,7,7,8,9,10,
              7,7,8,9,9,10,10,11,
              7,9,9,10,10,10,11,12,13),
	treatment = c(rep('A', 8), rep('B', 8), rep('P', 9))
	)
small <- data.frame(
	response = c(
	5.3, 6.0, 6.7, 
	5.5, 6.2, 6.4, 5.7,
	7.5, 7.2, 7.9),
	group = toupper(letters[c(1,1,1,2,2,2,2,3,3,3)])
	)

