getKnotType <- function( polynomial, invariant = 'HOMFLY', full.output = FALSE ) {
	invariant <- toupper( invariant )
	implemented <- c('ALEXANDER', 'JONES', 'HOMFLY')
	if( missing(invariant) || !( invariant %in% implemented) )
		return('invariant not specified')

	polyH <- c('1', '2*l^2 + l^2*m^2 - l^4', '-1/l^4 + 2/l^2 + m^2/l^2', '-1 + l^(-2) + l^2 - m^2',
			'l^2 + l^2*m^2 + l^4*m^2 + l^4 - l^6', 'l^(-4) - 1/l^2 + l^2 - m^2 - m^2/l^2')
	polyA <- c('1', '-1 + t1 + 1/t1', '-1 + t1 + 1/t1', '3 - t1 - 1/t1',
			'-3 + 2*t1 + 2/t1', '5 - 2*t1 - 2/t1')
	polyJ <- c('1', '1/t - 1/t^4 + t^(-3)', 't + t^3 - t^4', '1 - t - 1/t + t^(-2) + t^2',
			'1/t - 1/t^6 + t^(-5) - 1/t^4 + 2/t^3 - 1/t^2', '2 - 2*t - 1/t + t^(-2) + t^2 - t^3 + t^4')
	subject <- c('Unknot',
			'Left-handed Trefoil knot (3_1*)',
			'Right-handed Trefoil knot (3_1)',
			'Figure-eight knot (4_1)',
			'3-twist knot (5_2)',
			'Stevedore\'s knot (6_1)' )
	toDBN <- c('http://katlas.org/wiki/0_1', 'http://katlas.org/wiki/3_1',
			'http://katlas.org/wiki/3_1,', 'http://katlas.org/wiki/4_1',
			'http://katlas.org/wiki/5_2', 'http://katlas.org/wiki/6_1')
	knot.table <- cbind(Knot.type = subject, HOMFLY = polyH, Jones = polyJ, Alexander = polyA, link.Knot.Atlas = toDBN)
	
	switch(invariant,
			'ALEXANDER' = {
				pos <- which( polyA == polynomial)
				if( !full.output )	return( subject[pos] )
				else return(knot.table[pos, ])
			},
			'JONES' = {
				pos <- which( polyJ == polynomial)
				if( !full.output )	return( subject[pos] )
				else return(knot.table[pos, ])
			},
			'HOMFLY' = {
				pos <- which( polyH == polynomial)
				if( !full.output )	return( subject[pos] )
				else return(knot.table[pos, ])
			})
}
			