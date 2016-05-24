#
# a bunch of mixed edits 

# define datamodels
A %in% letters[1:6]
B %in% letters[7:10]
C %in% letters[11:20]

# some purely categorical stuff
if ( A %in% c('a','d') ) B %in% c('g','h')
if ( A == 'b' ) C %in% letters[12:18]

# some mixed numeric...
if ( C == 'k' ) x > 0
if ( y > 0  )  x + y >= z

# some pure numeric
u + v == w
2*v + 3*s == t







