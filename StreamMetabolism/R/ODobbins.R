`ODobbins` <-
function(vel, dep){
	#vel in m/sec and dep in m
	#changed 2013-11-12 from 3.74 to 3.93; this converts units to 1/d
	#vel in m/s and depth in m; k in 1/d
	(3.93*(vel^0.5))/(dep^1.5)
}

