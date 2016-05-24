xy.grid <-
function(rangex,rangey,nx,ny=nx)
	expand.grid(seq(rangex[1],rangex[2],length=nx ),seq(rangey[1],rangey[2],length=ny ))
