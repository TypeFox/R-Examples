library(randtoolbox)

try( randtoolbox:::.getrandtoolboxEnv(".torus.seed") )
!exists(".torus.seed", envir=randtoolbox:::.randtoolboxEnv, mode="list")
try( torus(1, init=FALSE) )
try( randtoolbox:::.getrandtoolboxEnv(".torus.seed") )


torus(1, init=TRUE)
randtoolbox:::.getrandtoolboxEnv(".torus.seed")
!exists(".torus.seed", envir=randtoolbox:::.randtoolboxEnv, mode="list")

torus(1, init=FALSE)
randtoolbox:::.getrandtoolboxEnv(".torus.seed")
!exists(".torus.seed", envir=randtoolbox:::.randtoolboxEnv, mode="list")

torus(1, init=FALSE)
randtoolbox:::.getrandtoolboxEnv(".torus.seed")
!exists(".torus.seed", envir=randtoolbox:::.randtoolboxEnv, mode="list")




try( randtoolbox:::.getrandtoolboxEnv(".halton.seed") ) #not found
!exists(".halton.seed", envir=randtoolbox:::.randtoolboxEnv, mode="list")
try( halton(1, init=FALSE) )

try( randtoolbox:::.getrandtoolboxEnv(".halton.seed") ) 
halton(7, init=TRUE)
randtoolbox:::.getrandtoolboxEnv(".halton.seed") #now initialized
!exists(".halton.seed", envir=randtoolbox:::.randtoolboxEnv, mode="list")

halton(1, init=FALSE)
randtoolbox:::.getrandtoolboxEnv(".halton.seed") #now initialized
halton(1, init=TRUE)
randtoolbox:::.getrandtoolboxEnv(".halton.seed") #reset
halton(6, init=FALSE)
randtoolbox:::.getrandtoolboxEnv(".halton.seed") 
halton(1, init=FALSE)

