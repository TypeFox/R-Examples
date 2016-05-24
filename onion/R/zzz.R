"O0"  <- octonion(Re=0)
"O1"  <- octonion(Re=1)
"Oi"  <- octonion( i=1)
"Oj"  <- octonion( j=1)
"Ok"  <- octonion( k=1)
"Ol"  <- octonion( l=1)
"Oil" <- octonion(il=1)
"Ojl" <- octonion(jl=1)
"Okl" <- octonion(kl=1)
"Oall" <- as.octonion(rep(1,8),single=TRUE)
"Oim" <- Oall-O1

"H0"  <- quaternion(Re=0)
"H1"  <- quaternion(Re=1)
"Hi"  <- quaternion( i=1)
"Hj"  <- quaternion( j=1)
"Hk"  <- quaternion( k=1)
"Hall" <- as.quaternion(rep(1,4),single=TRUE)
"Him" <- Hall-H1

