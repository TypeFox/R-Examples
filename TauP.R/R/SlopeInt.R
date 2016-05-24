SlopeInt <-
function(v,z){
g=(v[1]-v[2])/(z[1]-z[2])
v0=v[2]-z[2]*g
return(list(g=g,v0=v0))
}

