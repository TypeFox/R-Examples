`difflon` <-
function(LON1, LON2)
{

deg2rad = 0.01745329252
c1 = cos(LON1*deg2rad)
s1 = sin(LON1*deg2rad)
c2 = cos(LON2*deg2rad)
s2 = sin(LON2*deg2rad)

radd = acos(c1*c2+s1*s2)

sn = sign(c1*s2-s1*c2)

return(list(deg=radd/deg2rad, sn=sn))

}

