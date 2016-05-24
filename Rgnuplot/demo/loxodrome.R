# loxodrome tests
Gprun("set dummy u,v\nset angles degrees\nset parametric\nset urange [ -90.0000*15 : 90.0000*15 ] noreverse nowriteback\nset vrange [ 0.00000 : 360.000 ] noreverse nowriteback\nset view equal xyz\nset nokey\nmyAngle=120.0\ncosc(x)=cos(atan(x/myAngle))\nsinc(x)=sin(atan(x/myAngle))\nsplot cos(u)*cos(v),cos(u)*sin(v),sin(u) notitle w l lt 5, cos(u)*cosc(u),sin(u)*cosc(u),-sinc(u) w l lc rgb \"red\" notit\n\n", 
    TRUE)
 
