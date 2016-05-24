 lbinorm=function (xy, par) 
{
    m = par$m
    v = par$v
    x = xy[1]
    y = xy[2]
    zx = (x - m[1])/sqrt(v[1, 1])
    zy = (y - m[2])/sqrt(v[2, 2])
    r = v[1, 2]/sqrt(v[1, 1] * v[2, 2])
    return(-0.5/(1 - r^2) * (zx^2 - 2 * r * zx * zy + zy^2))
}
