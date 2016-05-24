gmat<-function(vec, p, alpha)
{
####   get a rotation matrix about an arbitrary vector
####     uses translation also
####  vec = 3D direction vector
####    p = point to translate to
####    angle of rotation about vec in degrees
  v=vec/sqrt(sum(vec^2))

  r1=roty4(v)
  r2=rotx4(v)
  r3=rotdelta4(alpha)

  t1=trans4(c(-p[1], -p[2], -p[3]))

  M<- t1 %*% r2 %*% r1
#### no translation    M<- r2 %*% r1
  MI<- solve(M)
  ans=M %*% r3 %*% MI
  return(ans)
}
