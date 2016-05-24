nose.rr <-
function(nos,cc)
{
tcol=ncol(nos)
trow=nrow(nos)
nature=matrix(0,trow,1,byrow=FALSE)
e=cc
f=1/e
for(m in 1:trow)
{
a=nos[m,1]
b=nos[m,2]
c=nos[m,3]
d=nos[m,4]
if(a>=0 && b>=0 && c>=0 && d>=0 && e>=0)
{
if(a>0 && b>0 && c>0 && d>0)
{nature[m]='atleast one of the four arguments must be zero'}
n=a+b
p=c+d
#case 1
if(a==0 && b>0 && c>0 && d>0)
{
if(d>(b+c))
{cl1= (c * n) / (d - c - b)
clr1=round(cl1,3)
if(e>clr1)
{nature[m]='SEVERE SPARSE'}
else
{nature[m]='MILD SPARSE'}
}
else
{nature[m]='MILD SPARSE'}
}
#case 2
if(c==0 && a>0 && b>0 && d>0)
{
if(b<=(d+a))
{nature[m]='MILD SPARSE'}
else
{cl2 = (a * p) / (b - a - d)
clr2=round(cl2,3)
if(e>clr2)
{nature[m]='SEVERE SPARSE'}
else
{nature[m]='MILD SPARSE'}
}

}
#CASE 3
if(b==0 && a>0 && c>0 && d>0)
{nature[m]='MILD SPARSE'}
#CASE 4
if(d==0 && a>0 && b>0 && c>0)
{nature[m]='MILD SPARSE'}

#CASE 5
if(a==0 && b==0 && c>0 && d>0)
{nature[m]='MILD SPARSE'}

#CASE 6
if(a==0 && c==0 && b>0 && d>0)
{nature[m]='MILD SPARSE'}
#CASE 7
if(a==0 && d==0 && b>0 && c>0)
{nature[m]='MILD SPARSE'}
#CASE 8
if(b==0 && c==0 && a>0 && d>0)
{nature[m]='MILD SPARSE'}
#CASE 9
if(b==0 && d==0 && a>0 && c>0)
{
if(a==c)
{nature[m]='MILD SPARSE'}
else
{nature[m]='SEVERE SPARSE'}
}
#CASE 10
if(c==0 && d==0 && a>0 && b>0)
{nature[m]='MILD SPARSE'}
#CASE 11
if(a==0 && b==0 && c==0 && d>0)
{nature[m]='MILD SPARSE'}
#CASE 12
if(a==0 && b==0 && d==0 && c>0)
{nature[m]='MILD SPARSE'}
#CASE 13
if(a==0 && c==0 && d==0 && b>0)
{nature[m]='MILD SPARSE'}
#CASE 14
if(b==0 && c==0 && d==0 && a>0)
{nature[m]='MILD SPARSE'}
#CASE 15
if(a==0 &&b==0 && c==0 && d==0 )
{nature[m]='MILD SPARSE'}
}
else
{
nature[m]='arguments must be non negative numbers'
}
}
nosrr=cbind(nos,nature)
nosrr
}
