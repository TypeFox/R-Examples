nose.or <-
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
if(a>=0 && b>=0 && c>=0 && d>=0)
{
if(a>0 && b>0 && c>0 && d>0)
{nature[m]='atleast one of the four arguments must be zero'
}

#case 1
if(a==0 && b>0 && c>0 && d>0 && e>0)
{    if (d <= (b + c))
   {nature[m]='MILD SPARSE'}
else
{if (d <= (b + c) + (f * b * c))
{nature[m]='MODERATE SPARSE'}
else
{nature[m]='SEVERE  SPARSE'}}
}
#case 2
if(c==0 && a>0 && b>0 && d>0)
{if (b <= (a + d))
    {nature[m]='MILD SPARSE'}
else
{if (b <= (a + d) + (f * a * d))
        {nature[m]='MODERATE SPARSE'}
else
{nature[m]='SEVERE  SPARSE'}}
}

#CASE 3
if(b==0 && a>0 && c>0 && d>0)
{if (c <= (a + d))
    {nature[m]='MILD SPARSE'}
else
{if (c <= (a + d) + (f * a * d))
        {nature[m]='MODERATE SPARSE'}
else
{nature[m]='SEVERE  SPARSE'}}
}
#CASE 4
if(d==0 && a>0 && b>0 && c>0)
{if (a <= (b + c))
    {nature[m]='MILD SPARSE'}
else
{if (a <= (b + c) + (f * a * d))
        {nature[m]='MODERATE SPARSE'}
else
{nature[m]='SEVERE  SPARSE'}}
}
#CASE 5
if(a==0 && b==0 && c>0 && d>0)
{if (d <= 2 * c * (c + 1))
{nature[m]='MILD SPARSE'}
else
    {nature[m]='SEVERE  SPARSE'}}
#CASE 6
if(a==0 && c==0 && b>0 && d>0)
{if (d <= 2 * b * (b + 1))
{nature[m]='MILD SPARSE'}
else
    {nature[m]='SEVERE  SPARSE'}}
#CASE 7
if(a==0 && d==0 && b>0 && c>0)
{nature[m]='MILD SPARSE'}
#CASE 8
if(b==0 && c==0 && a>0 && d>0)
{nature[m]='MILD SPARSE'}
#CASE 9
if(b==0 && d==0 && a>0 && c>0)
{if (a <= 2 * c * (c + 1))
{nature[m]='MILD SPARSE'}
else
    {nature[m]='SEVERE  SPARSE'}}

#CASE 10
if(c==0 && d==0 && a>0 && b>0)
{if (a <= 2 * b * (b + 1))
{nature[m]='MILD SPARSE'}
else
    {nature[m]='SEVERE  SPARSE'}}

#CASE 11
if(a==0 && b==0 && c==0 && d>0)
{nature[m]='SEVERE SPARSE'}
#CASE 12
if(a==0 && b==0 && d==0 && c>0)
{nature[m]='MODERATE SPARSE'}
#CASE 13
if(a==0 && c==0 && d==0 && b>0)
{nature[m]='SEVERE SPARSE'}
#CASE 14
if(b==0 && c==0 && d==0 && a>0)
{nature[m]='MODERATE SPARSE'}
#CASE 15
if(a==0 &&b==0 && c==0 && d==0 )
{nature[m]='MILD SPARSE'}
}
else
{
nature[m]='arguments must be non negative numbers'
}
}
nosor=cbind(nos,nature)
nosor
}
