"msep" <-
function(n,s,p.tr)
{

expected=0
for(y in 0:n)
  {
  expected=expected+( (1-(1-y/n)^(1/s)) * choose(n,y) * ((1-(1-p.tr)^s)^y) * ((1-p.tr)^(s*(n-y)) ) )
  }
expected


varsum=0
for(y in 0:n)
  {
  varsum=varsum+ ( ( (1-y/n)^(2/s) ) * choose(n,n-y) * ( ((1-p.tr)^s)^(n-y) ) * ((1-(1-p.tr)^s)^y) ) 
  }
varp=varsum - (1-expected)^2

expp=expected
bias=expected-p.tr
mse=varp + bias^2


list(varp=varp,
  mse=mse,
  bias=bias,
exp=expected)
}

