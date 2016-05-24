bin_ps <-
function(x,n,p)
{

random_needed    = 0
iteration_needed = numeric(x)
iter             = 0
Result           = numeric(x)
element          = 0


while (element < x)
{
coalescence = 0
T           = 1
V_up_down   = numeric(500)
Random_aux  = numeric(500)
l_value     = 0
u_value     = n
iterate     = 1
iter        = 0

while (iterate==1)
{
iter = iter +1

for (i in -T:-1)
{
if (i==-T)
{
l_value=0
u_value=n
}
else
{
if (i <= -T/2)
{
V              = draw_i()
V_up_down[-i]  = V
U              = runif(1)
Random_aux[-i] = U
random_needed  = random_needed+2
}
else (i > -T/2)
{
V              = V_up_down[-i]
U              = Random_aux[-i]
}

true_e1 = min(1,dbinom(l_value+V,n,p)/ dbinom(l_value,n,p))
true_e2 = min(1,dbinom(u_value+V,n,p)/ dbinom(u_value,n,p))

if (U < true_e1)
{
l_value = l_value + V
}
if (U < true_e2)
{
u_value  = u_value + V
}

if (l_value==u_value)
{
coalescence =1
}

}
}

if (coalescence == 1)
{
element                    = element+1
iterate                    = 0
Result[element]            = l_value
iteration_needed[element]  = iter
} 

if (coalescence == 0)
{
T=2*T
}
}
}
ret=list(values=Result,iters_until_coalescence=iteration_needed,rand_used=random_needed)
return(ret)
}

