adjust <-
function(age){
u=-1
v=1
n=length(age)
amin=min(age)
amax=max(age)

temp=u+(v-u)/(amax-amin)*(age-amin)
return(temp)
}

