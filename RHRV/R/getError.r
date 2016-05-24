getError=function(f,interval,type,delta,relative){


if (type=="lower")
{
   if (relative){
   error=abs((f-interval[1])/delta)*100
   }else{
     error=abs((f-interval[1]))
   }
}
else
{  if (relative){
    error=abs((f-interval[2])/delta)*100
    }else{
     error=abs((f-interval[2]))
   }
}

return(error)

}