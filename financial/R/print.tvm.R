"print.tvm" <-
function (x,...) 
{

cat("\nTime Value of Money model\n\n");

class(x)=NULL;

print(round(x,2));
cat("\n");

}

