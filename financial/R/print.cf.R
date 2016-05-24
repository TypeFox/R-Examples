"print.cf" <-
function (x,...) 
{

cat("\nCash Flow Model\n\n");

cat("Flows:\n");
print(x$cf);
cat("\n IRR%:",round(x$irr,2),"\n");
cat(" NPV Extremes at I%:",round(x$ext,2),"\n");
if (!is.null(x$mirr)) { cat("\n"); print(round(x$mirr,2)); }
cat("\n");

if (!is.null(x$tab)) { print(round(x$tab,2)); }

}

