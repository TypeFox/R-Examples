GetDetrendMethod = function(method,  n, nPerc, p){
Method= "No detrending"
if (method=="Mean") {Method = "Mean"}

if (method=="ModNegExp") {Method = "Exponential"}

if (method=="Spline") {Method = paste("Spline n=", n, ", p=", p,sep="")}

if (method=="Spline%") {Method = paste("Spline n%=", nPerc, ", p=", p, sep="")}
	
if (method=="No Detrending") {Method = paste("No detrending", sep="")}

return(Method)
}
