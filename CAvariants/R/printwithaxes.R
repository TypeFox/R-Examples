printwithaxes <-
function(x, thenames,digits=3) { 
names(x) <- thenames
print(round(x,digits=digits))
}
