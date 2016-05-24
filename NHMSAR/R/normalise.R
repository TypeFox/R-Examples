normalise <-
function(M){

c = sum(M)
# Set any zeros to one before dividing
d = c + (c==0);
M = M / d;
return(M)
}
