perickson <-
function(t.bar=NULL, s, f, d){
if(is.null(t.bar)) t.bar <- 1/(-log(s)) 
(t.bar*f/d*((exp(d/t.bar)-1)/(exp(d/t.bar)-1+f))) 
}

# changes:
# 3. 9. 2013: argument "t_quer" was renamed into "t.bar" and possibility to use s instead of tbar added
