listmake <- function (y=1:10,g=as.factor(c(1,1,1,2,2,2,2,2,2,3)),b=as.factor(c(1,1,1,1,3,2,2,2,3,1)))
{
data.sets <- lapply(levels(b), 
					function(blv){
						gg <- as.factor(as.numeric(unlist(g[b==blv])))
						lapply(levels(gg),
							function(glv){
						     return(y[b == blv & g == glv])
							}
						)
					}
				)
data.sets
}
