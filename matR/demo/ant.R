
####  what are people sharing in MG-RAST ?

dir.MGRAST (75, len=25)

search.MGRAST (metadata="antarctica")

####  we could have searched by functional or taxonomic annotation, as well.
####  now here is another search for "antarctica", but with more detail.

search.MGRAST (metadata="antarctica", detail=TRUE) [ , c("name","project1")]

####  here is how to retrieve and inspect some data for a project (it takes some time).

xx <- biomRequest('mgp10307', request='organism', source='Greengenes')

colnames(xx)

rownames(xx)

####  now we can begin to study it in detail.

xx_log <- transform (xx, t_Log)

columns (xx_log, "biome")

princomp (xx_log, map = c(col="biome"), label.cex=0.5)
