"mcp" <- function(xy, id, percent=95)
{
    ## Verifications
    if (length(id)!=nrow(xy)) stop("xy and id should be of the same length")
    if (percent>100) {
	warning("The MCP is estimated using all relocations (percent>100)")
	percent<-100
    }
    if (min(table(id))<5)
        stop("At least 5 relocations are required to fit an home range")

    ## First remove the missing values
    id<-id[!is.na(xy[,1])]
    id<-id[!is.na(xy[,2])]
    xy<-xy[!is.na(xy[,1]),]
    xy<-xy[!is.na(xy[,2]),]
    id<-factor(id)

    ## Computes the centroid of the relocations for each animal
    r<-split(xy, id)
    est.cdg<-function(xy) apply(xy, 2, mean)
    cdg<-lapply(r,est.cdg)
    levid<-levels(id)

    ## Preparation of outputs
    X<-0
    Y<-0
    ID<-"0"

    ## Then, for each animal...
    for (i in 1:nlevels(id)) {

	k<-levid[i]
	df.t<-r[[levid[i]]]
	cdg.t<-cdg[[levid[i]]]

        ## Distances from the relocations to the centroid: we keep
        ## the "percent" closest
	dist.cdg<-function(xyt) {
            d<-sqrt( ( (xyt[1]-cdg.t[1])^2 ) + ( (xyt[2]-cdg.t[2])^2 ) )
            return(d)
        }

	di<-apply(df.t, 1, dist.cdg)
	key<-c(1:length(di))

	acons<-key[di<=quantile(di,percent/100)]
	xy.t<-df.t[acons,]


        ## Coordinates of the MCP
	coords.t<-chull(xy.t[,1], xy.t[,2])
	xy.bord<-xy.t[coords.t,]

	X<-c(X,xy.bord[,1])
	Y<-c(Y,xy.bord[,2])
	ID<-c(ID, rep(as.character(levid[i]), nrow(xy.bord)))
    }

    ## Outputs: an object of class "area"
    ID<-as.data.frame(ID)
    res<-cbind.data.frame(ID,X,Y)
    res<-res[-1,]
    res[,1]<-factor(res[,1])
    names(res) <- c("ID","X","Y")
    res<-as.area(res)
    return(res)
}

