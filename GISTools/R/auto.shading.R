`auto.shading` <-
function(x,digits=2,cutter=quantileCuts,n=5,params=NA,cols=brewer.pal(n,'Reds')) {
	brk = cutter(x,n=n,params=params)
	if (!is.na(digits)) brk = signif(brk,digits=digits)
	brk = sort(brk)
	brk = brk[!duplicated(brk)]
	res=list(breaks=brk,cols=cols)
	class(res) = 'shading' 
	res}

poly.outer <- function (exo.object, input.poly, extend = 0) 
{
    box <- function(obj, extend) {
        xy = bbox(obj)
        x = xy[1, ] + c(-extend, extend)
        y = xy[2, ] + c(-extend, extend)
         SpatialPolygons(list(Polygons(list(Polygon(cbind(x[c(1, 1, 2, 2,1)], y[c(1, 2, 2, 1,1)]))), ID=1)))
    }
    choppo <- function(obj1, obj2, extend) {
        res = box(obj1, extend)
        res = gDifference(res, gUnaryUnion(obj2))
    }
    res = choppo(exo.object, input.poly, extend)
    proj4string(res) = CRS(proj4string(input.poly))
res
}

add.masking <- function (maskPoly, color = "white") 
  plot(maskPoly, col = color, border = color, usePolypath=TRUE, add = TRUE)
