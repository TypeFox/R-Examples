########## R function: polygonGrow ##########

# For growing a convex polygon by rolling
# a disc of radius "radius" (in standard
# deviation units) around the edge.

# Last changed: 03 APR 2009

polygonGrow <- function(sttPoly,r)
{
   sx <- sd(sttPoly$x) ; sy <- sd(sttPoly$y)
   sSttPoly <- list(x=sttPoly$x/sx,y=sttPoly$y/sy)

   nVert <- length(sttPoly$x)
   area <- areaPolygon(sSttPoly)
   centroid <- (c(sum((sSttPoly$x[-nVert]+sSttPoly$x[-1])*
	       	(sSttPoly$x[-1]*sSttPoly$y[-nVert]-sSttPoly$x[-nVert]*sSttPoly$y[-1])),
                 sum((sSttPoly$y[-nVert]+sSttPoly$y[-1])*
		(sSttPoly$x[-1]*sSttPoly$y[-nVert]-sSttPoly$x[-nVert]*sSttPoly$y[-1])))
                /(6*area))

   thetaVec <- atan(diff(sSttPoly$y)/diff(sSttPoly$x))
   trVec2x <- cos(thetaVec)*diff(sSttPoly$x) + sin(thetaVec)*diff(sSttPoly$y)
   ansPosx <- 0.5*cos(thetaVec)*trVec2x-2*sin(thetaVec)*r + sSttPoly$x[-nVert]
   ansPosy <- 0.5*sin(thetaVec)*trVec2x+2*cos(thetaVec)*r + sSttPoly$y[-nVert]
   ansNegx <- 0.5*cos(thetaVec)*trVec2x+2*sin(thetaVec)*r + sSttPoly$x[-nVert]
   ansNegy <- 0.5*sin(thetaVec)*trVec2x-2*cos(thetaVec)*r + sSttPoly$y[-nVert]
   distPos <- sqrt((ansPosx-centroid[1])^2+(ansPosy-centroid[2])^2)
   distNeg <- sqrt((ansNegx-centroid[1])^2+(ansNegy-centroid[2])^2)
   posWins <- as.numeric(distPos>distNeg)
   sEndPoly <- list(x=(posWins*ansPosx+(1-posWins)*ansNegx),
                 y=(posWins*ansPosy+(1-posWins)*ansNegy))
   sEndPoly$x <- c(sEndPoly$x,sEndPoly$x[1])
   sEndPoly$y <- c(sEndPoly$y,sEndPoly$y[1])

   grownPolygon <- list(x=sx*sEndPoly$x,y=sy*sEndPoly$y)
   
   return(list(grownPolygon=grownPolygon,volume=areaPolygon(grownPolygon)))
}

############ End of polygonGrow ############
