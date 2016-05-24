getETOPO<-function(topo, glat =c(-90,90), glon = c(0, 360) )
    {
### topo is either ETOPO5 or ETOPO2
###   they must have attributes lat and lon
####  glat and glon are coordinates of lower left and upper right
### lat-lon pairs

###  returns a matrix in the format suitable for
### functions image, persp, etc....
        glon = RPMG::fmod(glon, 360)
               
        glon  =   range(glon)
        glat  =   range(glat)

        lon = attr(topo, 'lon')
        lat = attr(topo, 'lat')
        
        z=topo[ , rev(1:length(lat) )  ]
        colz = lon>=glon[1] & lon<=glon[2]
        rowz = lat>=glat[1] & lat<=glat[2]
        
        b5 = z[colz, rowz ]

        attr(b5,'lat') = lat[lat>=glat[1] & lat<=glat[2]]
        attr(b5, 'lon') = lon[lon>=glon[1] & lon<=glon[2]]
        
        return(b5)
    }
