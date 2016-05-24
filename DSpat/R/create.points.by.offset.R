`create.points.by.offset` <-
function(lines,observations)
################################################################################
# For a set of observations with x,y locations on the line and a perpendicular
# distance, create a new observations dataframe with true x,y point locations.
#
# Arguments:
#  lines        - data frame of lines with the following structure
#                    label - unique label
#                    x0    - x coordinate of beginning of line
#                    y0    - y coordinate of beginning of line
#                    x1    - x coordinate of end of line
#                    y1    - y coordinate of end of line
#                    width - full width of transect around line (optional)
#                    ...   - any number of covariates
#
#  observations - data frame of observations with the following structure
#                    label - label linking it to a unique line
#                    x     - x coordinate
#                    y     - y coordinate
#                    ...   - any number of covariates
#
#  Value:  observations dataframe with true x,y locations
#
# 4/4/2008
# Jeff Laake
################################################################################
{
   number.of.lines=dim(lines)[1]
   obs=offset.points(lines[1,c("x0","x1","y0","y1")],observations[observations$label==lines$label[1]])
   for(i in 2:number.of.lines)
     obs=rbind(obs,offset.points(lines[i,c("x0","x1","y0","y1")],
             observations[observations$label==lines$label[i]]))
   return(obs)
}

