project2line=function(obs.ppp,lines.psp)
################################################################################
# Projects point process contained in strips to centerline. This is the
# inverse of the create.points.by.offset function.
#
# Arguments:
#
#   obs.ppp - point process for observations in a strip
#   line.psp   - line segment process with label field
#
# Value:  list with elements
#    projection- dataframe of projected locations on the lines
#
# Author: Devin Johnson
################################################################################
{
  proj.obs=function(label, obs.ppp, lines.psp)
  {
      pts = obs.ppp[obs.ppp$label==label]
      ends = lines.psp$ends[lines.psp$label==label,]
      return(dist2line(pts,ends)$projection)
  }
  return(do.call('rbind', lapply(lines.psp$label,proj.obs,
                   obs.ppp=obs.ppp,lines.psp=lines.psp)))
}

