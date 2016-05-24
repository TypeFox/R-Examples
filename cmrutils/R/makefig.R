`makefig` <-
function(plotfun, ..., device = c('default', 'pdf', 'svg'),
         width = 8, height = 6, scale = pointsize/12, pointsize = 12, 
         filename = 'Rplot', family = 'Helvetica')
{
  device.params <- list()
  device.params$width <- scale*width
  device.params$height <- scale*height
  device.params$pointsize <- pointsize

  switch(match.arg(device),
         default =
         {
           do.call(dev.new, device.params)
           result <- plotfun(...)
         },
         pdf =
         {
           filename <- sub('(\\.pdf)?$', '.pdf', filename)
           device.params$file <- filename
           device.params$family <- family

           do.call(pdf, device.params)
           tryCatch(result <- plotfun(...),
                    finally = dev.off())
           embedFonts(filename, fontpaths = .pfbpath)
         },
         svg =
         {
           filename <- sub('(\\.svg)?$', '.svg', filename)
           device.params$file <- filename

           do.call(svg, device.params)
           tryCatch(result <- plotfun(...),
                    finally = dev.off())
         })

  return(invisible(result))
}
