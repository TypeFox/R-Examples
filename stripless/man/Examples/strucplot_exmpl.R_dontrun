
 # Compare standard xyplot vs. strucplot of quakes data in datasets package.
 #  Cut depth into 6 groups and magitude into 5.

 require(datasets)
 # Note that as.table = TRUE is used to make depths increase
 # down the page. For strucplot(), this is the default.

 xyplot(lat ~ long|cut(mag,5)*cut(depth,6),  data = quakes, col="blue",
  as.table = TRUE,type = c("g","p"))

 # Compare to:
 strucplot(lat ~ long|cut(mag,5)*cut(depth,6), data = quakes, col="blue",
    type = c("g", "p"))

 # Visualizing designs:

 # A half fraction of a 2^5 (a 2^(5-1)) design

 # Build the design matrix
  ff <- do.call(expand.grid,rep(list(c(-1,1)),4))
  ff[[5]] <- do.call(mapply,c(FUN = prod,ff))
   names(ff) <- LETTERS[1:5]

 # Show the design
 strucplot(ff)

 # Plotting a 2 level design with a center point

 # Add a center point to ff and plot
  ffCenter <- rbind(ff,rep(0,5))
  \donttest{strucplot(ffCenter)}

 # Use center = TRUE for a more compact display and show legend below.
  print(strucplot(ffCenter, center = TRUE),legendLoc = "bottom")

 # The "npk" data. See help("npk") for details.
 # Visualize design with blocks the vertical factor and the rest horizontal
  strucplot(npk[,-5], xyLay = list(x = 2:4, y =1))

  # Plot the yield
  strucplot(~yield |., xyLay = list(x=2:4, y=1),data = npk, col = "darkblue",
    panel = function(...){
      panel.grid(h = -1, v = 0)
      panel.xyplot(...)}
    )

 # It may be more informative to plot bars instead of points.
 # See help(panel.bars) for details.
 #
 # Note also "shortcut" ways to specify the xyLayout
  strucplot(~yield |., xyLay = list(x=2:4),data = npk,
            panel = panel.bars)

 # Include a conditioning variable in the formula to reduce the
 # dimensionality of conditioning. Show legend on right of plot.
  print(strucplot(yield ~ N|., xyLay = 2:3, data = npk,
            panel = panel.bars), legendLoc = "right")

 # Use the horizontal = TRUE argument of panel.bars to plot the bars
 # horizontally. The left and right hand sides of the formula must also
 # be switched for 2-sided formulas (not for 1-sided).
  strucplot( N ~ yield |., xyLay = list(y=1), data = npk,
            panel = panel.bars, horizontal = TRUE)

 # Fit a linear model with all main effects and 2 factor interactions in N,P,K
 # and plot the fits, using the "newdata" argument to plot predictions at
 # non-design points).
  require("stats")
    npk.aov <- aov(yield ~ block + (N+P+K)^2, data = npk)
    full <- do.call(expand.grid,lapply(npk[,-5],levels))
    plot(strucplot(npk.aov,  xyLay = list(x = 2:4),panel = panel.bars,
            newdata = full),legendLoc = "bottom")

   # Compare to a grouped plot:
   ypred <- predict(npk.aov, new = full)
    plot(
      strucplot(ypred ~ N|K*block, groups = full$P, data = full,
        panel= function(x,y, groups, subscripts,cex=1.25,...){
          panel.grid(h=-1, v=0)
          panel.superpose(x,y,cex= cex, type = c("p","l"),...,
                    panel.groups = panel.xyplot,
                    groups=groups, subscripts =subscripts)},
        auto.key = list(points=FALSE,lines=TRUE, columns = 2,
                title = "P",cex.title=1), ylab = "Predicted Response" ),
      legendLoc = "right")

## Cleanup
rm(full, npk.aov, ypred,ff,ffCenter)

