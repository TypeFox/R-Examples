# Trying to convert to ggplot2

#  A work in progress!
require(ggplot2)
require(NNTbiomarker)
#qplot(tenYearDFS$RS, tenYearDFS$Recur, geom=geom_line())
ggplot(tenYearDFS, aes(x=RS, y=Recur, group=group)) +
  geom_line(aes(colour=group)) +
  geom_point(aes(colour=group), alpha=0.5, pch=15, fill="red", cex=10) +
  # theme_light()
  theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())

ROCplotsGG = function(data = data.frame(class=tenYearDFS$Recur, X=tenYearDFS$RS),
                    whichPlots=c("density", "raw", "ROC", "pv", "nnt", "nntRange"),
                    NNTlower=3, NNTupper=10,
                    N= 1000, prev=0.2, diffInSD=2,
                    ...) {
  seeThroughGrey = paste0("#404040", "88")
  seeThroughBlue =   paste0(rgb(0,0,.5), "22")
  seeThroughRed = paste0(rgb(0.9,0.1,.1), "22")
  Usr = function()par()$usr
  UsrX = function()Usr()[1:2]
  UsrY = function()Usr()[3:4]
  if(missing(data)) {  ## Simulate
    muD=0; muN=diffInSD; sd=1
    class = rbinom(N, 1, prev)
    X = rnorm(N, c(muD, muN)[1+class], sd)
    weights = rep(1, length(X))  ### Multiplicity
    data = data.frame(class=class, X=X, weights=weights)
  }
  else {
    if(is.null(data$weights))
      data$weights = rep(1, length(data$X))
  }
  data = data [order(data$X), ]
  class = data$class
  weights = data$weights
  X = data$X
  N = sum(weights)
  nD = sum(class * weights)
  nH = N - nD
  prevalence = nD/N  ### prevalence of BestToTreat.

  if(is.element(el = "density", set=whichPlots)) {
    ggplot(data=data) + aes(x=X, group=class) +
      geom_density(weights = weights/sum(weights)) +
      geom_rug(mapping = aes(col=class) )
  }

  cum1 = cumsum(class * weights)
  cum0 = cumsum((1-class) * weights)
  sensitivity = (nD - cum1)/ nD  # 1 - FN/nD = (TP/nD)
  specificity = cum0/nH       # TN/nH
  requiredSeSp = sesp.from.NNT(NNTlower, NNTupper, prev=prevalence)
  requiredSe = requiredSeSp[1]
  requiredSp = requiredSeSp[2]

  if(is.element(el = "ROC", set=whichPlots)) {
    ROCplot <- ggplot() + aes(x=1-specificity, y=sensitivity) +
      geom_line()
    if(!is.na(NNTlower)) {
      ROCplot <- ROCplot + guides(fill=
        guide_legend(title="ROC", label.position = "bottom",
                     label="acceptable\nregion",
                     box.col="yellow", bty="n", text.col="blue")
      )
      rect(xleft = UsrX()[1], xright = UsrX()[2], ybottom = UsrY()[1], ytop = requiredSe,
           col=seeThroughRed)
      rect(xleft = requiredSp, xright = UsrX()[2], ybottom = UsrY()[1], ytop = UsrY()[2],
           col=seeThroughRed)
      graphics::text(x=0, y=requiredSe, col="blue", labels = "required sensitivity",
                     xpd=NA, adj=c(0,0), cex=0.9)
      graphics::text(x=requiredSp, y=1, col="blue", labels = "required specificity",
                     xpd=NA, pos=3, cex=0.9)
    }
  }
}
