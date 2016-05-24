data(sampson)
set.seed(8675309)
require("ggplot2")
require("grid")
theme_set(theme_bw())

result <-
  mmsb.collapsed.gibbs.sampler(sampson$SAMPLK3,
                               K = 3, num.iterations=30,
                               alpha = 0.1,
                               burnin = 20L,
                               beta.prior = list(1, diag(5, 3) + 1))
print(result$document_expects)
cat("Time to plot the results...\n")
memberships <- with(result, t(document_sums) / colSums(document_sums))
colnames(memberships) <- paste("theta", 1:3, sep=".")
memberships <- cbind(data.frame(name=colnames(sampson$SAMPLK3)),
                     memberships)
memberships$colors <- with(memberships, rgb(theta.3, theta.1, theta.2))

center <- c(sqrt(3) / 4, sqrt(3) / 4)
angles <- with(memberships, atan2(theta.3 - center[2], theta.1 - center[1]))

angle.diffs <- tapply(angles, as.factor(angles), function(x) {
  pi/4 + seq(from=-(length(x) - 1) / 2,
             to=(length(x) - 1) / 2,
             length.out=length(x)) * pi / 6
})
angles[order(angles)] <- unlist(angle.diffs)

plot.1 <- ggplot(data = memberships) +
  geom_segment(aes(x = c(0, 0, 1),  y = c(0, 1, 0),
                   xend = c(0, 1, 0), yend = c(1, 0, 0))) +
  geom_point(aes(x = theta.1, y = theta.3, color = colors)) +
  scale_colour_manual(values = structure(memberships$colors, names = memberships$colors)) +
  scale_x_continuous(breaks=seq(0, 1, length.out=5),
                     limits = c(-0.25, 1.25)) +
  scale_y_continuous(breaks=seq(0, 1, length.out=5),
                     limits = c(-0.25, 1.25)) +
  geom_text(aes(x=theta.1, y=theta.3, label=name, colour = colors,
                angle=angles * 180 / pi), 
            data = memberships,
            size=2, hjust=-0.5) +
  ggtitle("Latent positions") +
  xlab(expression(theta[1])) +
  ylab(expression(theta[3])) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")
                   
  ## Block relations plot
ratio <- with(result, blocks.pos / (blocks.pos + blocks.neg))
total <- with(result, blocks.pos + blocks.neg)
data <- as.data.frame(cbind(Probability=as.numeric(ratio),
                            Count=as.numeric(total),
                            Column=rep(1:3, each=3),
                            Row=rep(1:3, times=3)))
plot.2 <- qplot(Column, Row, main="Block relations",
                size=Count, colour=Probability, data=data) +
  scale_size(range=c(7,15)) +
  scale_x_continuous(breaks=1:3, limits=c(0.5, 3.5)) +
  scale_y_reverse(breaks=1:3, limits=c(3.5, 0.5))

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(plot.1, vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(plot.2, vp=viewport(layout.pos.row=1,layout.pos.col=2))
