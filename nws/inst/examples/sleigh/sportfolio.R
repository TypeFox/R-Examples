library(grid)

# configurable parameters
numTasks <- 1000
numSamples <- 10000
numStocks <- 20
chunkSize <- 10

# randomly generate the mean and sd that describe each stock
set.seed(472)
smean <- rnorm(numStocks, mean=10.0, sd=1.0)
ssd <- rnorm(numStocks, mean=3.0, sd=0.5)
stocks <- data.frame(mean=smean, sd=ssd)

# this is the task function
fun <- function(numSamples, numStocks) {
  # generate the weights vector
  t <- runif(numStocks)
  w <- t / sum(t)

  # generate random stock returns matrix
  rnormWrapper <- function(i) rnorm(numSamples, mean=stocks$mean[[i]], sd=stocks$sd[[i]])
  s <- do.call(rbind, lapply(1:numStocks, rnormWrapper))

  # do the computation and return the results
  r <- drop(w %*% s)
  c(mean(r), var(r), SleighRank)
}

# create the plot window
grid.newpage()
vplay <- grid.layout(5, 3,
                     widths = unit(c(4, 1, 2),
                                   c('lines', 'null', 'lines')),
                     heights = unit(c(4, 12, 4, 1, 3),
                                    c('lines', 'null', 'lines', 'null', 'lines')))
pushViewport(viewport(layout=vplay))

pushViewport(viewport(layout.pos.col=1, layout.pos.row=2))
grid.text('Reward', x=unit(1, 'lines'), rot=90)

upViewport()
pushViewport(viewport(layout.pos.col=2, layout.pos.row=3))
grid.text('Risk', y=unit(1, 'lines'))

upViewport()
pushViewport(viewport(layout.pos.col=2, layout.pos.row=1, name='titleRegion'))
grid.text('Sequential Efficient Frontier', gp=gpar(fontsize=20, fontface='bold'))

upViewport()
pushViewport(viewport(yscale=c(9.0, 10.2), xscale=c(0.3, 1.3),
  layout.pos.col=2, layout.pos.row=2, name='plotRegion'))
grid.rect(gp=gpar(fill='light yellow'))
grid.segments(x0=unit(c(seq(0.4, by=0.2, length=5), rep(0.3, 5)), 'native'),
              y0=unit(c(rep(9.0, 5), seq(9.2, by=0.2, length=5)), 'native'),
              x1=unit(c(seq(0.4, by=0.2, length=5), rep(1.3, 5)), 'native'),
              y1=unit(c(rep(10.2, 5), seq(9.2, by=0.2, length=5)), 'native'),
              gp=gpar(col='gray', lty='dashed'))
grid.xaxis()
grid.yaxis()

upViewport()
pushViewport(viewport(layout.pos.col=2, layout.pos.row=4, name='barRegion'))
grid.rect(gp=gpar(fill='light yellow'))

seekViewport('barRegion')
bar <- rectGrob(x = unit(0, 'npc'), width=unit(0, 'npc'),
                gp=gpar(col='black', fill='red'), hjust=0)
grid.draw(bar)

upViewport()
pushViewport(viewport(layout.pos.col=2, layout.pos.row=5, name='subRegion'))
grid.rect(width=0.9, height=0.9, gp=gpar(col='white', fill='white'))
text <- textGrob(label=sprintf('Starting to execute %d tasks', numTasks))
grid.draw(text)

colors <- rainbow(3)

# do the work sequentially
reward <- vector()
risk <- vector()
rindx <- 1
SleighRank <- 1
while (rindx <= numTasks) {
  valueList <- lapply(rep(numSamples, min(chunkSize, numSamples - rindx + 1)), fun, numStocks)
  results <- unlist(valueList)
  numResults <- length(valueList)
  dim(results) <- c(length(valueList[[1]]), numResults)
  reward[rindx:(rindx + numResults - 1)] <- results[1,]
  risk[rindx:(rindx + numResults - 1)] <- results[2,]
  rindx <- rindx + numResults

  seekViewport('plotRegion')
  grid.points(results[2,], results[1,], pch=20, gp=gpar(cex=0.5, col=colors[1]))

  xlab <- if (rindx > numTasks) sprintf('Completed all %d Tasks', numTasks)
    else sprintf('Completed %d of %d Tasks', rindx-1, numTasks)

  seekViewport('subRegion')
  grid.rect(width=0.9, height=0.9, gp=gpar(col='white', fill='white'))
  text <- editGrob(text, NULL, label=xlab)
  grid.draw(text)

  seekViewport('barRegion')
  bar <- editGrob(bar, NULL, width=unit((rindx-1)/numTasks, 'npc'))
  grid.draw(bar)
}

# compute the efficient frontier
indx <- chull(risk, reward)
x <- risk[indx]
ix <- which.min(x)
y <- reward[indx]
iy <- which.max(y)
if (ix < iy) {
  i <- ix:iy
  x <- x[i]
  y <- y[i]
} else {
  i <- ix:(iy+length(indx))
  x <- c(x, x)[i]
  y <- c(y, y)[i]
}

seekViewport('plotRegion')
grid.points(x, y, gp=gpar(col='black'))
grid.lines(x=unit(x, 'native'), y=unit(y, 'native'), gp=gpar(col='black'))
