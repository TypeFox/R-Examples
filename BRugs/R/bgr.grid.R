"bgrGrid" <-
function(node, bins = 50)
# Calculate grid of points at which to evaluate bgr statistic
{
    sampleSize <- samplesSize(node)
    beg <- samplesGetBeg()
    end <- min(c(samplesGetEnd(), modelIteration()))
    numChains <- samplesGetLastChain() - samplesGetFirstChain() + 1
    sampleSize <- sampleSize %/% numChains
    beg <- end - (sampleSize * samplesGetThin() - 1)
    delta <- sampleSize %/% bins
    grid <- ((1 : (bins - 1)) * delta) + beg
    grid <- c(grid, end)
    grid
}
