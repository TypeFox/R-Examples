# Test Setup
 
if(FALSE) {
  # Not really needed, but can be handy when writing tests
  library("RUnit")
  library("Rclusterpp")
}

compare.hclust <- function(h1, h2) {
	checkEquals(h1$merge, h2$merge, msg="Agglomerations don't match")
	checkEquals(h1$height, h2$height, msg="Agglomeration heights are not equal")
	checkEquals(h1$labels, h2$labels, msg="Cluster labels do not match")
	checkEquals(h1$order, h2$order, msg="Cluster orders do not match")
}

# --- Test functions ---
 
test.hclust.ward <- function()
{
  d <- USArrests 
	h <- hclust((dist(d, method="euclidean")^2)/2.0, method="ward")
	r <- Rclusterpp.hclust(d, method="ward")
	compare.hclust(h, r)
}

test.hclust.average.euclidean <- function()
{
	d <- USArrests
	
	h <- hclust(dist(d, method="euclidean"), method="average")
	r <- Rclusterpp.hclust(d, method="average", distance="euclidean")
	compare.hclust(h, r)
}

test.hclust.average.manhattan <- function()
{
	d <- USArrests
	
	h <- hclust(dist(d, method="manhattan"), method="average")
	r <- Rclusterpp.hclust(d, method="average", distance="manhattan")
	compare.hclust(h, r)
}

test.hclust.average.maximum <- function()
{
	d <- USArrests
	
	h <- hclust(dist(d, method="maximum"), method="average")
	r <- Rclusterpp.hclust(d, method="average", distance="maximum")
	# USArrests clusters ambiguously under this metric, so only heights will match exactly
	checkEquals(r$height, h$height, msg="Agglomeration heights are not equal")
}

test.hclust.average.minkowski <- function()
{
	d <- USArrests
	
	h <- hclust(dist(d, method="minkowski"), method="average")
	r <- Rclusterpp.hclust(d, method="average", distance="minkowski")
	compare.hclust(h, r)
}

test.hclust.single.euclidean <- function()
{
	d <- USArrests
	
	h <- hclust(dist(d, method="euclidean"), method="single")
	r <- Rclusterpp.hclust(d, method="single", distance="euclidean")
	compare.hclust(h, r)
}

test.hclust.single.manhattan <- function()
{
	d <- USArrests
	
	h <- hclust(dist(d, method="manhattan"), method="single")
	r <- Rclusterpp.hclust(d, method="single", distance="manhattan")
	# USArrests clusters ambiguously under this metric, so only heights will match exactly
	checkEquals(r$height, h$height, msg="Agglomeration heights are not equal")
}

test.hclust.single.maximum <- function()
{
	d <- USArrests
	
	h <- hclust(dist(d, method="maximum"), method="single")
	r <- Rclusterpp.hclust(d, method="single", distance="maximum")
	# USArrests clusters ambiguously under this metric, so only heights will match exactly
	checkEquals(r$height, h$height, msg="Agglomeration heights are not equal")
}

test.hclust.single.minkowski <- function()
{
	d <- USArrests
	
	h <- hclust(dist(d, method="minkowski"), method="single")
	r <- Rclusterpp.hclust(d, method="single", distance="minkowski")
	compare.hclust(h, r)
}

test.hclust.complete.euclidean <- function()
{
	d <- USArrests
	
	h <- hclust(dist(d, method="euclidean"), method="complete")
	r <- Rclusterpp.hclust(d, method="complete", distance="euclidean")
	compare.hclust(h, r)
}

