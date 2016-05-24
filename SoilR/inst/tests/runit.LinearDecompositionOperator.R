#
# vim:set ff=unix expandtab ts=2 sw=2:
test.castTimeMap2BoundLinDecompOp <- function(){
	tm=TimeMap.new(0,1,function(t){t})
	BoundLinDecompOp(tm)
}
