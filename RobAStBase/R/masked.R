setMethod("start", "ANY", function(x,...) stats::start(x,...))
setMethod("clip", "ANY", function(x1, x2, y1, y2) graphics::clip(x1, x2, y1, y2))
