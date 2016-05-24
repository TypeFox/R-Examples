interaction.positioned <- function(..., drop = FALSE, sep = ".",
                                   b.offset=0,
                                   b.scale=.1) {
  dotdotdot <- list(...)
  if (length(dotdotdot) != 2)
    stop("interaction.positioned requires exactly two factors")
  a <- dotdotdot[[1]]
  b <- dotdotdot[[2]]
  ## interaction as distributed does not use Methods
  ## the 'get()' is here as defensive programming
  a.b <- if.R(r=get("interaction","package:base")(a, b, drop=drop, sep=sep),
              s=get("interaction","splus")(a, b, drop=drop))
  a.b <- ordered(a.b,
                 levels=matrix(
                   levels(a.b),
                   length(levels(b)),
                   length(levels(a)),
                   byrow=TRUE))
  new.pos.b <- (position(b)+b.offset)*b.scale
  new.pos.b <- new.pos.b - mean(new.pos.b)  ## center on 0
  position(a.b) <- as.vector(outer(new.pos.b, position(a), "+"))
  a.b
}
