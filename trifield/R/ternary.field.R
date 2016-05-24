ternary.field = function(grid, vals, dim.out = NULL)
{
with(grid,
{
if ( is.null(dim.out) ) dim.out = c(max(xi), max(yi))
q = array(NA, dim = dim.out)
q[cbind(xi, yi)] = vals
x.out = seq(0, 2 / sqrt(3), length = dim.out[1])
y.out = seq(0, 1, length = dim.out[2])
res = list(x = x.out, y = y.out, z = q)
class(res) = "trifield"
return(res)
})
}
