"print.gvlma" <-
function(x, ...)
{
  NextMethod("print", x,...)
  display.gvlmatests(x)
}

