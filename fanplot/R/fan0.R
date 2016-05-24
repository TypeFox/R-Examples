fan0 <-
  function(data = NULL, data.type = "simulations", style = "fan", type = "percentile",
           probs = if(type=="percentile") seq(0.01, 0.99, 0.01) else c(0.5, 0.8, 0.95), 
           start = 1, frequency = 1, anchor = NULL, anchor.time=NULL,
           fan.col = heat.colors, alpha = if (style == "spaghetti") 0.5 else 1, 
           n.fan = NULL,
           ln = NULL, ln.col = if(style=="spaghetti") "gray" else NULL, 
           med.ln = if(type=="interval") TRUE else FALSE, med.col= "orange",
           rlab = ln, rpos = 4, roffset = 0.1, rcex = 0.8, rcol = NULL, 
           llab = FALSE, lpos = 2, loffset = roffset, lcex = rcex, lcol = rcol, 
           upplab = "U", lowlab = "L", medlab=if(type == "interval") "M" else NULL,
           n.spag = 30, 
           space = if(style=="boxplot") 1/frequency else 0.9/frequency, 
           add = TRUE, ylim = range(data)*0.8,...){
    if(add==TRUE)
      plot(data[,1], type="n", ylim=ylim, ...)
    fan(data = data, data.type=data.type, style = style, type = type,
    probs = probs, 
    start = start, frequency = frequency, anchor = anchor, anchor.time=anchor.time,
    fan.col = fan.col, alpha = alpha, 
    n.fan = n.fan,
    ln = ln, ln.col = ln.col, 
    med.ln = med.ln, med.col= med.col,
    rlab = rlab, rpos = rpos, roffset = roffset, rcex = rcex, rcol = rcol, 
    llab = llab, lpos = lpos, loffset = loffset, lcex = lcex, lcol = lcol, 
    upplab = upplab, lowlab = lowlab, medlab=medlab,
    n.spag = n.spag, 
    space = space, 
    add = FALSE)
  }
    