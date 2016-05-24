# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


vplot = function (x, y, labels = NULL, colors = NULL, points = NULL, meansonly = FALSE, ellipsesd = 0, 
                  add = FALSE, alternateaxes = FALSE, xsampa = FALSE, logaxes = FALSE, ...){
  if (min(table (labels)) < 2 & ellipsesd>0) 
    stop ('At least 3 tokens per category are required to plot ellipses.')
  if (logaxes & min(x,y) <= 0) stop ('Log axes are incompatible with negative plotting values.')
  
  cl = match.call()
  matched = match(c('labels', 'meansonly', 'ellipsesd', 'add', 'colors', 
                    'alternateaxes', 'xsampa','points','logaxes'),names(cl),0)
  p = cl[-matched] ## for plot()
  args = sapply (2:length(p), function(x) p[[x]])
  names(args) = names(p)[-1]
  
  if (alternateaxes){tmp = x; x = y; y = tmp; 
                     tmp = args$x; args$x = args$y; args$y = tmp}
  allx = x; ally = y; alllabels = labels;
  
  if (meansonly) oranges = c(range(x), range(y))
  if (meansonly & !is.null(labels)){ 
    x = tapply (x, labels, mean); args$x = call('(', x);
    y = tapply (y, labels, mean); args$y = call('(', y)
    labels = names(y)
  }
  if (meansonly & is.null(labels)) 
    stop ('Mean vowel category plotting only possible if labels are given.')
  
  if (logaxes & !add) args$log = call('quote', 'xy')
  
  #######
  if (match("xlim", names(args), 0)>0) xlim = eval(args[[match("xlim", names(args), 0)]])
  if (match("ylim", names(args), 0)>0) ylim = eval(args[[match("ylim", names(args), 0)]])
  
  if (match("xlim", names(args), 0)==0 & !logaxes)
    xlim = range(allx)+c(-abs(diff(range(allx)))/20,abs(diff(range(allx)))/20)  
  if (match("ylim", names(args), 0)==0 & !logaxes)
    ylim = range(ally)+c(-abs(diff(range(ally)))/20,abs(diff(range(ally)))/20)
  
  if (match("xlim", names(args), 0)==0 & logaxes)
    xlim = range(allx)*c(.9,1.1)  
  if (match("ylim", names(args), 0)==0 & logaxes)
    ylim = range(ally)*c(.9,1.1)
  
  if (match("xlim", names(args), 0)==0 & !logaxes & meansonly & ellipsesd == 0)
    xlim = range(x)+c(-abs(diff(range(x)))/20,abs(diff(range(x)))/20)  
  if (match("ylim", names(args), 0)==0 & !logaxes & meansonly & ellipsesd == 0)
    ylim = range(y)+c(-abs(diff(range(y)))/20,abs(diff(range(y)))/20)
  
  if (match("xlim", names(args), 0)==0 & logaxes & meansonly & ellipsesd == 0)
    xlim = range(x)*c(.9,1.1)  
  if (match("ylim", names(args), 0)==0 & logaxes & meansonly & ellipsesd == 0)
    ylim = range(y)*c(.9,1.1)  
  
  if (alternateaxes){ xlim=rev(xlim); ylim=rev(ylim);}
  args$xlim = call('c', xlim[1],xlim[2]); args$ylim = call('c', ylim[1],ylim[2]);
  
  if (match("cex.axis", names(args),0)==0) args$cex.axis = call('(', 1.1)
  if (match("cex.lab", names(args),0)==0) args$cex.lab = call('(', 1.1)
  if (match("cex", names(args),0)==0 & meansonly) args$cex = call('(', 3)
  if (match("cex", names(args),0)==0 & !meansonly) args$cex = call('(', 1.2)
  
  if (match("lwd", names(args), 0)==0) lwd = 2
  if (match("lwd", names(args), 0)>0) lwd = args[[match("lwd", names(args), 0)]]
  
  if (match("xlab", names(args),0)==0 & !alternateaxes) args$xlab = call('quote', 'F1 (Hz)')
  if (match("ylab", names(args),0)==0 & !alternateaxes) args$ylab = call('quote', 'F2 (Hz)')
  if (match("xlab", names(args),0)==0 & alternateaxes) args$xlab = call('quote', 'F2 (Hz)')
  if (match("ylab", names(args),0)==0 & alternateaxes) args$ylab = call('quote', 'F1 (Hz)')
  #######
  
  vlevels = levels (as.factor (alllabels))
  vnums = as.numeric (as.factor(alllabels))
  
  if (is.null(colors)) colors = rep(colors()[c(24,506,118,610,30,124,556,258,290,151,84,657,404)],10)
  if (!meansonly) cols = colors[vnums]
  if (!meansonly & length(colors)==length(x)) cols = colors
  if (meansonly) cols = colors
  args$col = call ('[', cols)
  
  if (!is.null(points)){
    points = rep(points, 100)    
    args$pch = call ('(', quote(points[vnums]))
  }
  if (xsampa) args$pch = call ('xsampatoIPA', quote(labels))
  if (!add) do.call ('plot', args)
  if (add) do.call ('points', args)
  
  if (ellipsesd > 0){
    for (i in 1:length(vlevels)){  ##fix density so that you can specify number of points not spacing
      if (!logaxes)sdellipse (cbind (allx[alllabels==vlevels[i]],ally[alllabels==vlevels[i]]), 
                              stdev = ellipsesd, col = colors[i],lwd=lwd) 
      
      if (logaxes){ tmp = sdellipse (log(cbind (allx[alllabels==vlevels[i]],ally[alllabels==vlevels[i]])), 
                                     stdev = ellipsesd, show = F); lines (exp(tmp), col = colors[i],lwd=lwd)}
    }
  }   
}

