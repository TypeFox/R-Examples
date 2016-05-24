freqpolygon( ~ Exercise, data = StudentSurvey, breaks = seq(0, 45, by = 4), lwd = 3, 
            par.settings = col.whitebg(), 
            panel = function(x, ...) { 
  panel.xhistogram(x, ...);
  panel.freqpolygon(x, ...)}
             )

