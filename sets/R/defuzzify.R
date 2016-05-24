cset_defuzzify <-
gset_defuzzify <-
function(x, method = c("meanofmax", "smallestofmax", "largestofmax", "centroid"))
    switch(match.arg(method),
           centroid = mean(x),
           smallestofmax = min(gset_peak(x)),
           meanofmax = mean(gset_peak(x)),
           largestofmax = max(gset_peak(x))
           )

