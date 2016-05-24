###############################################################
#                                                             #
#   (c) Victor Maus <vwmaus1@gmail.com>                       #
#       Institute for Geoinformatics (IFGI)                   #
#       University of Muenster (WWU), Germany                 #
#                                                             #
#       Earth System Science Center (CCST)                    #
#       National Institute for Space Research (INPE), Brazil  #
#                                                             #
#                                                             #
#   R Package dtwSat - 2015-09-01                             #
#                                                             #
###############################################################


###############################################################
#### dtwSat TESTS


# Show query names
names(query.list)


# Perform twdtw alignment
alig = twdtw(query.list[["Soybean"]], waveletSmoothing(template), 
             weight = "logistic", alpha = 0.1, beta = 50, span=180, keep=TRUE)

is(alig, "dtwSat")
show(alig)
print(alig)


# Plot cost matrix paths
gp1 = plot(x=alig, type="path", normalize=TRUE, show.dist=TRUE)
gp1


# Plot alignment
gp2 = plot(alig, type="alignment", n=1, attr="evi", shift=0.5)
gp2


# Wavelet filter
sy = waveletSmoothing(x=template, frequency=8, wf = "la8", J=1, 
                      boundary = "periodic")


# Plot raw EVI and filtered EVI
gp = autoplot(sy, facets = NULL) + xlab("Time")
gp

# Plot all filtered bands
evi = merge(Raw=zoo(template$evi), Wavelet=zoo(sy$evi))
gp = autoplot(evi, facets = NULL) + xlab("Time")
gp


# Perform twdtw to query list 
malig = mtwdtw(query.list, template, weight = "logistic", 
               alpha = 0.1, beta = 100)
class(malig)
getAlignments(malig)
      
# Classify interval
best_class = classfyIntervals(x=malig, from=as.Date("2009-09-01"), 
                              to=as.Date("2013-09-01"), by = "6 month",
                              normalized=TRUE, overlap=.5, threshold=Inf)
best_class


malig = mtwdtw(query.list, template, weight = "logistic", 
               alpha = 0.1, beta = 100)
 
gp = plotClassify(x=malig, from=as.Date("2009-09-01"),  
              to=as.Date("2013-09-01"), by = "6 month",
              normalized=TRUE, overlap=.5)
gp
# ggsave("classify.png", plot=gp, width = 8.9, height=5.9/1.5, units="in",
#         family="Helvetica", type = "cairo-png")


# Plot cost matrix 
alig = twdtw(query.list[["Soybean"]], template, weight = "logistic", 
             alpha = 0.1, beta = 100, keep=TRUE)
 

gp = plotCostMatrix(x=alig, matrix.name="timeWeight")
gp

gp = plotCostMatrix(x=alig, matrix.name="localMatrix")
gp
 
gp = plotCostMatrix(x=alig, matrix.name="costMatrix")
gp

