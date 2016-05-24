#download shapefile
usMap <- readShapePoly(system.file("shapes/usMap.shp", package="mapStats")[1])

#create synthetic survey dataset



state_codes <- c('AL','AK','AZ','AR','CA','CO','CT','DE','DC','FL','GA','HI','ID','IL',
                 'IN','IA','KS','KY','LA','ME','MD','MA','MI','MN','MS','MO','MT','NE',
                 'NV','NH','NJ','NM','NY','NC','ND','OH','OK','OR','PA','RI','SC','SD',
                 'TN','TX','UT','VT','VA','WA','WV','WI','WY')

surveydata <- data.frame(state=factor(rep(rep(state_codes, 
                         times=sample(x=c(30,40,60), size=length(state_codes), replace=TRUE)), 
                         times=2)))
surveydata$year <- factor(rep(c(2009, 2010), each=nrow(surveydata)/2))
surveydata$educ <- factor(sample(c("HS","BA","Grad"), p=c(.4,.45,.15), replace=TRUE, size=nrow(surveydata)))
surveydata$sex <- factor(sample(c("Male","Female"), replace=TRUE, size=nrow(surveydata)))

surveydata$obs_weight <- runif(n=nrow(surveydata), min=0.8, max=1.5)

#two income distributions
surveydata$income <- 100000*rbeta(n=nrow(surveydata), 
                                  shape1=ifelse(surveydata$year==2009, 2, 1.5),
                                  shape2=ifelse(surveydata$year==2010, 10, 11))
surveydata$income_ge30k <- ifelse(surveydata$income >=30000, 100, 0)

surveydata$age <- round(100*rbeta(n=nrow(surveydata), shape1=2, shape2=4))

#these state and year combinations will not be shaded if they are missing entirely
surveydata[ surveydata$state == "NV" & surveydata$year == 2009, c("income","income_ge30k")] <- NA
surveydata[ surveydata$state == "OH" & surveydata$year == 2010, c("income","income_ge30k")] <- NA

#create some missing values in other variables
surveydata[ sample(1:nrow(surveydata), 15), c("income","income_ge30k")] <- NA
surveydata[ sample(1:nrow(surveydata), 15), "obs_weight"] <- NA
surveydata[ sample(1:nrow(surveydata), 15), "year"] <- NA
surveydata[ sample(1:nrow(surveydata), 15), "state"] <- NA


head(surveydata)

#Calculate weighted mean and 50th, and 60th quantiles of variable income,
#by state, survey year, and education.  Display using red sequential color palette with 6 groups.
#In the titles, rename 'income' with 'household income', and  'educ' with 'education'.
#Resize the strip text.

#print each statistic on a separate page, using separate color breaks for each.  Will print 2 pages.       


mapStats(d=surveydata, var="income", stat=c("mean", "quantile"), quantiles=.5, wt.var="obs_weight",
               d.geo.var="state", by.var=c("year","educ"), map.file=usMap, map.label=FALSE,
               par.strip.text = list(cex = .7))

#In the titles, rename 'income' with 'household income', and  'educ' with 'education'.
#Change the number of color breaks from 4 to 6, add state labels.               

mapStats(d=surveydata, var="income", d.geo.var="state", by.var=c("year","educ"),
         wt.var="obs_weight", map.file=usMap, map.geo.var="state",
         stat=c("mean","quantile"), quantiles=.5,
         paletteName="Reds", ngroups=6, var.pretty="household income",
         by.pretty=c("year","education"), map.label=TRUE,
         par.strip.text = list(cex = .7), cex.label=.5)



#print statistics on 2 pages, with separate color breaks, with 1 column and 3 rows, 
#and use layout argument to position panels within each plot properly .
#eliminate labels because too crowded.     

mapStats(d=surveydata, var="income", d.geo.var="state", by.var=c("year"),
         wt.var="obs_weight", map.file=usMap,
         map.geo.var="state", stat=c("mean","quantile"),
         quantiles=c(.25, .5, .75), paletteName="Reds", ngroups=6, 
         var.pretty="household income", map.label=FALSE,
         num.col=1, num.row=3, layout=c(2,1))

#plot all four statistics over all four years, on one page with 2 rows and 2 columns,
#using vertical order. Shrink labels so they fit better. 
#shrink plot titles so all fit. 


mapStats(d=surveydata, var="income", d.geo.var="state", by.var=c("year"),
         wt.var="obs_weight", map.file=usMap,
         map.geo.var="state", stat=c("mean","quantile"),
         quantiles=c(.25, .5, .75), paletteName="Reds", ngroups=6, 
         map.label=FALSE, var.pretty="household income",
         num.col=2, num.row=2, layout=c(2,1), cex.label=.5, horizontal.fill=FALSE)


#do the above example, except with a different palette for each to distinguish them 
#since separate breaks

mapStats(d=surveydata, var="income", d.geo.var="state", by.var=c("year"),
         wt.var="obs_weight", map.file=usMap,
         map.geo.var="state", stat=c("mean","quantile"),
         quantiles=c(.25, .5, .75), paletteName=c("Reds","Greens","Blues","Greys"), ngroups=6, 
         map.label=FALSE,  var.pretty="household income", 
         num.col=2, num.row=2, layout=c(2,1), cex.label=.5, horizontal.fill=FALSE)




#use the above with the same palette except with the same color breaks for all statistics,

mapStats(d=surveydata, var="income", d.geo.var="state", by.var=c("year"),
         wt.var="obs_weight", map.file=usMap,
         map.geo.var="state", stat=c("mean","quantile"),
         quantiles=c(.25, .5, .75), paletteName=c("Reds","Greens","Blues","Greys"), ngroups=6, 
         var.pretty="household income", map.label=FALSE, 
         num.col=2, num.row=2, layout=c(2,1), cex.label=.5, horizontal.fill=FALSE,
         separate=0)


#do the same except use jenks style to calculate the breaks (classIntervals)

mapStats(d=surveydata, var="income", d.geo.var="state", by.var=c("year"),
         wt.var="obs_weight", map.file=usMap,
         map.geo.var="state", stat=c("mean","quantile"),
         quantiles=c(.25, .5, .75), paletteName=c("Reds","Greens","Blues","Greys"), ngroups=6, 
         var.pretty="household income", map.label=FALSE, 
         num.col=2, num.row=2, layout=c(2,1), cex.label=.5, horizontal.fill=FALSE,
         separate=0, style="jenks")

#an example doing the above with HCL instead to override RColorBrewer
#differing numbers of groups, ngroups overridden
#three palettes will be reused over four plots

red_hcl <- rev(colorspace::sequential_hcl(n=5, h=0))
green_hcl <- rev(colorspace::sequential_hcl(n=7, h=120))
blue_hcl <- rev(colorspace::sequential_hcl(n=3, h=240))

mapStats(d=surveydata, var="income", d.geo.var="state", by.var=c("year"),
         wt.var="obs_weight", map.file=usMap,
         map.geo.var="state", stat=c("mean","quantile"),
         quantiles=c(.25, .75), colorVec=list(red_hcl, green_hcl, blue_hcl),
         var.pretty="household income", map.label=FALSE, 
         num.col=2, num.row=2, layout=c(2,1), cex.label=.5, horizontal.fill=FALSE,
         separate=1, style="jenks")


#Example with additional parameters to control the spacing of the shapefile within the panel

#plot bounds
bounds <- attributes(usMap)$bbox
ranges <- abs(bounds[,2]-bounds[,1])
new_bounds <- bounds
new_bounds[,1] <- new_bounds[,1] - 0.1*ranges
new_bounds[,2] <- new_bounds[,2] + 0.1*ranges
xbox <- new_bounds[1,]
ybox <- new_bounds[2,]


mapStats(d=surveydata, var="income", d.geo.var="state", by.var=c("year"),
         wt.var="obs_weight", map.file=usMap,
         map.geo.var="state", stat=c("mean"), paletteName="Reds", 
         ngroups=6, var.pretty="household income",
         map.label=TRUE, xlim=xbox, ylim=ybox)



#To calculate percentages of class variables, create an indicator variable, calculate its mean,
#and override the default titles to say you are calculating the percentage.  Here we illustrate by
#calculating the percent of respondents by state that have income above 30,000.


mapStats(d=surveydata, var="income_ge30k", 
         wt.var="obs_weight", map.file=usMap,
         d.geo.var="state", map.geo.var="state",
         stat="mean", paletteName="Reds", 
         titles="Percent of respondents with income\nat least $30,000")

#Do the same as above but add additional visual elements
#Normally you may want to overlay with a different shapefile showing other elements,
#such as rivers, or something else relevant to the question
#You can overlay with points to label specific areas on the map, or add a compass point

#subset the map, plot a few green
state_subset <- usMap[ usMap$state %in% c("TX", "AK", "CA"), ]
shaded_green <- list("sp.polygons", state_subset, fill="green")

#add a red triangle on the map
attributes(usMap)$bbox

red_triangle <- list("sp.points", matrix(c(0, 2000000), ncol=2), 
                      pch=2, col="red", cex=4, lwd=5)

#combine the two things above in a list

map_overlay <- list(shaded_green, red_triangle)


mapStats(d=surveydata, var="income_ge30k", 
         wt.var="obs_weight", map.file=usMap,
         d.geo.var="state", map.geo.var="state",
         stat="mean", paletteName="Reds", 
         titles="Percent of respondents with income\nat least $30,000",
         sp_layout.pars=map_overlay)



#color multiple regions:  create code for Census region

west <- c("AK","AZ","WA","MT","OR","ID","WY","CA","NV","UT","CO","NM","HI")
midwest <- c("ND","MN","SD","NE","IA","WI","MI","OH","IN","IL","MO","KS")
neast <- c("PA","NY","NJ","CT","RI","MA","VT","NH","ME")
south <- c("TX","OK","AR","LA","MS","KY","TN","AL","FL","GA","SC","NC","VA","WV","DC","MD","DE")

surveydata$region[ surveydata$state %in% west ] <- "West"
surveydata$region[ surveydata$state %in% midwest ] <- "Midwest"
surveydata$region[ surveydata$state %in% neast ] <- "Northeast"
surveydata$region[ surveydata$state %in% south ] <- "South"

usMap$region[ usMap$state %in% west ] <- "West"
usMap$region[ usMap$state %in% midwest ] <- "Midwest"
usMap$region[ usMap$state %in% neast ] <- "Northeast"
usMap$region[ usMap$state %in% south ] <- "South"


#calculate mean and median by Census region, using a combined scale for 
#all statistics (separate = 0), and stack vertically
#multiple statistics and no by variable

mapStats(d=surveydata, var=c("income","age"), wt.var="obs_weight", map.file=usMap,
         d.geo.var="region", map.geo.var="region",
         stat=c("mean", "quantile", "var"), quantiles=c(.5),
         paletteName=c("Reds"), ngroups=3, var.pretty=c("household income","age"),
         map.label=FALSE, cex.label=.6, separate=1,
         horizontal.fill=FALSE, num.row=3)


#use of separate with multiple statistics
#separate=1 allows you to compare the medians, so plotvars=FALSE to group by statistic
#Greens ignored because only 2 variables and separate=2
#plotbyvar can be changed

mapStats(d=surveydata, var=c("income","age"), wt.var="obs_weight", map.file=usMap,
         d.geo.var="region", map.geo.var="region", by.var=c("year"),
         stat=c("mean", "quantile"), quantiles=c(.25, .75),
         paletteName=c("Reds","Blues","Greens"), ngroups=c(3,4), var.pretty=c("household income","age"),
         map.label=FALSE, cex.label=.6, separate=2, 
         horizontal.fill=FALSE, num.row=3, layout=c(2,1))

#now Greens used because separate breaks by statistic and there are 3
#groups will be 3, 4, then 3
#here this doesn't make much sense because the age and income variables have different ranges
#so using the same colors is silly but it could


mapStats(d=surveydata, var=c("income","age"), wt.var="obs_weight", map.file=usMap,
         d.geo.var="region", map.geo.var="region", by.var=c("year"),
         stat=c("mean", "quantile"), quantiles=c(.25, .75),
         paletteName=c("Reds","Blues","Greens"), ngroups=c(3,4), var.pretty=c("household income","age"),
         map.label=FALSE, cex.label=.6, separate=1,
         horizontal.fill=FALSE, num.row=2, layout=c(2,1))

