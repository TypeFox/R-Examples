#Testing using package testthat for detailed error messages
library(testthat)

#Get the testing function applied to compare the two venn diagram objects
source("testFunction.R");

#load in the reference plot data
load("data/plotsFive.rda");

#Suppress plotting for sanity
options(device=pdf());

#initialize the testing list of venn diagrams
venn.test <- list();

#Colour

venn.test <- c(venn.test,list(draw.quintuple.venn(
    area1 = 301,
    area2 = 321,
    area3 = 311,
    area4 = 321,
    area5 = 301,
    n12 = 188,
    n13 = 191,
    n14 = 184,
    n15 = 177,
    n23 = 194,
    n24 = 197,
    n25 = 190,
    n34 = 190,
    n35 = 173,
    n45 = 186,
    n123 = 112,
    n124 = 108,
    n125 = 108,
    n134 = 111,
    n135 = 104,
    n145 = 104,
    n234 = 111,
    n235 = 107,
    n245 = 110,
    n345 = 100,
    n1234 = 61,
    n1235 = 60,
    n1245 = 59,
    n1345 = 58,
    n2345 = 57,
    n12345 = 31,
    category = c("A", "B", "C", "D", "E"),
    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    cat.cex = 2,
    margin = 0.05,
    cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 
    1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
    ind = TRUE
    )))

testNames <- c('colour');

#Strip the polygons of their x and y values. They have equivalent information in their params field

for(i in 1:length(venn.test)){
	for(j in 1:length(venn.test[[i]])){
		if(class(venn.test[[i]][[j]])[1] == "polygon"){
			venn.test[[i]][[j]]$x <- NULL;
			venn.test[[i]][[j]]$y <- NULL;
		}
	}
}

#Loop over all of the test cases
for(i in 1:length(venn.test)){
	test_that(paste("Case",testNames[i],"of five categories"),
	{
		for(j in 1:length(venn.test[[i]])){
			expect_that(venn.test[[i]][[j]],is_identical_without_name(venn.plot[[i]][[j]],maxLength=3));
		}
	})
}

#Reaches here only if error is not thrown beforehand
print("Five category tests complete. No discrepancies found");
