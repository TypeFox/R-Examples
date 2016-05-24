#Testing using package testthat for detailed error messages
library(testthat)

#Get the testing function applied to compare the two venn diagram objects
source("testFunction.R");

#load in the reference plot data
load("data/plotsThree.rda");

#Suppress plotting for sanity
options(device=pdf());

#initialize the testing list of venn diagrams
venn.test <- list();

#Default

venn.test <- c(venn.test,list(draw.triple.venn(65, 75, 85,
 35, 15, 25, 5, c("First", "Second", "Third"))))

#Default and Colour

venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 65,
    area2 = 75,
    area3 = 85,
    n12 = 35,
    n23 = 15,
    n13 = 25,
    n123 = 5,
    category = c("First", "Second", "Third"),
    fill = c("blue", "red", "green"),
    lty = "blank",
    cex = 2,
    cat.cex = 2,
    cat.col = c("blue", "red", "green")
    )))

#001
venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 4,
    area2 = 3,
    area3 = 4,
    n12 = 2,
    n23 = 2,
    n13 = 2,
    n123 = 1,
    category = c('A', 'B', 'C'),
    #category = c('C','B','A'),
    fill = c('red', 'blue', 'green'),
    cat.col = c('red', 'blue', 'green'),
    cex = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
    cat.cex = c(1,2,3),
    euler = TRUE,
    scaled = FALSE
    )))

#010
venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 3,
    area2 = 3,
    area3 = 4,
    n12 = 1,
    n23 = 2,
    n13 = 2,
    n123 = 1,
    category = c('A', 'B', 'C'),
    fill = c('red', 'blue', 'green'),
    cat.col = c('red', 'blue', 'green'),
    cex = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
    cat.cex = c(1,2,3),
    euler = TRUE,
    scaled = FALSE
    )))
    
#011A
venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 3,
    area2 = 2,
    area3 = 4,
    n12 = 1,
    n23 = 2,
    n13 = 2,
    n123 = 1,
    category = c('A', 'B', 'C'),
    fill = c('red', 'blue', 'green'),
    cat.col = c('red', 'blue', 'green'),
    cex = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
    cat.cex = c(1,2,3),
    euler = TRUE,
    scaled = FALSE
    )))
    
#011O
venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 3,
    area2 = 3,
    area3 = 3,
    n12 = 1,
    n23 = 2,
    n13 = 2,
    n123 = 1,
    category = c('A', 'B', 'C'),
    fill = c('red', 'blue', 'green'),
    cat.col = c('red', 'blue', 'green'),
    cex = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
    cat.cex = c(1,2,3),
    euler = TRUE,
    scaled = FALSE
    )))

#012AA
venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 2,
    area2 = 2,
    area3 = 4,
    n12 = 1,
    n23 = 2,
    n13 = 2,
    n123 = 1,
    category = c('A', 'B', 'C'),
    fill = c('red', 'blue', 'green'),
    cat.col = c('red', 'blue', 'green'),
    cex = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
    cat.cex = c(1,2,3),
    euler = TRUE,
    scaled = FALSE
    )))

#021AA
venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 3,
    area2 = 1,
    area3 = 3,
    n12 = 1,
    n23 = 1,
    n13 = 2,
    n123 = 1,
    category = c('A', 'B', 'C'),
    fill = c('red', 'blue', 'green'),
    cat.col = c('red', 'blue', 'green'),
    cex = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
    cat.cex = c(1,2,3),
    euler = TRUE,
    scaled = FALSE
    )))

#022AAAO
venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 2,
    area2 = 1,
    area3 = 3,
    n12 = 1,
    n23 = 1,
    n13 = 2,
    n123 = 1,
    category = c('A', 'B', 'C'),
    fill = c('red', 'blue', 'green'),
    cat.col = c('red', 'blue', 'green'),
    cex = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
    cat.cex = c(1,2,3),
    euler = TRUE,
    scaled = FALSE
    )))

#022AAOO
venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 2,
    area2 = 2,
    area3 = 2,
    n12 = 1,
    n23 = 1,
    n13 = 2,
    n123 = 1,
    category = c('A', 'B', 'C'),
    fill = c('red', 'blue', 'green'),
    cat.col = c('red', 'blue', 'green'),
    cex = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
    cat.cex = c(1,2,3),
    euler = TRUE,
    scaled = FALSE
    )))

#023
venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 2,
    area2 = 1,
    area3 = 2,
    n12 = 1,
    n23 = 1,
    n13 = 2,
    n123 = 1,
    category = c('A', 'B', 'C'),
    fill = c('red', 'blue', 'green'),
    cat.col = c('red', 'blue', 'green'),
    cex = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
    cat.cex = c(1,2,3),
    euler = TRUE,
    scaled = FALSE
    )))

#032
venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 2,
    area2 = 1,
    area3 = 1,
    n12 = 1,
    n23 = 1,
    n13 = 1,
    n123 = 1,
    category = c('A', 'B', 'C'),
    fill = c('red', 'blue', 'green'),
    cat.col = c('red', 'blue', 'green'),
    cex = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
    cat.cex = c(1,2,3),
    euler = TRUE,
    scaled = FALSE
    )))

#033
venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 1,
    area2 = 1,
    area3 = 1,
    n12 = 1,
    n23 = 1,
    n13 = 1,
    n123 = 1,
    category = c('A', 'B', 'C'),
    fill = c('red', 'blue', 'green'),
    cat.col = c('red', 'blue', 'green'),
    cex = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
    cat.cex = c(1,2,3),
    euler = TRUE,
    scaled = FALSE
    )))

#100
venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 3,
    area2 = 3,
    area3 = 3,
    n12 = 1,
    n23 = 1,
    n13 = 1,
    n123 = 0,
    category = c('A', 'B', 'C'),
    fill = c('red', 'blue', 'green'),
    cat.col = c('red', 'blue', 'green'),
    cex = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
    cat.cex = c(1,2,3),
    euler = TRUE,
    scaled = FALSE
    )))

#110
venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 2,
    area2 = 2,
    area3 = 3,
    n12 = 0,
    n23 = 1,
    n13 = 1,
    n123 = 0,
    category = c('A', 'B', 'C'),
    fill = c('red', 'blue', 'green'),
    cat.col = c('red', 'blue', 'green'),
    cex = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
    cat.cex = c(1,2,3),
    euler = TRUE,
    scaled = FALSE
    )))

#111A
venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 2,
    area2 = 1,
    area3 = 3,
    n12 = 0,
    n23 = 1,
    n13 = 1,
    n123 = 0,
    category = c('A', 'B', 'C'),
    fill = c('red', 'blue', 'green'),
    cat.col = c('red', 'blue', 'green'),
    cex = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
    cat.cex = c(1,2,3),
    euler = TRUE,
    scaled = FALSE
    )))

#112AA
venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 1,
    area2 = 1,
    area3 = 3,
    n12 = 0,
    n23 = 1,
    n13 = 1,
    n123 = 0,
    category = c('A', 'B', 'C'),
    fill = c('red', 'blue', 'green'),
    cat.col = c('red', 'blue', 'green'),
    cex = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
    cat.cex = c(1,2,3),
    euler = TRUE,
    scaled = FALSE
    )))

#120
venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 2,
    area2 = 1,
    area3 = 2,
    n12 = 0,
    n23 = 0,
    n13 = 1,
    n123 = 0,
    category = c('A', 'B', 'C'),
    fill = c('red', 'blue', 'green'),
    cat.col = c('red', 'blue', 'green'),
    cex = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
    cat.cex = c(1,2,3),
    euler = TRUE,
    scaled = FALSE
    )))

#121AO
venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 1,
    area2 = 1,
    area3 = 2,
    n12 = 0,
    n23 = 0,
    n13 = 1,
    n123 = 0,
    category = c('A', 'B', 'C'),
    fill = c('red', 'blue', 'green'),
    cat.col = c('red', 'blue', 'green'),
    cex = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
    cat.cex = c(1,2,3),
    euler = TRUE,
    scaled = FALSE
    )))


#122AAOO
venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 1,
    area2 = 1,
    area3 = 1,
    n12 = 0,
    n23 = 0,
    n13 = 1,
    n123 = 0,
    category = c('A', 'B', 'C'),
    fill = c('red', 'blue', 'green'),
    cat.col = c('red', 'blue', 'green'),
    cex = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
    cat.cex = c(1,2,3),
    euler = TRUE,
    scaled = FALSE
    )))

#130
venn.test <- c(venn.test,list(draw.triple.venn(
    area1 = 1,
    area2 = 1,
    area3 = 1,
    n12 = 0,
    n23 = 0,
    n13 = 0,
    n123 = 0,
    category = c('A', 'B', 'C'),
    fill = c('red', 'blue', 'green'),
    cat.col = c('red', 'blue', 'green'),
    cex = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
    cat.cex = c(1,2,3),
    euler = TRUE,
    scaled = FALSE
    )))

testNames <- c('default','colour-default','001','010','011A','011O','012AA','021AA','022AAAO','022AAOO','023','032','033','100','110','111A','112AA','120','121AO','122AAOO','130');

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
for(i in 1:length(venn.plot)){
	test_that(paste("Case",testNames[i],"of three categories"),
	{
		for(j in 1:length(venn.plot[[i]])){
			expect_that(venn.test[[i]][[j]],is_identical_without_name(venn.plot[[i]][[j]]));
		}
	})
}

#Reaches here only if error is not thrown beforehand
print("Three category tests complete. No discrepancies found");
