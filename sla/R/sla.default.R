sla.default <-
function (facxy, ...) 
{
## facxy contains 3 columns: 
	## (1) factor variable for groups,
	## (2) x, the covariate, and 
	## (3) y, the criterion variable.
	.call <- match.call()
	## Set up truth vectors to ensure proper input
	tv1 <- dim(facxy)[2] == 3
	tv2 <- is.factor(facxy[,1])
	tv3 <- is.numeric(facxy[,2])
	tv4 <- is.numeric(facxy[,3])
	tv.string <- c(tv1, tv2, tv3, tv4)
	tv.statement <- c('\n', 'Input data frame is incorrect.', '\n', 'It should contain 3 columns:', '\n', 'Column 1 ~ factor variable, Group', '\n', 'Column 2 ~ numeric variable, x', '\n', 'Column 3 ~ numeric variable, y')
	if(any(tv.string == FALSE)) stop(tv.statement) 
	summary.of.input.df <- summary(facxy)
	group <- facxy[,1] 
	x <- facxy[,2]
	y <- facxy[,3]
	new.df <- data.frame(group, x, y)
	## get lengths of each component in group
	n1 <- with(new.df, table(group))[1]
	n2 <- with(new.df, table(group))[2]
	## build design matrix for all variables, Model A, full model
	## Mod.A fits 4 parameters - intercepts and slopes for 2 groups
	i1 = c(rep(1, n1), rep(0, n2))
	i2 = c(rep(0, n1), rep(1, n2))
	i1x = i1*x
	i2x = i2*x
	X.A = cbind(i1, i1x, i2, i2x)
	## fit Model A - the full model with four parameters
	## Mod.A fits 4 parameters - intercepts and slopes for 2 groups
	Mod.A <- lm(y ~ X.A - 1) ## lm fit for full model
	## build design matrix for all variables, Model B, reduced model
 	X.B = cbind(1, x)
 	## Mod.B fits 2 parameters - an intercept and slope for all data
 	Mod.B <- lm(y ~ X.B - 1) ## linear model fit for reduced model
 	## build design matrix for Model C;  
 	## Mod.C fits 3 parameters - two intercepts and a common slope
 	X.C <- cbind(i1, i2, x)
 	Mod.C <- lm(y ~ X.C - 1) ## lm fit, 2 ints and 1 common slope
 	## build design matrix for Model D; 1 common int and 2 slopes
 	X.D <- cbind(1, i1x, i2x)
 	Mod.D <- lm(y ~ X.D - 1) ## lim fit, 1 int and 2 slopes
 	##	
 	## set up Top Table
 	##
 	TopTable.col.1 <- c('Mod A: Ind Ints, Ind Slopes', 'Mod B: Com Int,  Com Slope', 'Mod C: Ind Ints, Com Slope', 'Mod D: Com Int, Ind Slopes')
 	TopTable.col.2 <- n.parms <- c(4,2,3,3)
 	TopTable.col.3 <- c(anova(Mod.A)[2,1], anova(Mod.B)[2,1], anova(Mod.C)[2,1], anova(Mod.D)[2,1])
 	TopTable.col.4 <- c(anova(Mod.A)[2,2], anova(Mod.B)[2,2], anova(Mod.C)[2,2], anova(Mod.D)[2,2])
	TopTable.col.5 <- c(anova(Mod.A)[2,3], anova(Mod.B)[2,3], anova(Mod.C)[2,3], anova(Mod.D)[2,3])
	TopTable.df <- data.frame(TopTable.col.1, TopTable.col.2, TopTable.col.3, TopTable.col.4, TopTable.col.5)
	names(TopTable.df) = c('Description of Fit', 'Np', 'Resid df', 'Resid SS', 'Resid MS')	
	##
	## set up pretty TopTable
	##
	TopTable.Pretty.col.1 <- 
	c('Mod A: Ind I,Ind S', 
	  'Mod B: Com I,Com S', 
	  'Mod C: Ind I,Com S', 
	  'Mod D: Com I,Ind S')
	print.TopTable.col.2 <- round(TopTable.col.2, 0)
	print.TopTable.col.3 <- round(TopTable.col.3, 0)
	print.TopTable.col.4 <- round(TopTable.col.4, 2)
	print.TopTable.col.5 <- round(TopTable.col.5, 2)
	TopTable.Pretty.df <- data.frame(TopTable.Pretty.col.1, print.TopTable.col.2, print.TopTable.col.3, print.TopTable.col.4, print.TopTable.col.5)
 	names(TopTable.Pretty.df) = c('Description of Fit', 'Np', 'Res df', 'Res SS', 'Res MS')
 	##		
 	## set up Description of Test table
 	##
 	DoT.col.1 = c('Ho: Equiv D.Sets', 'Ho: Equiv Slopes', 'Ho: Equiv Inters')
 	DoT.col.3 = c(anova(Mod.B, Mod.A)[2,3], anova(Mod.C, Mod.A)[2,3], anova(Mod.D, Mod.A)[2,3])
 	DoT.col.4 = c(anova(Mod.B, Mod.A)[2,4], anova(Mod.C, Mod.A)[2,4], anova(Mod.D, Mod.A)[2,4])
 	DoT.col.5 = c(anova(Mod.B, Mod.A)[2,5], anova(Mod.C, Mod.A)[2,5], anova(Mod.D, Mod.A)[2,5])
 	DoT.col.6 = c(anova(Mod.B, Mod.A)[2,6], anova(Mod.C, Mod.A)[2,6], anova(Mod.D, Mod.A)[2,6])
 	DoT.df <- data.frame(DoT.col.1, DoT.col.3, DoT.col.4, DoT.col.5, DoT.col.6)
 	names(DoT.df) = c('Test', 'df', 'SS', 'F Stat', 'prob')
 	## 
 	print.df <- round(DoT.col.3, 0)
 	print.SS <- round(DoT.col.4, 2)
 	print.F <- round(DoT.col.5, 2)
 	print.prob <- round(DoT.col.6, 4)
 	print.data.frame <- data.frame(DoT.col.1, print.df, print.SS, print.F, print.prob)
 	## Print Pretty Fit Table
 	title.string.Pretty.Fit.table <- '    Description of Fits for 4 ANCOVA Models'
 	names(print.data.frame) = c('Test', 'df', 'SS', 'F Stat', 'prob')
	title.string <- '   ANCOVA Tests: Two Groups/Straight Line Fits'
	## 	
	## 		
	slaObj <- list(Call = .call, INPUT.df = facxy, 'Summary of Input Data Frame' = summary.of.input.df, Mod.A = Mod.A, Mod.B = Mod.B, Mod.C = Mod.C, Mod.D = Mod.D, Fit.Table = TopTable.df, 
Test.Table = DoT.df, Fit.Table.Pretty = TopTable.Pretty.df, Test.Table.Pretty = print.data.frame)
	## IMPORTANT: assign class to slaObj
	class(slaObj) <- 'sla'
	slaObj   
}
