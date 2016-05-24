# tests for formation and dissolution functions
require(tsna)
require(testthat)
require(networkDynamicData)

# ------ tEdgeDissolution -------
test<-network.initialize(4,directed = TRUE,loops = FALSE)
add.edges.active(test,1,2,onset = 0,terminus=3)
# check for ts class
expect_true(is.ts(tEdgeDissolution(test,start=0,end=3)))
# correct count
expect_equal(as.numeric(tEdgeDissolution(test,start=0,end=3)),c(0,0,0,1))
expect_equal(as.numeric(tEdgeDissolution(test,start=0,end=3,time.interval=0.5)),c(0,0,0,0,0,0, 1))

# check with fractions
add.edges.active(test,2,1,onset = 0,terminus=2)  
expect_equal(as.numeric(tEdgeDissolution(test,result.type = 'fraction')),c(0,0,0.5,1)) 
             
             
expect_warning(tEdgeDissolution(network.initialize(0)),regexp = 'network does not appear to have any time range information')





# -------- tEdgeFormation ---------
test<-network.initialize(4,directed = TRUE,loops = FALSE)
add.edges.active(test,1,2,onset = 0,terminus=3)
add.edges.active(test,2,3,onset = 1,terminus=3)

# test class
expect_true(is.ts(tEdgeFormation(test)))

# test count
expect_equal(as.numeric(tEdgeFormation(test)),c(1,1,0,0))

# test fraction
expect_equal(as.numeric(tEdgeFormation(test,result.type = 'fraction')),c(1/12, 1/11,0,0))


# test various network types
test%n%'directed'<-FALSE
expect_equal(as.numeric(tEdgeFormation(test,result.type = 'fraction')),c(1/6, 1/5,0,0))

test%n%'loops'<-TRUE
expect_equal(as.numeric(tEdgeFormation(test,result.type = 'fraction')),c(1/10, 1/9,0,0))

test%n%'directed'<-TRUE
expect_equal(as.numeric(tEdgeFormation(test,result.type = 'fraction')),c(1/16, 1/15,0,0))

# check bipartite
test<-network.initialize(5,directed = TRUE,loops = FALSE,bipartite=2)
add.edges.active(test,1,3,onset=0,terminus=2)
expect_equal(as.numeric(tEdgeFormation(test,result.type='fraction')),c(1/6,0,0))

# check multiplex 
expect_error(tEdgeFormation(network.initialize(0,multiple = TRUE),result.type='fraction'),regexp = 'can not compute possible number of free dyads for multiplex networks')
# 

# check censoring toggle
data(concurrencyComparisonNets)
expect_equal(as.numeric(tEdgeFormation(base,start=1,end=4,include.censored = FALSE)),c(0,0,24,18))
expect_equal(as.numeric(tEdgeFormation(base,start=1,end=4,include.censored = TRUE)),c(0,371,24,18))

expect_equal(as.numeric(tEdgeDissolution(base,start=100,end=103,include.censored = TRUE)),c(23,20,364,0))
expect_equal(as.numeric(tEdgeDissolution(base,start=100,end=103,include.censored = FALSE)),c(23,20,18,0))
