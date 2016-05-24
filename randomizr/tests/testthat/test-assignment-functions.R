library(testthat)
library(randomizr)

context("Complete Random Assignments")

expect_error(complete_ra(100, num_arms = 1))
expect_error(complete_ra(100, condition_names = c("Treatment")))
expect_error(complete_ra(100, m_each = c(100)))
expect_error(complete_ra(100, prob_each = c(.2)))
expect_error(complete_ra(100, m=50, condition_names = c("Control", "Treatment")))
expect_error(complete_ra(100, m_each = c(30, 70), prob_each = c(.3, .7)))
expect_error(complete_ra(2, 2))
expect_error(complete_ra(100, m_each=c(20, 20, 20)))
expect_error(complete_ra(100, m_each=c(20, 20, 60), condition_names=c(1,2)))
expect_error(complete_ra(100, prob_each = c(.2, .7)))

Z <- complete_ra(N=100)
expect_equal(sum(Z), 50)

Z <- complete_ra(N=101, prob = .34)
table(Z)

Z <- complete_ra(N=100, m=50)
expect_equal(sum(Z), 50)

Z <- complete_ra(N=100, m_each = c(30, 70), 
                 condition_names = c("control", "treatment"))

expect_equivalent(as.numeric(table(Z)), c(30, 70))

# Multi-arm Designs
Z <- complete_ra(N=100, num_arms=3)

expect_true(all(table(Z) %in% c(33, 34)))

Z <- complete_ra(N=100, m_each=c(30, 30, 40))

expect_equivalent(as.numeric(table(Z)), c(30, 30, 40))

Z <- complete_ra(N=100, m_each=c(30, 30, 40), 
                 condition_names=c("control", "placebo", "treatment"))
expect_equivalent(as.numeric(table(Z)), c(30, 30, 40))

Z <- complete_ra(N=100, condition_names=c("control", "placebo", "treatment"))

expect_true(all(table(Z) %in% c(33, 34)))

# Special Cases

expect_less_than(sum(replicate(1000, complete_ra(1))), 600)
expect_more_than(sum(replicate(1000, complete_ra(1))), 400)

context("Simple Random Assignments")

# Two Group Designs
Z <- simple_ra(N=100)
table(Z)

Z <- simple_ra(N=100, prob=0.5)
table(Z)

Z <- simple_ra(N=100, prob_each = c(0.3, 0.7), 
               condition_names = c("control", "treatment"))
table(Z)

# Multi-arm Designs
Z <- simple_ra(N=100, num_arms=3)
table(Z)

Z <- simple_ra(N=100, prob_each=c(0.3, 0.3, 0.4))
table(Z)

Z <- simple_ra(N=100, prob_each=c(0.3, 0.3, 0.4), 
               condition_names=c("control", "placebo", "treatment"))
table(Z)

Z <- simple_ra(N=100, condition_names=c("control", "placebo", "treatment"))
table(Z)

context("Block Random Assignments")

block_var <- rep(c("A", "B","C"), times=c(50, 100, 200))
Z <- block_ra(block_var=block_var)
table(block_var, Z)

block_m <- rbind(c(25, 25),
                 c(50, 50),
                 c(100, 100))

Z <- block_ra(block_var=block_var, block_m=block_m)
table(block_var, Z)

block_m <- rbind(c(10, 40),
                 c(30, 70),
                 c(50, 150))

Z <- block_ra(block_var=block_var, block_m=block_m, 
              condition_names=c("control", "treatment"))
table(block_var, Z)

# Multi-arm Designs
Z <- block_ra(block_var=block_var, num_arms=3)
table(block_var, Z)

block_m <- rbind(c(10, 20, 20),
                 c(30, 50, 20),
                 c(50, 75, 75))
Z <- block_ra(block_var=block_var, block_m=block_m )
table(block_var, Z)

Z <- block_ra(block_var=block_var, block_m=block_m, 
              condition_names=c("control", "placebo", "treatment"))
table(block_var, Z)

Z <- block_ra(block_var=block_var, prob_each=c(.1, .1, .8))
table(block_var, Z)

context("Cluster Random Assignments")

# Two Group Designs
clust_var <- rep(letters, times=1:26)
Z <- cluster_ra(clust_var=clust_var)
table(Z, clust_var)

Z <- cluster_ra(clust_var=clust_var, m=13)
table(Z, clust_var)
Z <- cluster_ra(clust_var=clust_var, m_each = c(10, 16), 
                condition_names = c("control", "treatment"))

table(Z, clust_var)

# Multi-arm Designs
Z <- cluster_ra(clust_var=clust_var, num_arms=3)
table(Z, clust_var)

Z <- cluster_ra(clust_var=clust_var, m_each=c(7, 7, 12))
table(Z, clust_var)

Z <- cluster_ra(clust_var=clust_var, m_each=c(7, 7, 12), 
                condition_names=c("control", "placebo", "treatment"))
table(Z, clust_var)

Z <- cluster_ra(clust_var=clust_var, 
                condition_names=c("control", "placebo", "treatment"))
table(Z, clust_var)

Z <- cluster_ra(clust_var=clust_var, prob_each = c(.1, .2, .7))
table(Z, clust_var)


# blocked and cluster -----------------------------------------------------

clust_var <- rep(letters, times=1:26)
block_var <- rep(NA, length(clust_var))
block_var[clust_var %in% letters[1:5]] <- "block_1"
block_var[clust_var %in% letters[6:10]] <- "block_2"
block_var[clust_var %in% letters[11:15]] <- "block_3"
block_var[clust_var %in% letters[16:20]] <- "block_4"
block_var[clust_var %in% letters[21:26]] <- "block_5"

table(block_var, clust_var)

Z <- block_and_cluster_ra(clust_var = clust_var, block_var = block_var)

table(Z, clust_var)
table(Z, block_var)

Z <- block_and_cluster_ra(clust_var = clust_var, block_var = block_var, num_arms = 3)

table(Z, clust_var)
table(Z, block_var)

Z <- block_and_cluster_ra(clust_var = clust_var, block_var = block_var, prob_each = c(.2, .5, .3))


block_m <- rbind(c(2, 3),
                 c(1, 4),
                 c(3, 2),
                 c(2, 3),
                 c(5, 1))

Z <- block_and_cluster_ra(clust_var = clust_var, block_var = block_var, block_m = block_m)

table(Z, clust_var)
table(Z, block_var)


