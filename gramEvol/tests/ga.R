# This program performs regression testing for package "gramEvol"

# Author: Farzad Noorian <farzad.noorian@sydney.edu.au>

# This program is free software, distributed under the terms of
# the GNU General Public License version 2.
# Refer to <http://www.gnu.org/licenses/> for full terms.
################################################################################
library("gramEvol")

set.seed(0)

odd <- seq(1, 20, 2)
even <- seq(2, 20, 2)
evalfunc <- function(l) {
    err <- sum(l[even]) - sum(l[odd]);

    stopifnot(!any(duplicated(l)))

    return (err)
}

x <- GeneticAlg.int(genomeLen = 20, codonMin = 0, codonMax = 20,
                allowrepeat = FALSE, terminationCost = -110,
                monitorFunc = NULL, evalFunc = evalfunc)

best.result <- x$best$genome

stopifnot(sort(best.result[odd]) == 11:20)
stopifnot(sort(best.result[even]) == 0:9)

