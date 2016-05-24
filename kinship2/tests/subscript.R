#
# Test out subscripting
#
require(kinship2)
data(minnbreast)

minnped <- with(minnbreast, pedigree(id, fatherid, motherid, sex,
                                     affected=cancer, famid=famid))
ped8 <- minnped['8']  # a modest sized family

# Subjects 150, 152, 154, 158 are children, and 143, 162, 149 are 
#  parents and a child
droplist <- c(150, 152, 154, 158, 143, 162, 149)

keep1 <- !(ped8$id %in% droplist)  #logical
keep2 <- which(keep1)              #numeric
keep3 <- as.character(ped8$id[keep1]) #character
keep4 <- factor(keep3)

test1 <- ped8[keep1]
test2 <- ped8[keep2]
test3 <- ped8[keep3]
test4 <- ped8[keep4]
all.equal(test1, test2)
all.equal(test1, test3)
all.equal(test1, test4)
