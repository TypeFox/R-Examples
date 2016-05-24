
## @knitr load_turner, message=FALSE
# load package turner
library(turner)


## @knitr data_matrix
# create a matrix
set.seed = 21
some_data = round(matrix(rnorm(90), 10, 9), 3)
rownames(some_data) = 1:10
colnames(some_data) = paste("X", 1:9, sep='')

# take a peek
head(some_data, n=5)


## @knitr blocks
# list of blocks
blocks = list(B1 = 1:3, B2 = 4:5, B3 = 6:9)
blocks


## @knitr indexify
# get indices of blocks
indices = indexify(blocks)
indices


## @knitr list_to_dummy
# get dummy matrix based on blocks
dummy = list_to_dummy(blocks)
dummy


## @knitr from_to
# get starting and ending positions
start_end = from_to(blocks)
start_end

# vectors from and to
from = start_end$from
to = start_end$to


## @knitr extract_first_block
# extract first block
some_data[,from[1]:to[1]]


## @knitr get_first_block
# get first block
some_data[,blocks[[1]]]


## @knitr str_list
# string list
str_list = list(c("a","b","c"), c("d", "e"), c("f","g","h","i"))


## @knitr failed_extraction, eval=FALSE
## # failed attempt
## some_data[,str_list[[1]]]


## @knitr success_extraction
# start-end position for 'str_list'
fromto_aux = from_to(str_list)
from1 = fromto_aux$from
to1 = fromto_aux$to

# successful attempt
some_data[,from1[1]:to1[1]]


## @knitr lenghts
# say you have some list
some_list = list(1:3, 4:5, 6:9)

# length of each vector (vector output)
lengths(some_list, out='vector')

# length of each vector (list output)
lengths(some_list, out='list')


## @knitr length
# compared to 'length()'
length(some_list)


## @knitr funlist
# sum of all elements in 'some_list'
funlist(some_list, sum)

# maixmum of all elements in 'some_list'
funlist(some_list, max)

# product of all elements in 'some_list'
funlist(some_list, prod)

# mean value of all elements in 'some_list'
funlist(some_list, mean)


## @knitr listsize
# number of elements in 'some_list'
listsize(some_list)


## @knitr listify
# vector of indices
number_elements = c(3, 1, 5)

# list of index vectors based on 'number_elements'
listify(number_elements)


