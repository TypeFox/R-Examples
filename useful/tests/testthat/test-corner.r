context('Ensure corner returns the correct part of the data.')

topL0 <- WhichCorner('topleft')
bottomL0 <- WhichCorner('bottomleft')
topR0 <- WhichCorner('topright')
bottomR0 <- WhichCorner('bottomright')

topLR6 <- WhichCorner('topleft', r=6)
bottomLR6 <- WhichCorner('bottomleft', r=6)
topRR6 <- WhichCorner('topright', r=6)
bottomRR6 <- WhichCorner('bottomright', r=6)

topLC7 <- WhichCorner('topleft', c=7)
bottomLC7 <- WhichCorner('bottomleft', c=7)
topRC7 <- WhichCorner('topright', c=7)
bottomRC7 <- WhichCorner('bottomright', c=7)

topLR8C3 <- WhichCorner('topleft', r=8, c=3)
bottomLR8C3 <- WhichCorner('bottomleft', r=8, c=3)
topRR8C3 <- WhichCorner('topright', r=8, c=3)
bottomRR8C3 <- WhichCorner('bottomright', r=8, c=3)

test_that('WhichCorner returns expression class', {
    expect_is(topL0, 'expression')
    expect_is(bottomL0, 'expression')
    expect_is(topR0, 'expression')
    expect_is(bottomR0, 'expression')
    
    expect_is(topLR6, 'expression')
    expect_is(bottomLR6, 'expression')
    expect_is(topRR6, 'expression')
    expect_is(bottomRR6, 'expression')
    
    expect_is(topLC7, 'expression')
    expect_is(bottomLC7, 'expression')
    expect_is(topRC7, 'expression')
    expect_is(bottomRC7, 'expression')
    
    expect_is(topLR8C3, 'expression')
    expect_is(bottomLR8C3, 'expression')
    expect_is(topRR8C3, 'expression')
    expect_is(bottomRR8C3, 'expression')
})

test_that('WhichCorner gets the right number of rows and columns', {
    expect_equal(as.character(topL0), 'list(rows = 1:5, cols = 1:5)')
    expect_equal(as.character(bottomL0), 'list(rows = (nrow(x) - 5 + 1):nrow(x), cols = 1:5)')
    expect_equal(as.character(topR0), 'list(rows = 1:5, cols = (ncol(x) - 5 + 1):ncol(x))')
    expect_equal(as.character(bottomR0), 'list(rows = (nrow(x) - 5 + 1):nrow(x), cols = (ncol(x) - 5 + 1):ncol(x))')
    
    expect_equal(as.character(topLR6), 'list(rows = 1:6, cols = 1:5)')
    expect_equal(as.character(bottomLR6), 'list(rows = (nrow(x) - 6 + 1):nrow(x), cols = 1:5)')
    expect_equal(as.character(topRR6), 'list(rows = 1:6, cols = (ncol(x) - 5 + 1):ncol(x))')
    expect_equal(as.character(bottomRR6), 'list(rows = (nrow(x) - 6 + 1):nrow(x), cols = (ncol(x) - 5 + 1):ncol(x))')
    
    expect_equal(as.character(topLC7), 'list(rows = 1:5, cols = 1:7)')
    expect_equal(as.character(bottomLC7), 'list(rows = (nrow(x) - 5 + 1):nrow(x), cols = 1:7)')
    expect_equal(as.character(topRC7), 'list(rows = 1:5, cols = (ncol(x) - 7 + 1):ncol(x))')
    expect_equal(as.character(bottomRC7), 'list(rows = (nrow(x) - 5 + 1):nrow(x), cols = (ncol(x) - 7 + 1):ncol(x))')
    
    expect_equal(as.character(topLR8C3), 'list(rows = 1:8, cols = 1:3)')
    expect_equal(as.character(bottomLR8C3), 'list(rows = (nrow(x) - 8 + 1):nrow(x), cols = 1:3)')
    expect_equal(as.character(topRR8C3), 'list(rows = 1:8, cols = (ncol(x) - 3 + 1):ncol(x))')
    expect_equal(as.character(bottomRR8C3), 'list(rows = (nrow(x) - 8 + 1):nrow(x), cols = (ncol(x) - 3 + 1):ncol(x))')
})

data('diamonds', package='ggplot2')

corn00 <- corner(diamonds)
corn40 <- corner(diamonds, r=4)
corn03 <- corner(diamonds, c=3)
corn37 <- corner(diamonds, r=3, c=7)

cornTL00 <- corner(diamonds, corner='topleft')
cornTL40 <- corner(diamonds, r=4, corner='topleft')
cornTL03 <- corner(diamonds, c=3, corner='topleft')
cornTL37 <- corner(diamonds, r=3, c=7, corner='topleft')

cornTR00 <- corner(diamonds, corner='topright')
cornTR40 <- corner(diamonds, corner='topright', r=4)
cornTR03 <- corner(diamonds, corner='topright', c=3)
cornTR37 <- corner(diamonds, corner='topright', r=3, c=7)

cornBL00 <- corner(diamonds, corner='bottomleft')
cornBL40 <- corner(diamonds, corner='bottomleft', r=4)
cornBL03 <- corner(diamonds, corner='bottomleft', c=3)
cornBL37 <- corner(diamonds, corner='bottomleft', r=3, c=7)

cornBR00 <- corner(diamonds, corner='bottomright')
cornBR40 <- corner(diamonds, corner='bottomright', r=4)
cornBR03 <- corner(diamonds, corner='bottomright', c=3)
cornBR37 <- corner(diamonds, corner='bottomright', r=3, c=7)

test_that('Corner returns proper classes', {
    expect_is(corn00, 'data.frame')
    expect_is(corn40, 'data.frame')
    expect_is(corn03, 'data.frame')
    expect_is(corn37, 'data.frame')
    
    expect_is(cornTL00, 'data.frame')
    expect_is(cornTL40, 'data.frame')
    expect_is(cornTL03, 'data.frame')
    expect_is(cornTL37, 'data.frame')
    
    expect_is(cornTR00, 'data.frame')
    expect_is(cornTR40, 'data.frame')
    expect_is(cornTR03, 'data.frame')
    expect_is(cornTR37, 'data.frame')
    
    expect_is(cornBL00, 'data.frame')
    expect_is(cornBL40, 'data.frame')
    expect_is(cornBL03, 'data.frame')
    expect_is(cornBL37, 'data.frame')
    
    expect_is(cornBR00, 'data.frame')
    expect_is(cornBR40, 'data.frame')
    expect_is(cornBR03, 'data.frame')
    expect_is(cornBR37, 'data.frame')
})

test_that('Corner returns proper sizes', {
    expect_equal(dim(corn00), c(5, 5))
    expect_equal(dim(corn40), c(4, 5))
    expect_equal(dim(corn03), c(5, 3))
    expect_equal(dim(corn37), c(3, 7))
    
    expect_equal(dim(cornTL00), c(5, 5))
    expect_equal(dim(cornTL40), c(4, 5))
    expect_equal(dim(cornTL03), c(5, 3))
    expect_equal(dim(cornTL37), c(3, 7))
    
    expect_equal(dim(cornTR00), c(5, 5))
    expect_equal(dim(cornTR40), c(4, 5))
    expect_equal(dim(cornTR03), c(5, 3))
    expect_equal(dim(cornTR37), c(3, 7))
    
    expect_equal(dim(cornBL00), c(5, 5))
    expect_equal(dim(cornBL40), c(4, 5))
    expect_equal(dim(cornBL03), c(5, 3))
    expect_equal(dim(cornBL37), c(3, 7))
    
    expect_equal(dim(cornBR00), c(5, 5))
    expect_equal(dim(cornBR40), c(4, 5))
    expect_equal(dim(cornBR03), c(5, 3))
    expect_equal(dim(cornBR37), c(3, 7))
})


topleft00 <- topleft(diamonds)
topleft40 <- topleft(diamonds, r=4)
topleft03 <- topleft(diamonds, c=3)
topleft37 <- topleft(diamonds, r=3, c=7)

topright00 <- topright(diamonds)
topright40 <- topright(diamonds, r=4)
topright03 <- topright(diamonds, c=3)
topright37 <- topright(diamonds, r=3, c=7)

bottomleft00 <- bottomleft(diamonds)
bottomleft40 <- bottomleft(diamonds, r=4)
bottomleft03 <- bottomleft(diamonds, c=3)
bottomleft37 <- bottomleft(diamonds, r=3, c=7)

bottomright00 <- bottomright(diamonds)
bottomright40 <- bottomright(diamonds, r=4)
bottomright03 <- bottomright(diamonds, c=3)
bottomright37 <- bottomright(diamonds, r=3, c=7)

test_that('Top, bottom, left, right functions return proper class', {
    expect_is(topleft00, 'data.frame')
    expect_is(topleft40, 'data.frame')
    expect_is(topleft03, 'data.frame')
    expect_is(topleft37, 'data.frame')
    
    expect_is(topright00, 'data.frame')
    expect_is(topright40, 'data.frame')
    expect_is(topright03, 'data.frame')
    expect_is(topright37, 'data.frame')
    
    expect_is(bottomleft00, 'data.frame')
    expect_is(bottomleft40, 'data.frame')
    expect_is(bottomleft03, 'data.frame')
    expect_is(bottomleft37, 'data.frame')
    
    expect_is(bottomright00, 'data.frame')
    expect_is(bottomright40, 'data.frame')
    expect_is(bottomright03, 'data.frame')
    expect_is(bottomright37, 'data.frame')
})

test_that('Top, bottom, left, right functions return proper sizes', {
    expect_equal(dim(topleft00), c(5, 5))
    expect_equal(dim(topleft40), c(4, 5))
    expect_equal(dim(topleft03), c(5, 3))
    expect_equal(dim(topleft37), c(3, 7))
    
    expect_equal(dim(topright00), c(5, 5))
    expect_equal(dim(topright40), c(4, 5))
    expect_equal(dim(topright03), c(5, 3))
    expect_equal(dim(topright37), c(3, 7))
    
    expect_equal(dim(bottomleft00), c(5, 5))
    expect_equal(dim(bottomleft40), c(4, 5))
    expect_equal(dim(bottomleft03), c(5, 3))
    expect_equal(dim(bottomleft37), c(3, 7))
    
    expect_equal(dim(bottomright00), c(5, 5))
    expect_equal(dim(bottomright40), c(4, 5))
    expect_equal(dim(bottomright03), c(5, 3))
    expect_equal(dim(bottomright37), c(3, 7))
})

left0 <- left(diamonds)
left3 <- left(diamonds, c=3)
right0 <- right(diamonds)
right3 <- right(diamonds, c=3)

test_that('Right and left return proper classes', {
    expect_is(left0, 'data.frame')
    expect_is(left3, 'data.frame')
    expect_is(right0, 'data.frame')
    expect_is(right0, 'data.frame')
})

test_that('Right and left return proper sizes', {
    expect_equal(dim(left0), c(53940, 5))
    expect_equal(dim(left3), c(53940, 3))
    expect_equal(dim(right0), c(53940, 5))
    expect_equal(dim(right3), c(53940, 3))
})

# test_that('', {
#     expect_is(eval(WhichCorner(corner="topleft", 5, 5, "testFrame")), 'list')
#     expect_is(eval(WhichCorner(corner="topright", 5, 5, "testFrame")), 'list')
#     expect_is(eval(WhichCorner(corner="bottomleft", 5, 5, "testFrame")), 'list')
#     expect_is(eval(WhichCorner(corner="bottomright", 5, 5, "testFrame")), 'list')
#     
#     eval(WhichCorner(corner="topleft", 5, 5, "testFrame"))
#     eval(WhichCorner(corner="topright", 5, 5, "testFrame"))
#     eval(WhichCorner(corner="bottomleft", 5, 5, "testFrame"))
#     eval(WhichCorner(corner="bottomright", 5, 5, "testFrame"))
# })