require(R.oo)

# CONSTRUCTOR

setConstructorS3("CellularAutomaton", function(n = 30, fun = NULL, k = 2, r = 1, t = 1, totalistic = 0, seed = 1, bg = 0)
{
    this = extend(Object(), "CellularAutomaton",
        .grid = NULL,                              # This field is the grid of cells.
        .ruleCode = NULL,                          # The .rule() method looks up the next state of a cell in this table.
        .rule = NULL,                              # This method determines the next state of a cell.
        .n = NULL,                                 # This is the rule number.
        .k = NULL,                                 # This is the number of colors.
        .r = NULL,                                 # This is the radius of a neighborhood.
        .t = NULL,                                 # This is the number of steps.
        .totalistic = NULL                         # 0 => general; 1 => totalistic
    )

    this$.validateOptions(n, k, r, t, totalistic)  # Make sure the options are valid.

    if (missing(fun))               # If the user doesn't supply a rule function,
    {
        if (this$.n >= 0)           # then he/she must supply a rule number,
            this$.generateRule()    # in which case we generate .ruleCode and .rule().
        else
            throw("You must specify a rule number (n) or provide your own rule (fun).") # Otherwise, complain.
    }
    else
    {
        this$.n = -1
        this$.rule = fun  # If the user supplied a rule function, point to it with the .rule field.
    }

    this$.initializeGrid(seed, bg)  # Initialize the grid.

    this$.run()  # Run for t generations.

    this
})

# HELPER METHOD   .validateOptions()

setMethodS3(".validateOptions", "CellularAutomaton", function(this, n, k, r, t, totalistic,...)
{
    n = as.integer(n)
    if (is.na(n) | length(n) == 0)  # n should be a nonnegative integer.
        n = -1
    this$.n = n

    k = as.integer(k)
    if (is.na(k) | length(k) == 0 | k < 2)  # k should be an integer >= 2.
        k = 2
    this$.k = k

    r = as.integer(r)
    if (is.na(r) | length(r) == 0 | r < 1)  # r should be an integer >= 1.
        r = 1
    this$.r = r

    t = as.integer(t)
    if (is.na(t) | length(t) == 0 | t < 1)  # t should be an integer >= 1.
        t = 1
    this$.t = t

    # totalistic must equal 0 or 1.

    totalistic = as.integer(totalistic)
    if (is.na(totalistic) | length(totalistic) == 0 | totalistic < 0 | totalistic > 1)
        totalistic = 0
    this$.totalistic = totalistic

}, protected = T)

# HELPER METHOD   .generateRule()

setMethodS3(".generateRule", "CellularAutomaton", function(this,...)
{
    if (this$.totalistic < 2)
    {
        n = this$.n
        k = this$.k
        ruleCode = rep(0, k^(2 * this$.r + 1))
        j = 1
        while (n > 0)
        {
            ruleCode[j] = n %% k
            j = j + 1
            n = n %/% k
        }
        this$.ruleCode = ruleCode
        if (this$.totalistic == 0)
            this$.rule = function(neighborhood)
            {
                index = 1
                for (j in 0:(2 * this$.r))
                    index = index + this$.k^j * neighborhood[length(neighborhood) - j]
                this$.ruleCode[index]
            }
        else
            this$.rule = function(neighborhood)
            {
                this$.ruleCode[sum(neighborhood) + 1]
            }
    }

}, protected = T)

# HELPER METHOD   .initializeGrid()

setMethodS3(".initializeGrid", "CellularAutomaton", function(this, seed, bg,...)
{
    if (bg == -1)
        time0 = seed
    else
    {
        time0 = rep(bg, (2 * this$.r * this$.t + 1) %/% length(bg))
        length(time0) = 2 * this$.r * this$.t + 1
        time0[is.na(time0)] = 0
        pos = length(time0) %/% 2 + 1 - length(seed) %/% 2
        for (j in 1:length(seed))
        {
            time0[pos] = seed[j]
            pos = pos + 1
        }
    }
    this$.grid = matrix(time0, nrow = 1, ncol = length(time0), byrow = TRUE)

}, protected = T)

# HELPER METHOD   .buildNeighborhood()

setMethodS3(".buildNeighborhood", "CellularAutomaton", function(this, current, j,...)
{
    neighborhood = numeric(2 * this$.r + 1)
    pos = j - this$.r
    if (pos <= 0)
        pos = length(current) + pos
    for (i in 1:length(neighborhood))
    {
        neighborhood[i] = current[pos]
        pos = (pos + 1) %% length(current)
        if (pos == 0)
            pos = length(current)
    }
    neighborhood

}, protected = T)

# HELPER METHOD   .run()

setMethodS3(".run", "CellularAutomaton", function(this,...)
{
    for (i in 1:this$.t)
    {
        current = this$.grid[i, ]
        numCells = ncol(this$.grid)
        nxt = numeric(numCells)
        for (j in 1:numCells)
        {
            neighborhood = this$.buildNeighborhood(current, j)
            nxt[j] = this$.rule(neighborhood)
        }
        this$.grid = rbind(this$.grid, nxt)
    }
    rownames(this$.grid) = paste("t", 0:this$.t, sep = "")

}, protected = T)

# PUBLIC METHOD   plot()

# This method plots the lattice using the image() function. The matrix is flipped
# horizontally and transposed before image() is called so that the plot shows the
# automaton's evolution from top to bottom. A list of colors may be supplied. If
# no colors are supplied, the plotting colors default to 0, 1,..., (k - 1).

setMethodS3("plot", "CellularAutomaton", function(x,..., col)
{
    if (missing(col))
        col = 0:(x$.k - 1)
    if (x$.n < 0)
        main = "user-defined rule"
    else
        main = paste("Rule", x$.n)
    main = c(main, paste("(k = ", x$.k, ", r = ", x$.r, ", ", x$.t, " steps)", sep = ""))
    image(t(x$.grid[nrow(x$.grid):1,]), main = main, col = col, axes = FALSE, asp = nrow(x$.grid) / ncol(x$.grid))
})

# PUBLIC METHOD   getRuleNumber()

setMethodS3("getRuleNumber", "CellularAutomaton", function(this,...)
{
    this$.n
})

# PUBLIC METHOD   getNumberOfColors()

setMethodS3("getNumberOfColors", "CellularAutomaton", function(this,...)
{
    this$.k
})

# PUBLIC METHOD   getRadius()

setMethodS3("getRadius", "CellularAutomaton", function(this,...)
{
    this$.r
})

# PUBLIC METHOD   getSteps()

setMethodS3("getSteps", "CellularAutomaton", function(this,...)
{
    this$.t
})

# PUBLIC METHOD   getTotalistic()

setMethodS3("getTotalistic", "CellularAutomaton", function(this,...)
{
    this$.totalistic
})

# PUBLIC METHOD   getLattice()

setMethodS3("getLattice", "CellularAutomaton", function(this,...)
{
    this$.grid
})
