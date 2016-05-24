#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: sienaeffects.r
# *
# * Description: This module contains utilities for updating an effects object
# *****************************************************************************/
##@includeEffect DataCreate
includeEffects <- function(myeff, ..., include=TRUE, name=myeff$name[1],
                           type="eval", interaction1="", interaction2="",
                           character=FALSE)
{
    if (character)
    {
        dots <- sapply(list(...), function(x)x)
    }
    else
    {
        dots <- substitute(list(...))[-1] ##first entry is the word 'list'
    }
    if (length(dots) == 0)
    {
        stop("This function needs some effect short names.")
    }
    if (!character)
    {
        effectNames <- sapply(dots, function(x)deparse(x))
    }
    else
    {
        effectNames <- dots
    }
    use <- myeff$shortName %in% effectNames &
    myeff$type==type &
    myeff$name==name &
    myeff$interaction1 == interaction1 &
    myeff$interaction2 == interaction2
    myeff[use, "include"] <- include
	if (sum(use) <= 0)
	{
		cat(paste("There is no effect with short name "))
		cat(paste(effectNames,", \n", sep=""))
		cat(paste("and with interaction1 = <",interaction1,">", sep=""))
		cat(paste(" and interaction2 = <",interaction2,">,\n", sep=""))
		cat(paste("for dependent variable",name,".\n"))
	}
	else
	{
		print.data.frame(myeff[use, c("name", "shortName", "type",
			"interaction1", "interaction2", "include")])
	}
    myeff
}
##@includeInteraction DataCreate
includeInteraction <- function(myeff, ...,
                               include=TRUE, name=myeff$name[1],
                        type="eval", interaction1=rep("", 3),
                               interaction2=rep("", 3), character=FALSE,
                               verbose=TRUE)
{
    if (character)
    {
        dots <- sapply(list(...), function(x)x)
    }
    else
    {
        ## check we have 2 or 3 short names
        dots <- substitute(list(...))[-1] ##first entry is the word 'list'
    }
    if (length(dots) == 0)
    {
        stop("need some effect short names")
    }
    if (length(dots) < 2 || length(dots) > 3)
    {
         stop("need exactly two or three effect short names")
    }
    if (!character)
    {
        shortNames <- sapply(dots, function(x)deparse(x))
    }
    else
    {
        shortNames <- dots
    }
    ## find the first underlying effect
    shortName <- shortNames[1]
    interact1 <- interaction1[1]
    interact2 <- interaction2[1]
    use <- myeff$shortName == shortName &
    myeff$type==type &
    myeff$name==name &
    myeff$interaction1 == interact1 &
    myeff$interaction2 == interact2
    if (sum(use) == 0)
    {
        stop("First effect not found")
    }
    if (sum(use) > 1)
    {
        stop("First effect not unique")
    }
    effect1 <- myeff[use, "effectNumber"]
    ## find the second underlying effect
    shortName <- shortNames[2]
    interact1 <- ifelse (length(interaction1) > 1, interaction1[2], "")
    interact2 <- ifelse (length(interaction2) > 1, interaction2[2], "")
    use <- myeff$shortName == shortName &
    myeff$type==type &
    myeff$name==name &
    myeff$interaction1 == interact1 &
    myeff$interaction2 == interact2
    if (sum(use) == 0)
    {
        stop("Second effect not found")
    }
    if (sum(use) > 1)
    {
        stop("Second effect not unique")
    }
    effect2 <- myeff[use, "effectNumber"]
    ## find the third underlying effect, if any
    if (length(shortNames) > 2)
    {
        shortName <- shortNames[3]
        interact1 <- ifelse (length(interaction1) > 2, interaction1[3], "")
        interact2 <- ifelse (length(interaction2) > 2, interaction2[3], "")
        use <- myeff$shortName == shortName &
        myeff$type==type &
        myeff$name==name &
        myeff$interaction1 == interact1 &
        myeff$interaction2 == interact2
        if (sum(use) == 0)
        {
            stop("Third effect not found")
        }
        if (sum(use) > 1)
        {
            stop("Third effect not unique")
        }
        effect3 <- myeff[use, "effectNumber"]
    }
    else
    {
        effect3 <- 0
    }
    ## if want to include, check that we have a spare row
    if (include)
    {
        ints <- myeff[myeff$name == name & myeff$shortName  %in%
                      c("unspInt", "behUnspInt") &
                      (is.na(myeff$effect1) | myeff$effect1 == 0)&
                      myeff$type == type, ]
        if (nrow(ints) == 0)
        {
            baseEffect<- myeff[myeff$name == name, ][1, ]
            if (baseEffect$netType != "behavior")
            {
                tmprow <- createEffects("unspecifiedNetInteraction", name=name,
                                        netType=baseEffect$netType,
                                        groupName=baseEffect$groupName,
                                        group=baseEffect$group)
            }
            else
            {
                tmprow <- createEffects("unspecifiedBehaviorInteraction",
                                        name=name,
                                        netType=baseEffect$netType,
                                        groupName=baseEffect$groupName,
                                        group=baseEffect$group)
            }
            tmprow$include <- TRUE
            tmprow <- tmprow[tmprow$type==type, ]
            tmprow$effectNumber <- max(myeff$effectNumber) + 1
            rownames(tmprow) <-
                paste(name, "obj", "type", tmprow$effectNumber, sep='.')
           # if ('requested' %in% names(myeff))
           # {
            #    tmprow$requested <- TRUE
           # }
            myeff <- rbind(myeff, tmprow)
            ints <- tmprow
        }
        ints <- ints[1, ]
    }
    if (include)
    {
        intn <- myeff$effectNumber == ints$effectNumber
        myeff[intn, "include"] <- include
        myeff[intn, c("effect1", "effect2", "effect3")] <-
            c(effect1, effect2, effect3)
    }
    else
    {
        intn <- (myeff$effect1 == effect1) & (myeff$effect2 == effect2)
        if (effect3 > 0)
        {
            intn <- intn & (myeff$effect3 == effect3)
        }
        myeff[intn, "include"] <- FALSE
    }
    if (verbose)
    {
        print.data.frame(myeff[intn, c("name", "shortName", "type",
                                       "interaction1", "interaction2",
                                       "include", "effect1", "effect2",
                                       "effect3")])
    }
    myeff
}

##@setEffect DataCreate
setEffect <- function(myeff, shortName, parameter=0,
                      fix=FALSE, test=FALSE, initialValue=0,
                      timeDummy=",",
                      include=TRUE, name=myeff$name[1],
                      type="eval", interaction1="", interaction2="",
                       period=1, group=1, character=FALSE)
{
    if (!character)
    {
        shortName <- deparse(substitute(shortName))
    }
    use <- myeff$shortName == shortName &
    myeff$name == name &
    myeff$type == type &
    myeff$interaction1 == interaction1 &
    myeff$interaction2 == interaction2 &
    (is.na(myeff$period) | myeff$period == period) &
    myeff$group == group
    if (sum(use) == 0)
    {
        stop("Effect not found")
    }
    if (sum(use) > 1)
    {
        stop("Effect not unique")
    }
    myeff[use, "parm"] <- parameter
    myeff[use, "include"] <- include
    myeff[use, "fix"] <- fix
    myeff[use, "test"] <- test
    myeff[use, "initialValue"] <- initialValue
    myeff[use, "timeDummy"] <- timeDummy
    print.data.frame(myeff[use, c("name", "shortName", "type", "interaction1",
                       "interaction2", "include", "parm", "fix", "test",
                       "initialValue", "timeDummy", "period", "group")])
    myeff
}
