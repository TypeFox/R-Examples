## rats.R --- 
## Author          : Claus Dethlefsen
## Created On      : Mon Mar 11 15:22:48 2002
## Last Modified By: Claus Dethlefsen
## Last Modified On: Tue Apr 19 07:25:25 2005
## Update Count    : 51
## Status          : Unknown, Use with caution!
###############################################################################
##
##    Copyright (C) 2002  Susanne Gammelgaard Bottcher, Claus Dethlefsen
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program; if not, write to the Free Software
##    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
######################################################################

## op <- par(ask = interactive(), mfrow = c(1,1))

data(rats)
rats.df      <- rats

cat("Draw the prior DAG.\n",
"To insert an arrow from node 'A' to node 'B',\n",
"first click node 'A' and then click node 'B'.\n",
"When the DAG is finished, click 'stop'\n",
"\n",
"Then, inspect the local probability distributions\n",
"by clicking on the nodes. Finish by clicking 'stop'\n")
rats  <- network(rats.df,specifygraph=TRUE,inspectprob=TRUE)


#save this rats object
rats.orig <- rats
rats.prior <- jointprior(rats,12)

rats <- getnetwork(learn(rats,rats.df,rats.prior))
rats.empty <- getnetwork(learn(network(rats.df),rats.df,rats.prior))
banlist(rats.empty) <- banlist(rats)
## transfer node positions
for (i in 1:size(rats)) nodes(rats.empty)[[i]]$position <- nodes(rats)[[i]]$position

printline()
cat("Now, draw your favorite network. Notice how the\n",
    "network score changes. When bored, click stop\n",
    "and see how the search tries to find the network\n",
    "with highest score. The search algorithm is greedy\n",
    "search with random restart.\n")
newrat  <- getnetwork(drawnetwork(rats.empty,rats.df,rats.prior))


hiscorelist <- heuristic(newrat,rats.df,rats.prior,restart=10,degree=7,trace=TRUE)

op <- par(ask=TRUE)
cat("Now, we have tried out several networks\n")
cat("Ready to see the Hiscorelist?\n")

print(getnetwork(hiscorelist))
plot(getnetwork(hiscorelist))

par(op)
banlist(rats.empty) <- banlist(newrat)
for (i in 1:size(rats)) nodes(rats.empty)[[i]]$position <- nodes(newrat)[[i]]$position
allrats <- networkfamily(rats.df,rats.empty,rats.prior)
op <- par(ask=TRUE)
cat("We have now generated all",numbermixed(2,2),"networks\n")

print(getnetwork(allrats))
plot(nwfsort(getnetwork(allrats)))

par(op)
