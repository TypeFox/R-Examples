# aRpsDCA
# Copyright (C) 2016 dwt | terminus data science, LLC
# <dwt [at] terminusdatascience.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.

# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
# USA

arps.t.el <- function(decl, q.limit)
{
    nlminb(10,   # initial guess 10 periods
            function (t) ((arps.q(decl, t) - q.limit) ^ 2), # cost function
            lower=0,   # minimum 0
            upper=1e6  # 1 million [time units]
    )$par
}

arps.eur <- function(decl, q.limit)
{
    arps.Np(decl, arps.t.el(decl, q.limit))
}
