## Copyright 2013-2014 Stefan Widgren and Maria Noremark,
## National Veterinary Institute, Sweden
##
## Licensed under the EUPL, Version 1.1 or - as soon they
## will be approved by the European Commission - subsequent
## versions of the EUPL (the "Licence");
## You may not use this work except in compliance with the
## Licence.
## You may obtain a copy of the Licence at:
##
## http://ec.europa.eu/idabc/eupl
##
## Unless required by applicable law or agreed to in
## writing, software distributed under the Licence is
## distributed on an "AS IS" basis,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
## express or implied.
## See the Licence for the specific language governing
## permissions and limitations under the Licence.

library(EpiContactTrace)

##
## Argument checking
##
## Check that the functions stop if not the expected arguments are
## given
##

##
## Missing parameters in call to IngoingContactChain
##
tools::assertError(IngoingContactChain(data.frame(source=1L,
                                                  destination=2L,
                                                  t=as.Date('2011-08-10')),
                                       tEnd='2005-10-31',
                                       days=90))

##
## Use either tEnd and days or inBegin, inEnd, outBegin and outEnd in
## call to NetworkSummary
##
tools::assertError(IngoingContactChain(data.frame(source=1L,
                                                  destination=2L,
                                                  t=as.Date('2011-08-10')),
                                       root=1,
                                       days=90))

##
## Use either tEnd and days or inBegin, inEnd, outBegin and outEnd in
## call to IngoingContactChain
##
tools::assertError(IngoingContactChain(data.frame(source=1L,
                                                  destination=2L,
                                                  t=as.Date('2011-08-10')),
                                       root=1,
                                       tEnd='2005-10-31'))

##
## Missing parameters in call to OutgoingContactChain
##
tools::assertError(OutgoingContactChain(data.frame(source=1L,
                                                   destination=2L,
                                                   t=as.Date('2011-08-10')),
                                        tEnd='2005-10-31',
                                        days=90))

##
## Use either tEnd and days or inBegin, inEnd, outBegin and outEnd in
## call to NetworkSummary
##
tools::assertError(OutgoingContactChain(data.frame(source=1L,
                                                   destination=2L,
                                                   t=as.Date('2011-08-10')),
                                        root=1,
                                        days=90))

##
## Use either tEnd and days or inBegin, inEnd, outBegin and outEnd in
## call to NetworkSummary
##
tools::assertError(OutgoingContactChain(data.frame(source=1L,
                                                   destination=2L,
                                                   t=as.Date('2011-08-10')),
                                        root=1,
                                        tEnd='2005-10-31'))

##
## Missing parameters in call to InDegree
##
tools::assertError(InDegree(data.frame(source=1L,
                                       destination=2L,
                                       t=as.Date('2011-08-10')),
                            tEnd='2005-10-31',
                            days=90))

##
## Use either tEnd and days or inBegin, inEnd, outBegin and outEnd in
## call to NetworkSummary
##
tools::assertError(InDegree(data.frame(source=1L,
                                       destination=2L,
                                       t=as.Date('2011-08-10')),
                            root=1,
                            days=90))

##
## Use either tEnd and days or inBegin, inEnd, outBegin and outEnd in
## call to NetworkSummary
##
tools::assertError(InDegree(data.frame(source=1L,
                                       destination=2L,
                                       t=as.Date('2011-08-10')),
                            root=1,
                            tEnd='2005-10-31'))

##
## Missing parameters in call to OutDegree
##
tools::assertError(OutDegree(data.frame(source=1L,
                                        destination=2L,
                                        t=as.Date('2011-08-10')),
                             tEnd='2005-10-31',
                             days=90))

##
## Use either tEnd and days or inBegin, inEnd, outBegin and outEnd in
## call to NetworkSummary
##
tools::assertError(OutDegree(data.frame(source=1L,
                                        destination=2L,
                                        t=as.Date('2011-08-10')),
                             root=1,
                             days=90))

##
## Use either tEnd and days or inBegin, inEnd, outBegin and outEnd in
## call to NetworkSummary
##
tools::assertError(OutDegree(data.frame(source=1L,
                                        destination=2L,
                                        t=as.Date('2011-08-10')),
                             root=1,
                             tEnd='2005-10-31'))

##
## Missing parameters in call to Trace
##
tools::assertError(Trace())

##
## Missing parameters in call to Trace
##
tools::assertError(Trace(movement=data.frame(
                             source=1L,
                             destination=2L,
                             t=as.Date('2011-08-10')),
                         tEnd='2005-10-31',
                         days=90))

##
## Use either tEnd and days or inBegin, inEnd, outBegin and outEnd in
## call to Trace
##
tools::assertError(Trace(movement=data.frame(
                             source=1L,
                             destination=2L,
                             t=as.Date('2011-08-10')),
                         root=1,
                         days=90))
##
## Use either tEnd and days or inBegin, inEnd, outBegin and outEnd in call to Trace
##
tools::assertError(Trace(movement=data.frame(
                             source=1L,
                             destination=2L,
                             t=as.Date('2011-08-10')),
                         root=1,
                         tEnd='2005-10-31'))

##
## movements must be a data.frame
##
tools::assertError(Trace(movements=1:3,
                         root=1L,
                         inBegin=as.Date('2011-08-10'),
                         inEnd=as.Date('2011-08-10'),
                         outBegin=as.Date('2011-08-10'),
                         outEnd=as.Date('2011-08-10')))

##
## movements must contain the columns source, destination and t
##
tools::assertError(Trace(movements=data.frame(destination=1, t=1),
                         root=1L,
                         inBegin=as.Date('2011-08-10'),
                         inEnd=as.Date('2011-08-10'),
                         outBegin=as.Date('2011-08-10'),
                         outEnd=as.Date('2011-08-10')))

##
## movements must contain the columns source, destination and t
##
tools::assertError(Trace(movements=data.frame(source=1, t=1),
                         root=1L,
                         inBegin=as.Date('2011-08-10'),
                         inEnd=as.Date('2011-08-10'),
                         outBegin=as.Date('2011-08-10'),
                         outEnd=as.Date('2011-08-10')))

##
## movements must contain the columns source, destination and t
##
tools::assertError(Trace(movements=data.frame(source=1, destination=1),
                         root=1L,
                         inBegin=as.Date('2011-08-10'),
                         inEnd=as.Date('2011-08-10'),
                         outBegin=as.Date('2011-08-10'),
                         outEnd=as.Date('2011-08-10')))

##
## invalid class of column t in movements
##
tools::assertError(Trace(movements=data.frame(source=1L, destination=2L, t=1),
                         root=1L,
                         inBegin=as.Date('2011-08-10'),
                         inEnd=as.Date('2011-08-10'),
                         outBegin=as.Date('2011-08-10'),
                         outEnd=as.Date('2011-08-10')))

##
## source in movements contains NA
##
tools::assertError(Trace(movements=data.frame(
                             source=c(1L,NA),
                             destination=c(2L,3L),
                             t=c('2011-08-10', '2011-08-10')),
                         root=1L,
                         inBegin=as.Date('2011-08-10'),
                         inEnd=as.Date('2011-08-10'),
                         outBegin=as.Date('2011-08-10'),
                         outEnd=as.Date('2011-08-10')))

##
## destination in movements contains NA
##
tools::assertError(Trace(movements=data.frame(
                             source=c(1L,2L),
                             destination=c(2L,NA),
                             t=c('2011-08-10', '2011-08-10')),
                         root=1L,
                         inBegin=as.Date('2011-08-10'),
                         inEnd=as.Date('2011-08-10'),
                         outBegin=as.Date('2011-08-10'),
                         outEnd=as.Date('2011-08-10')))

##
## t in movements contains NA
##
tools::assertError(Trace(movements=data.frame(
                             source=c(1L,2L),
                             destination=c(2L,3L),
                             t=c('2011-08-10', NA)),
                         root=1L,
                         inBegin=as.Date('2011-08-10'),
                         inEnd=as.Date('2011-08-10'),
                         outBegin=as.Date('2011-08-10'),
                         outEnd=as.Date('2011-08-10')))

##
## t in movements contains NA
##
tools::assertError(Trace(movements=data.frame(
                             source=c(1L,2L),
                             destination=c(2L,3L),
                             t=c('2011-08-10', NA)),
                         root=1L,
                         inBegin=as.Date('2011-08-10'),
                         inEnd=as.Date('2011-08-10'),
                         outBegin=as.Date('2011-08-10'),
                         outEnd=as.Date('2011-08-10')))

##
## invalid class of column n in movements
##
tools::assertError(Trace(movements=data.frame(
                             source=1L,
                             destination=2L,
                             t=as.Date('2011-08-10'),
                             n='3'),
                         root=1L,
                         inBegin=as.Date('2011-08-10'),
                         inEnd=as.Date('2011-08-10'),
                         outBegin=as.Date('2011-08-10'),
                         outEnd=as.Date('2011-08-10')))

##
## invalid class of column id in movements
##
tools::assertError(Trace(movements=data.frame(
                             source=1L,
                             destination=2L,
                             t=as.Date('2011-08-10'),
                             n=3L,
                             id=4),
                         root=1L,
                         inBegin=as.Date('2011-08-10'),
                         inEnd=as.Date('2011-08-10'),
                         outBegin=as.Date('2011-08-10'),
                         outEnd=as.Date('2011-08-10')))

##
## invalid class of column category in movements
##
tools::assertError(Trace(movements=data.frame(
                             source=1L,
                             destination=2L,
                             t=as.Date('2011-08-10'),
                             n=3L,
                             id=4L,
                             category=3),
                         root=1L,
                         inBegin=as.Date('2011-08-10'),
                         inEnd=as.Date('2011-08-10'),
                         outBegin=as.Date('2011-08-10'),
                         outEnd=as.Date('2011-08-10')))

##
## 'root' must be an integer or character
##
tools::assertError(Trace(movements=data.frame(
                             source=1L,
                             destination=2L,
                             t=as.Date('2011-08-10')),
                         root=1.1,
                         inBegin=as.Date('2011-08-10'),
                         inEnd=as.Date('2011-08-10'),
                         outBegin=as.Date('2011-08-10'),
                         outEnd=as.Date('2011-08-10')))

##
## 'inBegin' must be a Date vector
##
tools::assertError(Trace(movements=data.frame(
                             source=1L,
                             destination=2L,
                             t=as.Date('2011-08-10')),
                         root=1L,
                         inBegin=2011,
                         inEnd=as.Date('2011-08-10'),
                         outBegin=as.Date('2011-08-10'),
                         outEnd=as.Date('2011-08-10')))

##
## 'inEnd' must be a Date vector
##
tools::assertError(Trace(movements=data.frame(
                             source=1L,
                             destination=2L,
                             t=as.Date('2011-08-10')),
                         root=1L,
                         inBegin=as.Date('2011-08-10'),
                         inEnd=2011,
                         outBegin=as.Date('2011-08-10'),
                         outEnd=as.Date('2011-08-10')))

##
## 'outBegin' must be a Date vector
##
tools::assertError(Trace(movements=data.frame(
                             source=1L,
                             destination=2L,
                             t=as.Date('2011-08-10')),
                         root=1L,
                         inBegin=as.Date('2011-08-10'),
                         inEnd=as.Date('2011-08-10'),
                         outBegin=2011,
                         outEnd=as.Date('2011-08-10')))

##
## 'outEnd' must be a Date vector
##
tools::assertError(Trace(movements=data.frame(
                             source=1L,
                             destination=2L,
                             t=as.Date('2011-08-10')),
                         root=1L,
                         inBegin=as.Date('2011-08-10'),
                         inEnd=as.Date('2011-08-10'),
                         outBegin=as.Date('2011-08-10'),
                         outEnd=2011))

##
## inEnd < inBegin
##
tools::assertError(Trace(movements=data.frame(
                             source=1L,
                             destination=2L,
                             t=as.Date('2011-08-10')),
                         root=1L,
                         inBegin=as.Date('2011-08-10'),
                         inEnd=as.Date('2011-07-10'),
                         outBegin=as.Date('2011-08-10'),
                         outEnd=as.Date('2011-08-10')))

##
## outEnd < outBegin
##
tools::assertError(Trace(movements=data.frame(
                             source=1L,
                             destination=2L,
                             t=as.Date('2011-08-10')),
                         root=1L,
                         inBegin=as.Date('2011-08-10'),
                         inEnd=as.Date('2011-08-10'),
                         outBegin=as.Date('2011-08-10'),
                         outEnd=as.Date('2011-07-10')))
