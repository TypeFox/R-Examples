# /*
# R-Code
# S3 atomic 64bit integers for R
# (c) 2011 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2011-12-11
# Last changed:  2011-12-11
#*/

#! \name{bit64-package}
#! \alias{bit64-package}
#! \alias{bit64}
#! \alias{integer64}
#! \alias{is.integer64}
#! \alias{is.integer.integer64}
#! \alias{is.vector.integer64}
#! %as.vector.integer64 removed as requested by the CRAN maintainer \alias{as.vector.integer64}
#! \alias{length<-.integer64}
#! \alias{print.integer64}
#! \docType{package}
#! \title{
#!    A S3 class for vectors of 64bit integers
#! }
#! \description{
#! Package 'bit64' provides fast serializable S3 atomic 64bit (signed) integers 
#! that can be used in vectors, matrices, arrays and data.frames. Methods are 
#! available for coercion from and to logicals, integers, doubles, characters  
#! and factors as well as many elementwise and summary functions. 
#! \cr
#! \bold{Version 0.8}
#! With 'integer64' vectors you can store very large integers at the expense
#! of 64 bits, which is by factor 7 better than 'int64' from package 'int64'.
#! Due to the smaller memory footprint, the atomic vector architecture and  
#! using only S3 instead of S4 classes, most operations are one to three orders 
#! of magnitude faster: Example speedups are 4x for serialization, 250x for 
#! adding, 900x for coercion and 2000x for object creation. Also 'integer64' 
#! avoids an ongoing (potentially infinite) penalty for garbage collection
#! observed during existence of 'int64' objects (see code in example section). 
#! \cr
#! \bold{Version 0.9}
#! Package 'bit64' - which extends R with fast 64-bit integers - now has fast
#! (single-threaded) implementations the most important univariate algorithmic 
#! operations (those based on hashing and sorting). We now have methods for 
#! 'match', '%in%', 'duplicated', 'unique', 'table', 'sort', 'order', 'rank', 
#! 'quantile', 'median' and 'summary'. Regarding data management we also have 
#! novel generics 'unipos' (positions of the unique values), 'tiepos' (
#! positions of ties), 'keypos' (positions of foreign keys in a sorted 
#! dimension table) and derived methods 'as.factor' and 'as.ordered'. This 64-
#! bit functionality is implemented carefully to be not slower than the 
#! respective 32-bit operations in Base R and also to avoid outlying waiting 
#! times observed with 'order', 'rank' and 'table' (speedup factors 20/16/200 
#! respective). This increases the dataset size with wich we can work truly 
#! interactive. The speed is achieved by simple heuristic optimizers in high-
#! level functions choosing the best from multiple low-level algorithms and 
#! further taking advantage of a novel caching if activated. In an example R 
#! session using a couple of these operations the 64-bit integers performed 22x
#!  faster than base 32-bit integers, hash-caching improved this to 24x, 
#! sortorder-caching was most efficient with 38x (caching hashing and sorting 
#! is not worth it with 32x at duplicated RAM consumption).
#! }
#! \usage{
#!  integer64(length)
#!  \method{is}{integer64}(x)
#!  \method{length}{integer64}(x) <- value
#!  \method{print}{integer64}(x, quote=FALSE, \dots)
#! }
#! \arguments{
#!   \item{length}{ length of vector using \code{\link{integer}} }
#!   \item{x}{ an integer64 vector }
#!   \item{value}{ an integer64 vector of values to be assigned }
#!   \item{quote}{ logical, indicating whether or not strings should be printed with surrounding quotes. }
#!   \item{\dots}{ further arguments to the \code{\link{NextMethod}} }
#! }
#! \details{
#! \tabular{ll}{
#!    Package: \tab bit64\cr
#!    Type: \tab Package\cr
#!    Version: \tab 0.5.0\cr
#!    Date: \tab 2011-12-12\cr
#!    License: \tab GPL-2\cr
#!    LazyLoad: \tab yes\cr
#!    Encoding: \tab latin1\cr
#! }
#! }
#! \section{Design considerations}{
#!   64 bit integers are related to big data: we need them to overcome address space limitations. 
#!   Therefore performance of the 64 bit integer type is critical. 
#!   In the S language -- designed in 1975 -- atomic objects were defined to be vectors for a couple of good reasons:
#!   simplicity, option for implicit parallelization, good cache locality. 
#!   In recent years many analytical databases have learnt that lesson: column based data bases provide superior performance
#!   for many applications, the result are products such as MonetDB, Sybase IQ, Vertica, Exasol, Ingres Vectorwise.
#!   If we introduce 64 bit integers not natively in Base R but as an external package, we should at least strive to 
#!   make them as 'basic' as possible. Therefore the design choice of bit64 not only differs from \code{int64}, it is obvious: 
#!   Like the other atomic types in Base R, we model data type 'integer64' as a contiguous \code{\link{atomic}} vector in memory, 
#!   and we use the more basic \code{\link{S3}} class system, not \code{\link{S4}}. Like package \code{int64} we want our 'integer64' to be \code{\link{serialize}able}, 
#!   therefore we also use an existing data type as the basis. Again the choice is obvious: R has only one 64 bit data type: doubles.
#!   By using \code{\link{double}s}, \code{integer64} \code{\link{inherits}} some functionality such as \code{\link{is.atomic}}, \code{\link{length}}, 
#!   \code{\link{length<-}}, \code{\link{names}}, \code{\link{names<-}}, \code{\link{dim}}, \code{\link{dim<-}}, \code{\link{dimnames}}, \code{\link{dimnames}}.
#!   \cr
#!   Our R level functions strictly follow the functional programming paragdim: 
#!   no modification of arguments or other sideffects. Before version 0.93  we internally deviated from the strict paradigm
#!   in order to boost performance. Our C functions do not create new return values, 
#!   instead we pass-in the memory to be returned as an argument. This gives us the freedom to apply the C-function 
#!   to new or old vectors, which helps to avoid unnecessary memory allocation, unnecessary copying and unnessary garbage collection.
#!   Prior to 0.93 \emph{within} our R functions we also deviated from conventional R programming by not using \code{\link{attr<-}} and \code{\link{attributes<-}} 
#!   because they always did new memory allocation and copying in older R versions. If we wanted to set attributes of return values that we have freshly created,
#!   we instead used functions \code{\link[bit]{setattr}} and \code{\link[bit]{setattributes}} from package \code{\link[bit]{bit}}. 
#!   From version 0.93 \code{\link[bit]{setattr}} is only used for manipulating \code{\link{cache}} objects, in \code{\link{ramsort.integer64}} and \code{\link{sort.integer64}} and in \code{\link{as.data.frame.integer64}}.
#! }
#! \section{Arithmetic precision and coercion}{
#!   The fact that we introduce 64 bit long long integers -- without introducing 128-bit long doubles -- creates some subtle challenges:
#!   Unlike 32 bit \code{\link{integer}s}, the \code{integer64} are no longer a proper subset of \code{\link{double}}. 
#!   If a binary arithmetic operation does involve a \code{double} and a \code{integer}, it is a no-brainer to return \code{double} 
#!   without loss of information. If an \code{integer64} meets a \code{double}, it is not trivial what type to return. 
#!   Switching to \code{integer64} limits our ability to represent very large numbers, switching to \code{integer64} limits our ability 
#!   to distinguish \code{x} from \code{x+1}. Since the latter is purpose of introducing 64 bit integers, we usually return \code{integer64} 
#!   from functions involving \code{integer64}, for example in \code{\link[=c.integer64]{c}}, \code{\link[=cbind.integer64]{cbind}} 
#!   and \code{\link[=rbind.integer64]{rbind}}. 
#!   \cr
#!   Different from Base R, our operators \code{\link[=+.integer64]{+}}, 
#!   \code{\link[=-.integer64]{-}}, \code{\link[=\%/\%.integer64]{\%/\%}} and \code{\link[=\%\%.integer64]{\%\%}} coerce their arguments to 
#!   \code{integer64} and always return \code{integer64}. 
#!   \cr
#!   The multiplication operator \code{\link[=*.integer64]{*}} coerces its first argument to \code{integer64} 
#!   but allows its second argument to be also \code{double}: the second argument is internaly coerced to 'long double' 
#!   and the result of the multiplication is returned as \code{integer64}. 
#!   \cr
#!   The division \code{\link[=/.integer64]{/}} and power \code{\link[=^.integer64]{^}} operators also coerce their first argument to \code{integer64} 
#!   and coerce internally their second argument to 'long double', they return as \code{double}, like \code{\link[=sqrt.integer64]{sqrt}}, 
#!   \code{\link[=log.integer64]{log}}, \code{\link[=log2.integer64]{log2}} and \code{\link[=log10.integer64]{log10}} do. 
#! }
#! \section{Creating and testing S3 class 'integer64'}{
#!   Our creator function \code{integer64} takes an argument \code{length}, creates an atomic double vector of this length,
#!   attaches an S3 class attribute 'integer64' to it, and that's it. We simply rely on S3 method dispatch and interpret those 
#!   64bit elements as 'long long int'. 
#!   \cr
#!  \code{\link{is.double}} currently returns TRUE for \code{integer64} and might return FALSE in a later release.
#!  Consider \code{is.double} to have undefined behaviour and do query \code{is.integer64} \emph{before} querying \code{is.double}.
#! %As a second line of defense against misinterpretation we make \code{\link{is.double}}
#! %return \code{FALSE} by making it S3 generic and adding a method \code{\link{as.double.integer64}}. 
#!   The methods \code{\link{is.integer64}} and \code{\link{is.vector}} both return \code{TRUE} for \code{integer64}. 
#!  Note that we did not patch \code{\link{storage.mode}} and \code{\link{typeof}}, which both continue returning 'double' 
#!  Like for 32 bit \code{\link{integer}}, \code{\link{mode}} returns 'numeric' and \code{\link{as.double}}) tries coercing to \code{\link{double}}).
#!  It is likely that 'integer64' becomes a \code{\link[ff]{vmode}} in package \code{\link[ff]{ff}}. 
#!  \cr
#!  Further methods for creating \code{integer64} are \code{\link[=range.integer64]{range}} which returns the range of the data type if calles without arguments,
#!  \code{\link[=rep.integer64]{rep}}, \code{\link[=seq.integer64]{seq}}. 
#!  \cr
#!  For all available methods on \code{integer64} vectors see the index below and the examples.
#! }
#! \section{Index of implemented methods}{
#! \tabular{rrl}{
#!    \bold{creating,testing,printing} \tab \bold{see also}          \tab \bold{description} \cr
#!    \code{NA_integer64_} \tab \code{\link{NA_integer_}} \tab NA constant \cr
#!    \code{integer64} \tab \code{\link{integer}} \tab create zero atomic vector \cr
#!    \code{\link{rep.integer64}} \tab \code{\link{rep}} \tab  \cr
#!    \code{\link{seq.integer64}} \tab \code{\link{seq}} \tab  \cr
#!    \code{\link{is.integer64}} \tab \code{\link{is}} \tab  \cr
#!                                      \tab \code{\link{is.integer}} \tab inherited from Base R \cr
#!    %\code{\link{is.double.integer64}} \tab \code{\link{is.double}} \tab  \cr
#!    \code{\link{is.vector.integer64}} \tab \code{\link{is.vector}} \tab  \cr
#!    \code{\link{identical.integer64}} \tab \code{\link{identical}} \tab  \cr
#!    \code{\link{length<-.integer64}} \tab \code{\link{length<-}} \tab  \cr
#!                                      \tab \code{\link{length}} \tab inherited from Base R \cr
#!                                      \tab \code{\link{names<-}} \tab inherited from Base R \cr
#!                                      \tab \code{\link{names}} \tab inherited from Base R \cr
#!                                      \tab \code{\link{dim<-}} \tab inherited from Base R \cr
#!                                      \tab \code{\link{dim}} \tab inherited from Base R \cr
#!                                      \tab \code{\link{dimnames<-}} \tab inherited from Base R \cr
#!                                      \tab \code{\link{dimnames}} \tab inherited from Base R \cr
#!                                     \tab \code{\link{str}} \tab inherited from Base R, does not print values correctly \cr
#!    \code{\link{print.integer64}} \tab \code{\link{print}} \tab  \cr
#!  \cr
#!    \bold{coercing to integer64} \tab \bold{see also}          \tab \bold{description} \cr
#!    \code{\link{as.integer64}} \tab   \tab generic \cr
#!    \code{\link{as.integer64.character}} \tab \code{\link{character}} \tab  \cr
#!    \code{\link{as.integer64.double}} \tab \code{\link{double}} \tab  \cr
#!    \code{\link{as.integer64.integer}} \tab \code{\link{integer}} \tab  \cr
#!    \code{\link{as.integer64.integer64}} \tab \code{integer64} \tab  \cr
#!    \code{\link{as.integer64.logical}} \tab \code{\link{logical}} \tab  \cr
#!    \code{\link{as.integer64.NULL}} \tab \code{\link{NULL}} \tab  \cr
#!  \cr
#!    \bold{coercing from integer64} \tab \bold{see also}          \tab \bold{description} \cr
#!    \code{\link{as.bitstring}} \tab \code{\link{as.bitstring}} \tab generic \cr
#!    \code{\link{as.bitstring.integer64}} \tab  \tab  \cr
#!    \code{\link{as.character.integer64}} \tab \code{\link{as.character}} \tab  \cr
#!    \code{\link{as.double.integer64}} \tab \code{\link{as.double}} \tab  \cr
#!    \code{\link{as.integer.integer64}} \tab \code{\link{as.integer}} \tab  \cr
#!    \code{\link{as.logical.integer64}} \tab \code{\link{as.logical}} \tab  \cr
#!    %as.vector.integer64 removed as requested by the CRAN maintainer \code{\link{as.vector.integer64}} \tab \code{\link{as.vector}} \tab  \cr
#!  \cr
#!    \bold{data structures} \tab \bold{see also}          \tab \bold{description} \cr
#!    \code{\link{c.integer64}} \tab \code{\link{c}} \tab vector concatenate \cr
#!    \code{\link{cbind.integer64}} \tab \code{\link{cbind}} \tab column bind \cr
#!    \code{\link{rbind.integer64}} \tab \code{\link{rbind}} \tab row bind \cr
#!    \code{\link{as.data.frame.integer64}} \tab \code{\link{as.data.frame}} \tab coerce atomic object to data.frame \cr
#!                                          \tab \code{\link{data.frame}} \tab inherited from Base R since we have coercion \cr
#!  \cr
#!    \bold{subscripting} \tab \bold{see also}          \tab \bold{description} \cr
#!    \code{\link{[.integer64}} \tab \code{\link{[}} \tab vector and array extract \cr
#!    \code{\link{[<-.integer64}} \tab \code{\link{[<-}} \tab vector and array assign \cr
#!    \code{\link{[[.integer64}} \tab \code{\link{[[}} \tab scalar extract \cr
#!    \code{\link{[[<-.integer64}} \tab \code{\link{[[<-}} \tab scalar assign \cr
#!  \cr
#!    \bold{binary operators} \tab \bold{see also}          \tab \bold{description} \cr
#!    \code{\link{+.integer64}} \tab \code{\link{+}} \tab returns integer64 \cr
#!    \code{\link{-.integer64}} \tab \code{\link{-}} \tab returns integer64 \cr
#!    \code{\link{*.integer64}} \tab \code{\link{*}} \tab returns integer64 \cr
#!    \code{\link{^.integer64}} \tab \code{\link{^}} \tab returns double \cr
#!    \code{\link{/.integer64}} \tab \code{\link{/}} \tab returns double \cr
#!    \code{\link{\%/\%.integer64}} \tab \code{\link{\%/\%}} \tab returns integer64 \cr
#!    \code{\link{\%\%.integer64}} \tab \code{\link{\%\%}} \tab returns integer64 \cr
#!  \cr
#!    \bold{comparison operators} \tab \bold{see also}          \tab \bold{description} \cr
#!    \code{\link{==.integer64}} \tab \code{\link{==}} \tab  \cr
#!    \code{\link{!=.integer64}} \tab \code{\link{!=}} \tab  \cr
#!    \code{\link{<.integer64}} \tab \code{\link{<}} \tab  \cr
#!    \code{\link{<=.integer64}} \tab \code{\link{<=}} \tab  \cr
#!    \code{\link{>.integer64}} \tab \code{\link{>}} \tab  \cr
#!    \code{\link{>=.integer64}} \tab \code{\link{>=}} \tab  \cr
#!  \cr
#!    \bold{logical operators} \tab \bold{see also}          \tab \bold{description} \cr
#!    \code{\link{!.integer64}} \tab \code{\link{!}} \tab  \cr
#!    \code{\link{&.integer64}} \tab \code{\link{&}} \tab  \cr
#!    \code{\link{|.integer64}} \tab \code{\link{|}} \tab  \cr
#!    \code{\link{xor.integer64}} \tab \code{\link{xor}} \tab  \cr
#!  \cr
#!    \bold{math functions} \tab \bold{see also}          \tab \bold{description} \cr
#!    \code{\link{is.na.integer64}} \tab \code{\link{is.na}} \tab returns logical \cr
#!    \code{\link{format.integer64}} \tab \code{\link{format}} \tab returns character \cr
#!    \code{\link{abs.integer64}} \tab \code{\link{abs}} \tab returns integer64 \cr
#!    \code{\link{sign.integer64}} \tab \code{\link{sign}} \tab returns integer64 \cr
#!    \code{\link{log.integer64}} \tab \code{\link{log}} \tab returns double \cr
#!    \code{\link{log10.integer64}} \tab \code{\link{log10}} \tab  returns double \cr
#!    \code{\link{log2.integer64}} \tab \code{\link{log2}} \tab  returns double \cr
#!    \code{\link{sqrt.integer64}} \tab \code{\link{sqrt}} \tab  returns double \cr
#!    \code{\link{ceiling.integer64}} \tab \code{\link{ceiling}} \tab dummy returning its argument \cr
#!    \code{\link{floor.integer64}} \tab \code{\link{floor}} \tab dummy returning its argument \cr
#!    \code{\link{trunc.integer64}} \tab \code{\link{trunc}} \tab dummy returning its argument \cr
#!    \code{\link{round.integer64}} \tab \code{\link{round}} \tab dummy returning its argument \cr
#!    \code{\link{signif.integer64}} \tab \code{\link{signif}} \tab dummy returning its argument \cr
#!  \cr
#!    \bold{cumulative functions} \tab \bold{see also}          \tab \bold{description} \cr
#!    \code{\link{cummin.integer64}} \tab \code{\link{cummin}} \tab \cr
#!    \code{\link{cummax.integer64}} \tab \code{\link{cummax}} \tab \cr
#!    \code{\link{cumsum.integer64}} \tab \code{\link{cumsum}} \tab \cr
#!    \code{\link{cumprod.integer64}} \tab \code{\link{cumprod}} \tab \cr
#!    \code{\link{diff.integer64}} \tab \code{\link{diff}} \tab \cr
#!  \cr
#!    \bold{summary functions} \tab \bold{see also}          \tab \bold{description} \cr
#!    \code{\link{range.integer64}} \tab \code{\link{range}} \tab \cr
#!    \code{\link{min.integer64}} \tab \code{\link{min}} \tab  \cr
#!    \code{\link{max.integer64}} \tab \code{\link{max}} \tab  \cr
#!    \code{\link{sum.integer64}} \tab \code{\link{sum}} \tab  \cr
#!    \code{\link{mean.integer64}} \tab \code{\link{mean}} \tab  \cr
#!    \code{\link{prod.integer64}} \tab \code{\link{prod}} \tab  \cr
#!    \code{\link{all.integer64}} \tab \code{\link{all}} \tab  \cr
#!    \code{\link{any.integer64}} \tab \code{\link{any}} \tab  \cr
#!  \cr
#!    \bold{algorithmically complex functions} \tab \bold{see also}          \tab \bold{description (caching)}  \cr
#!    \code{\link{match.integer64}} \tab \code{\link{match}} \tab position of x in table (h//o/so) \cr
#!    \code{\link{\%in\%.integer64}} \tab \code{\link{\%in\%}} \tab is x in table? (h//o/so) \cr
#!    \code{\link{duplicated.integer64}} \tab \code{\link{duplicated}} \tab is current element duplicate of previous one? (h//o/so) \cr
#!    \code{\link{unique.integer64}} \tab \code{\link{unique}} \tab (shorter) vector of unique values only (h/s/o/so) \cr
#!    \code{\link{unipos.integer64}} \tab \code{\link{unipos}} \tab positions corresponding to unique values (h/s/o/so) \cr
#!    \code{\link{tiepos.integer64}} \tab \code{\link{tiepos}} \tab positions of values that are tied (//o/so) \cr
#!    \code{\link{keypos.integer64}} \tab \code{\link{keypos}} \tab position of current value in sorted list of unique values (//o/so) \cr
#!    \code{\link{as.factor.integer64}} \tab \code{\link{as.factor}} \tab convert to (unordered) factor with sorted levels of previous values (//o/so) \cr
#!    \code{\link{as.ordered.integer64}} \tab \code{\link{as.ordered}} \tab convert to ordered factor with sorted levels of previous values (//o/so) \cr
#!    \code{\link{table.integer64}} \tab \code{\link{table}} \tab unique values and their frequencies (h/s/o/so) \cr
#!    \code{\link{sort.integer64}} \tab \code{\link{sort}} \tab sorted vector (/s/o/so) \cr
#!    \code{\link{order.integer64}} \tab \code{\link{order}} \tab positions of elements that would create sorted vector (//o/so) \cr
#!    \code{\link{rank.integer64}} \tab \code{\link{rank}} \tab (average) ranks of non-NAs, NAs kept in place (/s/o/so) \cr
#!    \code{\link{quantile.integer64}} \tab \code{\link{quantile}} \tab (existing) values at specified percentiles (/s/o/so) \cr
#!    \code{\link{median.integer64}} \tab \code{\link{median}} \tab (existing) value at percentile 0.5 (/s/o/so) \cr
#!    \code{\link{summary.integer64}} \tab \code{\link{summary}} \tab  (/s/o/so) \cr
#!  \cr
#!    \bold{helper functions} \tab \bold{see also}          \tab \bold{description} \cr
#!    \code{\link{minusclass}} \tab \code{\link{minusclass}} \tab removing class attritbute \cr
#!    \code{\link{plusclass}} \tab \code{\link{plusclass}} \tab inserting class attribute \cr
#!    \code{\link{binattr}} \tab \code{\link{binattr}} \tab define binary op behaviour \cr
#!  \cr
#!    \bold{tested I/O functions} \tab \bold{see also}          \tab \bold{description} \cr
#!                                \tab \code{\link{read.table}} \tab inherited from Base R \cr
#!                                \tab \code{\link{write.table}} \tab inherited from Base R \cr
#!                                \tab \code{\link{serialize}} \tab inherited from Base R \cr
#!                                \tab \code{\link{unserialize}} \tab inherited from Base R \cr
#!                                \tab \code{\link{save}} \tab inherited from Base R \cr
#!                                \tab \code{\link{load}} \tab inherited from Base R \cr
#!                                \tab \code{\link{dput}} \tab inherited from Base R \cr
#!                                \tab \code{\link{dget}} \tab inherited from Base R \cr
#! }
#! }
#! \section{Limitations inherited from implementing 64 bit integers via an external package}{
#!   \itemize{
#!     \item \bold{vector size} of atomic vectors is still limited to \code{\link{.Machine}$integer.max}. 
#!     However, external memory extending packages such as \code{\link[ff]{ff}} or \code{bigmemory} 
#!     can extend their address space now with \code{integer64}. Having 64 bit integers also help 
#!     with those not so obvious address issues that arise once we exchange data with SQL databases 
#!     and datawarehouses, which use big integers as surrogate keys, e.g. on indexed primary key columns.
#!     This puts R into a relatively strong position compared to certain commercial statistical 
#!     softwares, which sell database connectivity but neither have the range of 64 bit integers, 
#!     nor have integers at all, nor have a single numeric data type in their macro-glue-language.
#!
#!     \item \bold{literals} such as \code{123LL} would require changes to Base R, up to then we need to write (and call) 
#!     \code{as.integer64(123L)} or \code{as.integer64(123)} or \code{as.integer64('123')}. 
#!     Only the latter allows to specify numbers beyond Base R's numeric data types and therefore is the recommended
#!     way to use -- using only one way may facilitate migrating code to literals at a later stage.
#!
#!   }
#! }
#! \section{Limitations inherited from Base R, Core team, can you change this?}{
#!   \itemize{
#!     \item \bold{\code{\link{identical}}} with default parameters does not distinguish all bit-patterns of doubles. 
#!     For testing purposes we provide a wrapper \code{\link{identical.integer64}} that will distinguish all bit-patterns.
#!     It would be desireable to have a single call of \code{\link{identical}} handle both, \code{\link{double}} and \code{integer64}.
#! 
#!     \item the \bold{colon} operator \code{\link{:}} officially does not dispatches S3 methods, however, we have made it generic
#!      \preformatted{
#!      from <- lim.integer64()[1]
#!      to <- from+99
#!      from:to
#!    }
#!    As a limitation remains: it will only dispatch at its first argument \code{from} but not at its second \code{to}.
#! 
#!     \item \bold{\code{\link{is.double}}} does not dispatches S3 methods, However, we have made it generic 
#!		and it will return \code{FALSE} on \code{integer64}.
#!
#!     \item \bold{\code{\link{c}}} only dispatches \code{\link{c.integer64}} if the first argument is \code{integer64}
#!     and it does not recursively dispatch the proper method when called with argument \code{recursive=TRUE}
#!     Therefore \preformatted{
#!       c(list(integer64,integer64))
#!     }
#!      does not work and for now you can only call \preformatted{
#!        c.integer64(list(x,x))
#!      }
#!
#!     \item \bold{generic binary operators} fail to dispatch *any* user-defined S3 method 
#!     if the two arguments have two different S3 classes. For example we have two classes 
#!     \code{\link{bit}} and \code{\link{bitwhich}} sparsely representing boolean vectors 
#!     and we have methods \code{\link{&.bit}} and \code{\link{&.bitwhich}}. For an expression
#!     involving both as in \code{ bit & bitwhich}, none of the two methods is dispatched. 
#!     Instead a standard method is dispatched, which neither handles \code{\link{bit}} 
#!     nor \code{\link{bitwhich}}. Although it lacks symmetry, the better choice would be to 
#!     dispatch simply the method of the class of the first argument in case of class conflict. 
#!     This choice would allow authors of extension packages providing coherent behaviour 
#!     at least within their contributed classes. But as long as none of the package authors 
#!     methods is dispatched, he cannot handle the conflicting classes at all.
#!
#!     \item \bold{\code{\link{unlist}}} is not generic and if it were, we would face similar problems as with \code{c()}
#!
#!     \item \bold{\code{\link{vector}}} with argument \code{mode='integer64'} cannot work without adjustment of Base R
#!     \item \bold{\code{\link{as.vector}}} with argument \code{mode='integer64'} cannot work without adjustment of Base R
#!
#!     \item \bold{\code{\link{is.vector}}} does not dispatch its method \code{\link{is.vector.integer64}}
#!
#!     \item \bold{\code{\link{mode<-}}} drops the class 'integer64' which is returned from \code{as.integer64}.
#!        Also it does not remove an existing class 'integer64' when assigning mode 'integer'. 
#!
#!     \item \bold{\code{\link{storage.mode<-}}} does not support external data types such as \code{as.integer64}
#!
#!     \item \bold{\code{\link{matrix}}} does drop the 'integer64' class attribute.
#!
#!     \item \bold{\code{\link{array}}}  does drop the 'integer64' class attribute. 
#!            In current R versions (1.15.1) this can be circumvented by activating the function 
#!						\code{as.vector.integer64} further down this file.
#!						However, the CRAN maintainer has requested to remove \code{as.vector.integer64}, 
#! 						even at the price of breaking previously working functionality of the package. 
#!
#!     \item \bold{\code{\link{str}}} does not print the values of \code{integer64} correctly
#!
#!   }
#! }
#! \section{further limitations}{
#!   \itemize{
#!     \item \bold{subscripting} non-existing elements and subscripting with \code{NA}s is currently not supported. 
#!     Such subscripting currently returns \code{9218868437227407266} instead of \code{NA} (the \code{NA} value of the underlying double code).
#!     Following the full R behaviour here would either destroy performance or require extensive C-coding. 
#!   }
#! }
#! \value{
#!   \code{integer64} returns a vector of 'integer64', 
#!    i.e. a vector of \code{\link{double}} decorated with class 'integer64'.
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! Maintainer: Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \keyword{ package }
#! \keyword{ classes }
#! \keyword{ manip }
#! \seealso{ \code{\link{integer}} in base R }
#! \examples{
#! message("Using integer64 in vector")
#! x <- integer64(8)    # create 64 bit vector
#! x
#! is.atomic(x)         # TRUE
#! is.integer64(x)      # TRUE
#! is.numeric(x)        # TRUE
#! is.integer(x)        # FALSE - debatable
#! is.double(x)         # FALSE - might change
#! x[] <- 1:2           # assigned value is recycled as usual
#! x[1:6]               # subscripting as usual
#! length(x) <- 13      # changing length as usual
#! x
#! rep(x, 2)            # replicate as usual
#! seq(as.integer64(1), 10)     # seq.integer64 is dispatched on first given argument
#! seq(to=as.integer64(10), 1)  # seq.integer64 is dispatched on first given argument
#! seq.integer64(along.with=x)  # or call seq.integer64 directly
#! # c.integer64 is dispatched only if *first* argument is integer64 ...
#! x <- c(x,runif(length(x), max=100)) 
#! # ... and coerces everything to integer64 - including double
#! x                                   
#! names(x) <- letters  # use names as usual
#! x
#! 
#! message("Using integer64 in array - note that 'matrix' currently does not work")
#! message("as.vector.integer64 removed as requested by the CRAN maintainer")
#! message("as consequence 'array' also does not work anymore")
#! %y <- array(as.integer64(NA), dim=c(3,4), dimnames=list(letters[1:3], LETTERS[1:4]))
#! message("we still can create a matrix or array by assigning 'dim'")
#! y <- rep(as.integer64(NA), 12)
#! dim(y) <- c(3,4)
#! dimnames(y) <- list(letters[1:3], LETTERS[1:4])
#! y["a",] <- 1:2       # assigning as usual
#! y
#! y[1:2,-4]            # subscripting as usual
#! # cbind.integer64 dispatched on any argument and coerces everything to integer64
#! cbind(E=1:3, F=runif(3, 0, 100), G=c("-1","0","1"), y)
#! 
#! message("Using integer64 in data.frame")
#! str(as.data.frame(x))
#! str(as.data.frame(y))
#! str(data.frame(y))
#! str(data.frame(I(y)))
#! d <- data.frame(x=x, y=runif(length(x), 0, 100))
#! d
#! d$x
#! 
#! message("Using integer64 with csv files")
#! fi64 <- tempfile()
#! write.csv(d, file=fi64, row.names=FALSE)
#! e <- read.csv(fi64, colClasses=c("integer64", NA))
#! unlink(fi64)
#! str(e)
#! identical.integer64(d$x,e$x)
#! 
#! message("Serializing and unserializing integer64")
#! dput(d, fi64)
#! e <- dget(fi64)
#! identical.integer64(d$x,e$x)
#! e <- d[,]
#! save(e, file=fi64)
#! rm(e)
#! load(file=fi64)
#! identical.integer64(d,e)
#! 
#! ### A couple of unit tests follow hidden in a dontshow{} directive ###
#!   \dontshow{
#! message("Testing identical.integer64")
#! i64 <- as.double(NA); class(i64) <- "integer64"
#! stopifnot(identical(unclass(i64-1), unclass(i64+1)))
#! stopifnot(identical(i64-1, i64+1))
#! stopifnot(!identical.integer64(i64-1, i64+1))
#! 
#! message("Testing dispatch of 'c' method")
#! stopifnot(identical.integer64(c(integer64(0), NA), as.integer64(NA)))
#! message("Dispatch on the second argument fails and we want to be notified once that changes")
#! stopifnot(!identical.integer64(c(NA, integer64(0)), as.integer64(NA)))
#! 
#! message("Testing minus and plus")
#! d64 <- c(-.Machine$double.base^.Machine$double.digits, -.Machine$integer.max, -1, 0, 1, .Machine$integer.max, .Machine$double.base^.Machine$double.digits)
#! i64 <- as.integer64(d64)
#! stopifnot(identical.integer64(i64-1+1,i64))
#! stopifnot(identical.integer64(i64+1-1,i64))
#! 
#! message("Testing minus and plus edge cases and 'rev'\n")
#! stopifnot(identical.integer64(lim.integer64()+1-1, c(lim.integer64()[1], NA)))
#! stopifnot(identical.integer64(rev(lim.integer64())-1+1, c(lim.integer64()[2], NA)))
#! 
#! message("Testing 'range.integer64', multiplication and integer division")
#! i64 <- integer64(63)
#! i64[1] <- 1
#! for (i in 2:63)
#! 	i64[i] <- 2*i64[i-1]
#! stopifnot(identical.integer64(i64 * rev(i64), rep(i64[63], 63)))
#! for (i in 63:2)
#! 	i64[i-1] <- i64[i]\%/\%2
#! stopifnot(identical.integer64(i64 * rev(i64), rep(i64[63], 63)))
#! for (i in 63:2)
#! 	i64[i-1] <- i64[i]/2
#! stopifnot(identical.integer64(i64 * rev(i64), rep(i64[63], 63)))
#! stopifnot(identical.integer64(c( -i64[63] - (i64[63]-1), i64[63]+(i64[63]-1) ), lim.integer64()))
#! 
#! stopifnot(identical.integer64(i64[-1]\%/\%2*as.integer64(2), i64[-1]))
#! stopifnot(identical.integer64(i64[-1]\%/\%2L*as.integer64(2), i64[-1]))
#! stopifnot(identical.integer64(i64[-1]/2*as.integer64(2), i64[-1]))
#! stopifnot(identical.integer64(i64[-1]/2*as.integer64(2), i64[-1]))
#! 
#! stopifnot(identical.integer64(i64[-63]*2\%/\%2, i64[-63]))
#! stopifnot(identical.integer64(i64[-63]*2L\%/\%2L, i64[-63]))
#! stopifnot(identical.integer64(as.integer64(i64[-63]*2/2), i64[-63]))
#! stopifnot(identical.integer64(as.integer64(i64[-63]*2L/2L), i64[-63]))
#! 
#! message("Testing sqrt, power and log")
#! stopifnot(identical.integer64( as.integer64(sqrt(i64[-1][c(FALSE, TRUE)])*sqrt(i64[-1][c(FALSE, TRUE)])), i64[-1][c(FALSE, TRUE)] ))
#! 
#! stopifnot(identical.integer64(as.integer64(2)^(0:62), i64))
#! stopifnot(identical.integer64(as.integer64(0:62), as.integer64(round(log2(i64)))))
#! stopifnot(identical.integer64(as.integer64(round(log(as.integer64(2)^(0:62), 2))), as.integer64(0:62)))
#! stopifnot(identical.integer64(as.integer64(round(log(as.integer64(3)^(0:39), 3))), as.integer64(0:39)))
#! stopifnot(identical.integer64(as.integer64(round(log(as.integer64(10)^(0:18), 10))), as.integer64(0:18)))
#! stopifnot(identical.integer64(as.integer64(round(log10(as.integer64(10)^(0:18)))), as.integer64(0:18)))
#! 
#! stopifnot(identical.integer64((as.integer64(2)^(1:62))^(1/1:62), as.integer64(rep(2, 62))))
#! stopifnot(identical.integer64((as.integer64(3)^(1:39))^(1/1:39), as.integer64(rep(3, 39))))
#! stopifnot(identical.integer64((as.integer64(10)^(1:18))^(1/1:18), as.integer64(rep(10, 18))))
#! 
#! message("Testing c and rep")
#! stopifnot(identical.integer64( as.integer64(rep(1:3, 1:3)), rep(as.integer64(1:3), 1:3)))
#! stopifnot(identical.integer64( as.integer64(rep(1:3, 3)), rep(as.integer64(1:3), 3)))
#!  
#! x <- as.double(c(NA,NA,NA))
#! class(x) <- "integer64"
#! x <- x + -1:1
#! stopifnot(identical.integer64(rep(x, 3), c(x,x,x) ))
#! stopifnot(identical.integer64(c.integer64(list(x,x,x), recursive=TRUE), c(x,x,x) ))
#! 
#! message("Testing seq")
#! stopifnot(identical.integer64(seq(as.integer64(1), 10, 2), as.integer64(seq(1, 10, 2)) ))
#! stopifnot(identical.integer64(seq(as.integer64(1), by=2, length.out=5), as.integer64(seq(1, by=2, length.out=5)) ))
#! stopifnot(identical.integer64(seq(as.integer64(1), by=2, length.out=6), as.integer64(seq(1, by=2, length.out=6)) ))
#! stopifnot(identical.integer64(seq.integer64(along.with=3:5), as.integer64(seq(along.with=3:5)) ))
#! stopifnot(identical.integer64(seq(as.integer64(1), to=-9), as.integer64(seq(1, to=-9)) ))
#! 
#! message("Testing cbind and rbind")
#! stopifnot(identical.integer64( cbind(as.integer64(1:3), 1:3), {x <- rep(as.integer64(1:3), 2); dim(x)<-c(3,2);x}))
#! stopifnot(identical.integer64( rbind(as.integer64(1:3), 1:3), t({x <- rep(as.integer64(1:3), 2); dim(x)<-c(3,2);x})))
#! 
#! message("Testing coercion")
#! stopifnot(identical( as.double(as.integer64(c(NA, seq(0, 9, 0.25)))), as.double(as.integer(c(NA, seq(0, 9, 0.25))))))
#! stopifnot(identical( as.character(as.integer64(c(NA, seq(0, 9, 0.25)))), as.character(as.integer(c(NA, seq(0, 9, 0.25))))))
#! stopifnot(identical( as.integer(as.integer64(c(NA, seq(0, 9, 0.25)))), as.integer(c(NA, seq(0, 9, 0.25)))))
#! stopifnot(identical( as.logical(as.integer64(c(NA, seq(0, 9, 0.25)))), as.logical(as.integer(c(NA, seq(0, 9, 0.25))))))
#! stopifnot(identical( as.integer(as.integer64(c(NA, FALSE, TRUE))), as.integer(c(NA, FALSE, TRUE))))
#! stopifnot(identical( as.integer64(as.integer(as.integer64(-9:9))), as.integer64(-9:9)))
#! stopifnot(identical( as.integer64(as.double(as.integer64(-9:9))), as.integer64(-9:9)))
#! stopifnot(identical( as.integer64(as.character(as.integer64(-9:9))), as.integer64(-9:9)))
#! stopifnot(identical( as.integer64(as.character(lim.integer64())), lim.integer64()))
#!
#! message("-- testing logical operators --")
#! stopifnot(identical.integer64(!c(NA, -1:1), !c(as.integer64(NA), -1:1)))
#! stopifnot(identical.integer64(rep(c(NA, -1:1), 4)&rep(c(NA, -1:1), rep(4, 4)), as.integer64(rep(c(NA, -1:1), 4))&as.integer64(rep(c(NA, -1:1), rep(4, 4)))))
#! stopifnot(identical.integer64(rep(c(NA, -1:1), 4)|rep(c(NA, -1:1), rep(4, 4)), as.integer64(rep(c(NA, -1:1), 4))|as.integer64(rep(c(NA, -1:1), rep(4, 4)))))
#! stopifnot(identical.integer64(xor(rep(c(NA, -1:1), 4),rep(c(NA, -1:1), rep(4, 4))), xor(as.integer64(rep(c(NA, -1:1), 4)),as.integer64(rep(c(NA, -1:1), rep(4, 4))))))
#! 
#! message("-- testing comparison operators --")
#! stopifnot(identical.integer64(rep(c(NA, -1:1), 4)==rep(c(NA, -1:1), rep(4, 4)), as.integer64(rep(c(NA, -1:1), 4))==as.integer64(rep(c(NA, -1:1), rep(4, 4)))))
#! stopifnot(identical.integer64(rep(c(NA, -1:1), 4)!=rep(c(NA, -1:1), rep(4, 4)), as.integer64(rep(c(NA, -1:1), 4))!=as.integer64(rep(c(NA, -1:1), rep(4, 4)))))
#! stopifnot(identical.integer64(rep(c(NA, -1:1), 4)>rep(c(NA, -1:1), rep(4, 4)), as.integer64(rep(c(NA, -1:1), 4))>as.integer64(rep(c(NA, -1:1), rep(4, 4)))))
#! stopifnot(identical.integer64(rep(c(NA, -1:1), 4)>=rep(c(NA, -1:1), rep(4, 4)), as.integer64(rep(c(NA, -1:1), 4))>=as.integer64(rep(c(NA, -1:1), rep(4, 4)))))
#! stopifnot(identical.integer64(rep(c(NA, -1:1), 4)<rep(c(NA, -1:1), rep(4, 4)), as.integer64(rep(c(NA, -1:1), 4))<as.integer64(rep(c(NA, -1:1), rep(4, 4)))))
#! stopifnot(identical.integer64(rep(c(NA, -1:1), 4)<=rep(c(NA, -1:1), rep(4, 4)), as.integer64(rep(c(NA, -1:1), 4))<=as.integer64(rep(c(NA, -1:1), rep(4, 4)))))
#!
#! message("-- testing vector functions --")
#! stopifnot(identical.integer64( is.na(as.integer64(c(NA, -1:1))), is.na(c(NA, -1:1)) ))
#! stopifnot(identical.integer64( format(as.integer64(c(NA, -1:1))), format(c(NA, -1:1)) ))
#! stopifnot(identical.integer64( abs(as.integer64(c(NA, -1:1))), as.integer64(abs(c(NA, -1:1))) ))
#! stopifnot(identical.integer64( sign(as.integer64(c(NA, -1:1))), as.integer64(sign(c(NA, -1:1))) ))
#! stopifnot(identical.integer64( ceiling(as.integer64(c(NA, -1:1))), as.integer64(ceiling(c(NA, -1:1))) ))
#! stopifnot(identical.integer64( floor(as.integer64(c(NA, -1:1))), as.integer64(floor(c(NA, -1:1))) ))
#! stopifnot(identical.integer64( trunc(as.integer64(c(NA, -1:1))), as.integer64(trunc(c(NA, -1:1))) ))
#! stopifnot(identical.integer64( signif(as.integer64(c(NA, -1:1))), as.integer64(c(NA, -1:1)) ))
#!
#! message("Testing summary functions")
#! stopifnot(identical(all(as.integer(1)), all(as.integer64(1))))
#! stopifnot(identical(all(as.integer(0)), all(as.integer64(0))))
#! stopifnot(identical(all(as.integer(NA)), all(as.integer64(NA))))
#! stopifnot(identical(all(as.integer(NA), na.rm=TRUE), all(as.integer64(NA), na.rm=TRUE)))
#! stopifnot(identical(all(as.integer(1), NA), all(as.integer64(1), NA)))
#! stopifnot(identical(all(as.integer(0), NA), all(as.integer64(0), NA)))
#! stopifnot(identical(all(as.integer(1), NA, na.rm=TRUE), all(as.integer64(1), NA, na.rm=TRUE)))
#! stopifnot(identical(all(as.integer(0), NA, na.rm=TRUE), all(as.integer64(0), NA, na.rm=TRUE)))
#! stopifnot(identical(all(as.integer(c(1, NA))), all(as.integer64(c(1, NA)))))
#! stopifnot(identical(all(as.integer(c(0, NA))), all(as.integer64(c(0, NA)))))
#! stopifnot(identical(all(as.integer(c(1, NA)), na.rm=TRUE), all(as.integer64(c(1, NA)), na.rm=TRUE)))
#! stopifnot(identical(all(as.integer(c(0, NA)), na.rm=TRUE), all(as.integer64(c(0, NA)), na.rm=TRUE)))
#! 
#! stopifnot(identical(any(as.integer(1)), any(as.integer64(1))))
#! stopifnot(identical(any(as.integer(0)), any(as.integer64(0))))
#! stopifnot(identical(any(as.integer(NA)), any(as.integer64(NA))))
#! stopifnot(identical(any(as.integer(NA), na.rm=TRUE), any(as.integer64(NA), na.rm=TRUE)))
#! stopifnot(identical(any(as.integer(1), NA), any(as.integer64(1), NA)))
#! stopifnot(identical(any(as.integer(0), NA), any(as.integer64(0), NA)))
#! stopifnot(identical(any(as.integer(1), NA, na.rm=TRUE), any(as.integer64(1), NA, na.rm=TRUE)))
#! stopifnot(identical(any(as.integer(0), NA, na.rm=TRUE), any(as.integer64(0), NA, na.rm=TRUE)))
#! stopifnot(identical(any(as.integer(c(1, NA))), any(as.integer64(c(1, NA)))))
#! stopifnot(identical(any(as.integer(c(0, NA))), any(as.integer64(c(0, NA)))))
#! stopifnot(identical(any(as.integer(c(1, NA)), na.rm=TRUE), any(as.integer64(c(1, NA)), na.rm=TRUE)))
#! stopifnot(identical(any(as.integer(c(0, NA)), na.rm=TRUE), any(as.integer64(c(0, NA)), na.rm=TRUE)))
#! 
#! stopifnot(identical.integer64(as.integer64(sum(c(2, 3, NA))), sum(as.integer64(c(2, 3, NA)))))
#! stopifnot(identical.integer64(as.integer64(sum(c(2, 3, NA), na.rm=TRUE)), sum(as.integer64(c(2, 3, NA)), na.rm=TRUE)))
#! stopifnot(identical.integer64(as.integer64(sum(c(2, 3, NA))), sum(as.integer64(c(2, 3, NA)))))
#! stopifnot(identical.integer64(as.integer64(sum(c(2, 3, NA), na.rm=TRUE)), sum(as.integer64(c(2, 3, NA)), na.rm=TRUE)))
#! stopifnot(identical.integer64(as.integer64(sum(2, 3, NA)), sum(as.integer64(2), 3, NA)))
#! stopifnot(identical.integer64(as.integer64(sum(2, 3, NA, na.rm=TRUE)), sum(as.integer64(2), 3, NA, na.rm=TRUE)))
#! stopifnot(identical.integer64(as.integer64(sum(2, 3, NA)), sum(as.integer64(2), 3, NA)))
#! stopifnot(identical.integer64(as.integer64(sum(2, 3, NA, na.rm=TRUE)), sum(as.integer64(2), 3, NA, na.rm=TRUE)))
#! 
#! stopifnot(identical.integer64(as.integer64(prod(c(2, 3, NA))), prod(as.integer64(c(2, 3, NA)))))
#! stopifnot(identical.integer64(as.integer64(prod(c(2, 3, NA), na.rm=TRUE)), prod(as.integer64(c(2, 3, NA)), na.rm=TRUE)))
#! stopifnot(identical.integer64(as.integer64(prod(c(2, 3, NA))), prod(as.integer64(c(2, 3, NA)))))
#! stopifnot(identical.integer64(as.integer64(prod(c(2, 3, NA), na.rm=TRUE)), prod(as.integer64(c(2, 3, NA)), na.rm=TRUE)))
#! stopifnot(identical.integer64(as.integer64(prod(2, 3, NA)), prod(as.integer64(2), 3, NA)))
#! stopifnot(identical.integer64(as.integer64(prod(2, 3, NA, na.rm=TRUE)), prod(as.integer64(2), 3, NA, na.rm=TRUE)))
#! stopifnot(identical.integer64(as.integer64(prod(2, 3, NA)), prod(as.integer64(2), 3, NA)))
#! stopifnot(identical.integer64(as.integer64(prod(2, 3, NA, na.rm=TRUE)), prod(as.integer64(2), 3, NA, na.rm=TRUE)))
#! 
#! stopifnot(identical.integer64(as.integer64(min(c(2, 3, NA))), min(as.integer64(c(2, 3, NA)))))
#! stopifnot(identical.integer64(as.integer64(min(c(2, 3, NA), na.rm=TRUE)), min(as.integer64(c(2, 3, NA)), na.rm=TRUE)))
#! stopifnot(identical.integer64(as.integer64(min(c(2, 3, NA))), min(as.integer64(c(2, 3, NA)))))
#! stopifnot(identical.integer64(as.integer64(min(c(2, 3, NA), na.rm=TRUE)), min(as.integer64(c(2, 3, NA)), na.rm=TRUE)))
#! stopifnot(identical.integer64(as.integer64(min(2, 3, NA)), min(as.integer64(2), 3, NA)))
#! stopifnot(identical.integer64(as.integer64(min(2, 3, NA, na.rm=TRUE)), min(as.integer64(2), 3, NA, na.rm=TRUE)))
#! stopifnot(identical.integer64(as.integer64(min(2, 3, NA)), min(as.integer64(2), 3, NA)))
#! stopifnot(identical.integer64(as.integer64(min(2, 3, NA, na.rm=TRUE)), min(as.integer64(2), 3, NA, na.rm=TRUE)))
#! 
#! stopifnot(identical.integer64(as.integer64(max(c(2, 3, NA))), max(as.integer64(c(2, 3, NA)))))
#! stopifnot(identical.integer64(as.integer64(max(c(2, 3, NA), na.rm=TRUE)), max(as.integer64(c(2, 3, NA)), na.rm=TRUE)))
#! stopifnot(identical.integer64(as.integer64(max(c(2, 3, NA))), max(as.integer64(c(2, 3, NA)))))
#! stopifnot(identical.integer64(as.integer64(max(c(2, 3, NA), na.rm=TRUE)), max(as.integer64(c(2, 3, NA)), na.rm=TRUE)))
#! stopifnot(identical.integer64(as.integer64(max(2, 3, NA)), max(as.integer64(2), 3, NA)))
#! stopifnot(identical.integer64(as.integer64(max(2, 3, NA, na.rm=TRUE)), max(as.integer64(2), 3, NA, na.rm=TRUE)))
#! stopifnot(identical.integer64(as.integer64(max(2, 3, NA)), max(as.integer64(2), 3, NA)))
#! stopifnot(identical.integer64(as.integer64(max(2, 3, NA, na.rm=TRUE)), max(as.integer64(2), 3, NA, na.rm=TRUE)))
#! 
#! stopifnot(identical.integer64(as.integer64(range(c(2, 3, NA))), range(as.integer64(c(2, 3, NA)))))
#! stopifnot(identical.integer64(as.integer64(range(c(2, 3, NA), na.rm=TRUE)), range(as.integer64(c(2, 3, NA)), na.rm=TRUE)))
#! stopifnot(identical.integer64(as.integer64(range(c(2, 3, NA))), range(as.integer64(c(2, 3, NA)))))
#! stopifnot(identical.integer64(as.integer64(range(c(2, 3, NA), na.rm=TRUE)), range(as.integer64(c(2, 3, NA)), na.rm=TRUE)))
#! stopifnot(identical.integer64(as.integer64(range(2, 3, NA)), range(as.integer64(2), 3, NA)))
#! stopifnot(identical.integer64(as.integer64(range(2, 3, NA, na.rm=TRUE)), range(as.integer64(2), 3, NA, na.rm=TRUE)))
#! stopifnot(identical.integer64(as.integer64(range(2, 3, NA)), range(as.integer64(2), 3, NA)))
#! stopifnot(identical.integer64(as.integer64(range(2, 3, NA, na.rm=TRUE)), range(as.integer64(2), 3, NA, na.rm=TRUE)))
#! 
#! message("-- testing cummulative functions --")
#! stopifnot(identical.integer64(as.integer64(cumsum(c(2, 3, NA, 1, 4))), cumsum(as.integer64(c(2, 3, NA, 1, 4)))))
#! stopifnot(identical.integer64(as.integer64(cumprod(c(2, 3, NA, 1, 4))), cumprod(as.integer64(c(2, 3, NA, 1, 4)))))
#! stopifnot(identical.integer64(as.integer64(cummin(c(2, 3, NA, 1, 4))), cummin(as.integer64(c(2, 3, NA, 1, 4)))))
#! stopifnot(identical.integer64(as.integer64(cummax(c(2, 3, NA, 1, 4))), cummax(as.integer64(c(2, 3, NA, 1, 4)))))
#! 
#! message("testing diff")
#! d64 <- diffinv(rep(.Machine$integer.max, 100), lag=2, differences=2)
#! i64 <- as.integer64(d64)
#! identical(diff(d64, lag=2, differences=2), as.double(diff(i64, lag=2, differences=2)))
#!
#!   }
#!
#!   \dontrun{
#! message("== Differences between integer64 and int64 ==")
#! require(bit64)
#! require(int64)
#! 
#! message("-- integer64 is atomic --")
#! is.atomic(integer64())
#! #is.atomic(int64())
#! str(integer64(3))
#! #str(int64(3))
#! 
#! message("-- The following performance numbers are measured under RWin64  --")
#! message("-- under RWin32 the advantage of integer64 over int64 is smaller --")
#!
#! message("-- integer64 needs 7x/5x less RAM than int64 under 64/32 bit OS 
#! (and twice the RAM of integer as it should be) --")
#! #as.vector(object.size(int64(1e6))/object.size(integer64(1e6)))
#! as.vector(object.size(integer64(1e6))/object.size(integer(1e6)))
#! 
#! message("-- integer64 creates 2000x/1300x faster than int64 under 64/32 bit OS
#! (and 3x the time of integer) --")
#! t32 <- system.time(integer(1e8))
#! t64 <- system.time(integer64(1e8))
#! #T64 <- system.time(int64(1e7))*10  # using 1e8 as above stalls our R on an i7 8 GB RAM Thinkpad
#! #T64/t64
#! t64/t32
#! 
#! i32 <- sample(1e6)
#! d64 <- as.double(i32)
#! 
#! message("-- the following timings are rather conservative since timings
#!  of integer64 include garbage collection -- due to looped calls")
#! message("-- integer64 coerces 900x/100x faster than int64 
#!  under 64/32 bit OS (and 2x the time of coercing to integer) --")
#! t32 <- system.time(for(i in 1:1000)as.integer(d64))
#! t64 <- system.time(for(i in 1:1000)as.integer64(d64))
#! #T64 <- system.time(as.int64(d64))*1000
#! #T64/t64
#! t64/t32
#! td64 <- system.time(for(i in 1:1000)as.double(i32))
#! t64 <- system.time(for(i in 1:1000)as.integer64(i32))
#! #T64 <- system.time(for(i in 1:10)as.int64(i32))*100
#! #T64/t64
#! t64/td64
#! 
#! message("-- integer64 serializes 4x/0.8x faster than int64 
#!  under 64/32 bit OS (and less than 2x/6x the time of integer or double) --")
#! t32 <- system.time(for(i in 1:10)serialize(i32, NULL))
#! td64 <- system.time(for(i in 1:10)serialize(d64, NULL))
#! i64 <- as.integer64(i32); 
#! t64 <- system.time(for(i in 1:10)serialize(i64, NULL))
#! rm(i64); gc()
#! #I64 <- as.int64(i32); 
#! #T64 <- system.time(for(i in 1:10)serialize(I64, NULL))
#! #rm(I64); gc()
#! #T64/t64
#! t64/t32
#! t64/td64
#! 
#! 
#! message("-- integer64 adds 250x/60x faster than int64
#!  under 64/32 bit OS (and less than 6x the time of integer or double) --")
#! td64 <- system.time(for(i in 1:100)d64+d64)
#! t32 <- system.time(for(i in 1:100)i32+i32)
#! i64 <- as.integer64(i32); 
#! t64 <- system.time(for(i in 1:100)i64+i64)
#! rm(i64); gc()
#! #I64 <- as.int64(i32); 
#! #T64 <- system.time(for(i in 1:10)I64+I64)*10
#! #rm(I64); gc()
#! #T64/t64
#! t64/t32
#! t64/td64
#! 
#! message("-- integer64 sums 3x/0.2x faster than int64 
#! (and at about 5x/60X the time of integer and double) --")
#! td64 <- system.time(for(i in 1:100)sum(d64))
#! t32 <- system.time(for(i in 1:100)sum(i32))
#! i64 <- as.integer64(i32); 
#! t64 <- system.time(for(i in 1:100)sum(i64))
#! rm(i64); gc()
#! #I64 <- as.int64(i32); 
#! #T64 <- system.time(for(i in 1:100)sum(I64))
#! #rm(I64); gc()
#! #T64/t64
#! t64/t32
#! t64/td64
#! 
#! message("-- integer64 diffs 5x/0.85x faster than integer and double
#! (int64 version 1.0 does not support diff) --")
#! td64 <- system.time(for(i in 1:10)diff(d64, lag=2L, differences=2L))
#! t32 <- system.time(for(i in 1:10)diff(i32, lag=2L, differences=2L))
#! i64 <- as.integer64(i32); 
#! t64 <- system.time(for(i in 1:10)diff(i64, lag=2L, differences=2L))
#! rm(i64); gc()
#! t64/t32
#! t64/td64
#! 
#! 
#! message("-- integer64 subscripts 1000x/340x faster than int64
#! (and at the same speed / 10x slower as integer) --")
#! ts32 <- system.time(for(i in 1:1000)sample(1e6, 1e3))
#! t32<- system.time(for(i in 1:1000)i32[sample(1e6, 1e3)])
#! i64 <- as.integer64(i32); 
#! t64 <- system.time(for(i in 1:1000)i64[sample(1e6, 1e3)])
#! rm(i64); gc()
#! #I64 <- as.int64(i32); 
#! #T64 <- system.time(for(i in 1:100)I64[sample(1e6, 1e3)])*10
#! #rm(I64); gc()
#! #(T64-ts32)/(t64-ts32)
#! (t64-ts32)/(t32-ts32)
#! 
#! message("-- integer64 assigns 200x/90x faster than int64
#! (and 50x/160x slower than integer) --")
#! ts32 <- system.time(for(i in 1:100)sample(1e6, 1e3))
#! t32 <- system.time(for(i in 1:100)i32[sample(1e6, 1e3)] <- 1:1e3)
#! i64 <- as.integer64(i32); 
#! i64 <- system.time(for(i in 1:100)i64[sample(1e6, 1e3)] <- 1:1e3)
#! rm(i64); gc()
#! #I64 <- as.int64(i32); 
#! #I64 <- system.time(for(i in 1:10)I64[sample(1e6, 1e3)] <- 1:1e3)*10
#! #rm(I64); gc()
#! #(T64-ts32)/(t64-ts32)
#! (t64-ts32)/(t32-ts32)
#! 
#! 
#! tdfi32 <- system.time(dfi32 <- data.frame(a=i32, b=i32, c=i32))
#! tdfsi32 <- system.time(dfi32[1e6:1,])
#! fi32 <- tempfile()
#! tdfwi32 <- system.time(write.csv(dfi32, file=fi32, row.names=FALSE))
#! tdfri32 <- system.time(read.csv(fi32, colClasses=rep("integer", 3)))
#! unlink(fi32)
#! rm(dfi32); gc()
#! 
#! i64 <- as.integer64(i32); 
#! tdfi64 <- system.time(dfi64 <- data.frame(a=i64, b=i64, c=i64))
#! tdfsi64 <- system.time(dfi64[1e6:1,])
#! fi64 <- tempfile()
#! tdfwi64 <- system.time(write.csv(dfi64, file=fi64, row.names=FALSE))
#! tdfri64 <- system.time(read.csv(fi64, colClasses=rep("integer64", 3)))
#! unlink(fi64)
#! rm(i64, dfi64); gc()
#! 
#! #I64 <- as.int64(i32); 
#! #tdfI64 <- system.time(dfI64<-data.frame(a=I64, b=I64, c=I64))
#! #tdfsI64 <- system.time(dfI64[1e6:1,])
#! #fI64 <- tempfile()
#! #tdfwI64 <- system.time(write.csv(dfI64, file=fI64, row.names=FALSE))
#! #tdfrI64 <- system.time(read.csv(fI64, colClasses=rep("int64", 3)))
#! #unlink(fI64)
#! #rm(I64, dfI64); gc()
#! 
#! message("-- integer64 coerces 40x/6x faster to data.frame than int64
#! (and factor 1/9 slower than integer) --")
#! #tdfI64/tdfi64
#! tdfi64/tdfi32
#! message("-- integer64 subscripts from data.frame 20x/2.5x faster than int64
#!  (and 3x/13x slower than integer) --")
#! #tdfsI64/tdfsi64
#! tdfsi64/tdfsi32
#! message("-- integer64 csv writes about 2x/0.5x faster than int64
#! (and about 1.5x/5x slower than integer) --")
#! #tdfwI64/tdfwi64
#! tdfwi64/tdfwi32
#! message("-- integer64 csv reads about 3x/1.5 faster than int64
#! (and about 2x slower than integer) --")
#! #tdfrI64/tdfri64
#! tdfri64/tdfri32
#! 
#! rm(i32, d64); gc()
#! 
#! 
#! message("-- investigating the impact on garbage collection: --")
#! message("-- the fragmented structure of int64 messes up R's RAM --")
#! message("-- and slows down R's gargbage collection just by existing --")
#! 
#! td32 <- double(21)
#! td32[1] <- system.time(d64 <- double(1e7))[3]
#! for (i in 2:11)td32[i] <- system.time(gc(), gcFirst=FALSE)[3]
#! rm(d64)
#! for (i in 12:21)td32[i] <- system.time(gc(), gcFirst=FALSE)[3]
#! 
#! t64 <- double(21)
#! t64[1] <- system.time(i64 <- integer64(1e7))[3]
#! for (i in 2:11)t64[i] <- system.time(gc(), gcFirst=FALSE)[3]
#! rm(i64)
#! for (i in 12:21)t64[i] <- system.time(gc(), gcFirst=FALSE)[3]
#! 
#! #T64 <- double(21)
#! #T64[1] <- system.time(I64 <- int64(1e7))[3]
#! #for (i in 2:11)T64[i] <- system.time(gc(), gcFirst=FALSE)[3]
#! #rm(I64)
#! #for (i in 12:21)T64[i] <- system.time(gc(), gcFirst=FALSE)[3]
#! 
#! #matplot(1:21, cbind(td32, t64, T64), pch=c("d","i","I"), log="y")
#! matplot(1:21, cbind(td32, t64), pch=c("d","i"), log="y")
#!   }
#!
#! }
#! \name{identical.integer64}
#! \alias{identical.integer64}
#! \title{
#!    Identity function for class 'integer64'
#! }
#! \description{
#!   This will discover any deviation between objects containing integer64 vectors. 
#! }
#! \usage{
#!  identical.integer64(x, y, num.eq = FALSE, single.NA = FALSE
#! , attrib.as.set = TRUE, ignore.bytecode = TRUE)
#! }
#! \arguments{
#!   \item{x}{ atomic vector of class 'integer64' }
#!   \item{y}{ atomic vector of class 'integer64' }
#!   \item{num.eq}{ see \code{\link{identical}} }
#!   \item{single.NA}{ see \code{\link{identical}} }
#!   \item{attrib.as.set}{ see \code{\link{identical}} }
#!   \item{ignore.bytecode}{ see \code{\link{identical}} }
#! }
#! \details{
#!   This is simply a wrapper to \code{\link{identical}} with default arguments \code{num.eq = FALSE, single.NA = FALSE}.
#! }
#! \value{
#!   A single logical value, \code{TRUE} or \code{FALSE}, never \code{NA} and never anything other than a single value. 
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \keyword{ classes }
#! \keyword{ manip }
#! \seealso{ \code{\link{==.integer64}} \code{\link{identical}} \code{\link{integer64}}  }
#! \examples{
#!   i64 <- as.double(NA); class(i64) <- "integer64"
#!   identical(i64-1, i64+1)
#!   identical.integer64(i64-1, i64+1)
#! }


#! \name{as.character.integer64}
#! \alias{as.character.integer64}
#! \alias{as.double.integer64}
#! \alias{as.integer.integer64}
#! \alias{as.logical.integer64}
#! \alias{as.bitstring}
#! \alias{as.bitstring.integer64}
#! \alias{as.factor.integer64}
#! \alias{as.ordered.integer64}
#! \title{
#!    Coerce from integer64
#! }
#! \description{
#!   Methods to coerce integer64 to other atomic types. 
#!   'as.bitstring' coerces to a human-readable bit representation (strings of zeroes and ones). 
#!   The methods \code{\link{format}}, \code{\link{as.character}}, \code{\link{as.double}},
#!   \code{\link{as.logical}}, \code{\link{as.integer}} do what you would expect.
#! }
#! \usage{
#!  as.bitstring(x, \dots)
#!  \method{as.bitstring}{integer64}(x, \dots)
#!  \method{as.character}{integer64}(x, \dots)
#!  \method{as.double}{integer64}(x, \dots)
#!  \method{as.integer}{integer64}(x, \dots)
#!  \method{as.logical}{integer64}(x, \dots)
#!  \method{as.factor}{integer64}(x)
#!  \method{as.ordered}{integer64}(x)
#! }
#! \arguments{
#!   \item{x}{ an integer64 vector }
#!   \item{\dots}{ further arguments to the \code{\link{NextMethod}} }
#! }
#! \value{
#!   \code{as.bitstring} returns a string of . \cr
#!   The other methods return atomic vectors of the expected types
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \keyword{ classes }
#! \keyword{ manip }
#! \seealso{ \code{\link{as.integer64.character}} \code{\link{integer64}}  }
#! \examples{
#!   as.character(lim.integer64())
#!   as.bitstring(lim.integer64())
#! }

#! \name{as.integer64.character}
#! \alias{as.integer64}
#! \alias{as.integer64.integer64}
#! \alias{as.integer64.NULL}
#! \alias{as.integer64.character}
#! \alias{as.integer64.double}
#! \alias{as.integer64.integer}
#! \alias{as.integer64.logical}
#! \alias{as.integer64.factor}
#! \alias{NA_integer64_}
#! \title{
#!    Coerce to integer64
#! }
#! \description{
#!   Methods to coerce from other atomic types to integer64. 
#! }
#! \usage{
#!  NA_integer64_
#!  as.integer64(x, \dots)
#!  \method{as.integer64}{integer64}(x, \dots)
#!  \method{as.integer64}{NULL}(x, \dots)
#!  \method{as.integer64}{character}(x, \dots)
#!  \method{as.integer64}{double}(x, \dots)
#!  \method{as.integer64}{integer}(x, \dots)
#!  \method{as.integer64}{logical}(x, \dots)
#!  \method{as.integer64}{factor}(x, \dots)
#! }
#! \arguments{
#!   \item{x}{ an atomic vector }
#!   \item{\dots}{ further arguments to the \code{\link{NextMethod}} }
#! }
#! \details{
#!   \code{as.integer64.character} is realized using C function \code{strtoll} which does not support scientific notation. 
#!   Instead of '1e6' use '1000000'.
#! }
#! \value{
#!   The other methods return atomic vectors of the expected types
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \keyword{ classes }
#! \keyword{ manip }
#! \seealso{ \code{\link{as.character.integer64}} \code{\link{integer64}}  }
#! \examples{
#!   as.integer64(as.character(lim.integer64()))
#! }


#! \name{extract.replace.integer64}
#! \alias{[.integer64}
#! \alias{[[.integer64}
#! \alias{[[<-.integer64}
#! \alias{[<-.integer64}
#! \title{
#!    Extract or Replace Parts of an integer64 vector
#! }
#! \description{
#!   Methods to extract and replace parts of an integer64 vector.
#! }
#! \usage{
#!  \method{[}{integer64}(x, \dots)
#!  \method{[}{integer64}(x, \dots) <- value 
#!  \method{[[}{integer64}(x, \dots)
#!  \method{[[}{integer64}(x, \dots) <- value
#! }
#! \arguments{
#!   \item{x}{ an atomic vector }
#!   \item{value}{ an atomic vector with values to be assigned }
#!   \item{\dots}{ further arguments to the \code{\link{NextMethod}} }
#! }
#! \note{
#!   You should not subscript non-existing elements and not use \code{NA}s as subscripts.
#!   The current implementation returns \code{9218868437227407266} instead of \code{NA}.
#! }
#! \value{
#!   A vector or scalar of class 'integer64'
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \keyword{ classes }
#! \keyword{ manip }
#! \seealso{ \code{\link{[}} \code{\link{integer64}}  }
#! \examples{
#!   as.integer64(1:12)[1:3]
#!   x <- as.integer64(1:12)
#!   dim(x) <- c(3,4)
#!   x
#!   x[]
#!   x[,2:3]
#! }

#! \name{format.integer64}
#! \alias{format.integer64}
#! \alias{is.na.integer64}
#! \alias{!.integer64}
#! \alias{sign.integer64}
#! \alias{abs.integer64}
#! \alias{sqrt.integer64}
#! \alias{log.integer64}
#! \alias{log2.integer64}
#! \alias{log10.integer64}
#! \alias{floor.integer64}
#! \alias{ceiling.integer64}
#! \alias{trunc.integer64}
#! \alias{round.integer64}
#! \alias{signif.integer64}
#! \title{
#!    Unary operators and functions for integer64 vectors
#! }
#! \description{
#!   Unary operators and functions for integer64 vectors.
#! }
#! \usage{
#! \method{format}{integer64}(x, justify="right", \dots)
#! \method{is.na}{integer64}(x)
#! \method{!}{integer64}(x)
#! \method{sign}{integer64}(x)
#! \method{abs}{integer64}(x)
#! \method{sqrt}{integer64}(x)
#! \method{log}{integer64}(x, base)
#! \method{log2}{integer64}(x)
#! \method{log10}{integer64}(x)
#! \method{floor}{integer64}(x)
#! \method{ceiling}{integer64}(x)
#! \method{trunc}{integer64}(x, \dots)
#! \method{round}{integer64}(x, digits=0)
#! \method{signif}{integer64}(x, digits=6)
#! }
#! \arguments{
#!   \item{x}{ an atomic vector of class 'integer64'}
#!   \item{base}{ an atomic scalar (we save 50\% log-calls by not allowing a vector base) }
#!   \item{digits}{ integer indicating the number of decimal places (round) or significant digits (signif) to be used. 
#!                  Negative values are allowed (see \code{\link{round}}) }
#!   \item{justify}{ should it be right-justified (the default), left-justified, centred or left alone. }
#!   \item{\dots}{ further arguments to the \code{\link{NextMethod}} }
#! }
#! \value{
#!   \code{\link{format}} returns a character vector \cr
#!   \code{\link{is.na}} and \code{\link{!}} return a logical vector \cr
#!   \code{\link{sqrt}}, \code{\link{log}}, \code{\link{log2}} and \code{\link{log10}} return a double vector \cr
#!   \code{\link{sign}}, \code{\link{abs}}, \code{\link{floor}}, \code{\link{ceiling}}, \code{\link{trunc}} and 
#!   \code{\link{round}} return a vector of class 'integer64' \cr
#!   \code{\link{signif}} is not implemented 
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \keyword{ classes }
#! \keyword{ manip }
#! \seealso{ \code{\link{xor.integer64}} \code{\link{integer64}}  }
#! \examples{
#!   sqrt(as.integer64(1:12))
#! }


#! \name{xor.integer64}
#! \alias{&.integer64}
#! \alias{|.integer64}
#! \alias{xor.integer64}
#! \alias{!=.integer64}
#! \alias{==.integer64}
#! \alias{<.integer64}
#! \alias{<=.integer64}
#! \alias{>.integer64}
#! \alias{>=.integer64}
#! \alias{+.integer64}
#! \alias{-.integer64}
#! \alias{*.integer64}
#! \alias{^.integer64}
#! \alias{/.integer64}
#! \alias{\%/\%.integer64}
#! \alias{\%\%.integer64}
#! \alias{binattr}
#! \title{
#!    Binary operators for integer64 vectors
#! }
#! \description{
#!   Binary operators for integer64 vectors.
#! }
#! \usage{
#! \method{&}{integer64}(e1,e2)
#! \method{|}{integer64}(e1,e2)
#! \method{xor}{integer64}(x,y)
#! \method{!=}{integer64}(e1,e2)
#! \method{==}{integer64}(e1,e2)
#! \method{<}{integer64}(e1,e2)
#! \method{<=}{integer64}(e1,e2)
#! \method{>}{integer64}(e1,e2)
#! \method{>=}{integer64}(e1,e2)
#! \method{+}{integer64}(e1,e2)
#! \method{-}{integer64}(e1,e2)
#! \method{*}{integer64}(e1,e2)
#! \method{^}{integer64}(e1,e2)
#! \method{/}{integer64}(e1,e2)
#! \method{\%/\%}{integer64}(e1,e2)
#! \method{\%\%}{integer64}(e1,e2)
#! binattr(e1,e2) # for internal use only
#! }
#! \arguments{
#!   \item{e1}{ an atomic vector of class 'integer64'}
#!   \item{e2}{ an atomic vector of class 'integer64'}
#!   \item{x}{ an atomic vector of class 'integer64'}
#!   \item{y}{ an atomic vector of class 'integer64'}
#! }
#! \value{
#!   \code{\link{&}}, \code{\link{|}}, \code{\link{xor}}, \code{\link{!=}}, \code{\link{==}}, 
#!   \code{\link{<}}, \code{\link{<=}}, \code{\link{>}}, \code{\link{>=}} return a logical vector \cr
#!   \code{\link{^}} and \code{\link{/}} return a double vector\cr
#!   \code{\link{+}}, \code{\link{-}}, \code{\link{*}}, \code{\link{\%/\%}}, \code{\link{\%\%}}
#!    return a vector of class 'integer64'
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \keyword{ classes }
#! \keyword{ manip }
#! \seealso{ \code{\link{format.integer64}} \code{\link{integer64}}  }
#! \examples{
#!   as.integer64(1:12) - 1
#! }


#! \name{sum.integer64}
#! \alias{all.integer64}
#! \alias{any.integer64}
#! \alias{min.integer64}
#! \alias{max.integer64}
#! \alias{range.integer64}
#! \alias{lim.integer64}
#! \alias{sum.integer64}
#! \alias{prod.integer64}
#! \title{
#!    Summary functions for integer64 vectors
#! }
#! \description{
#!   Summary functions for integer64 vectors. 
#!   Function 'range' without arguments returns the smallest and largest value of the 'integer64' class.
#! }
#! \usage{
#! \method{all}{integer64}(\dots, na.rm = FALSE)
#! \method{any}{integer64}(\dots, na.rm = FALSE)
#! \method{min}{integer64}(\dots, na.rm = FALSE)
#! \method{max}{integer64}(\dots, na.rm = FALSE)
#! \method{range}{integer64}(\dots, na.rm = FALSE)
#! lim.integer64()
#! \method{sum}{integer64}(\dots, na.rm = FALSE)
#! \method{prod}{integer64}(\dots, na.rm = FALSE)
#! }
#! \arguments{
#!   \item{\dots}{ atomic vectors of class 'integer64'}
#!   \item{na.rm}{ logical scalar indicating whether to ignore NAs }
#! }
#! \details{
#!   The numerical summary methods always return \code{integer64}. 
#!   Therefor the methods for \code{min},\code{max} and \code{range} do not return \code{+Inf,-Inf}
#!   on empty arguments, but \code{+9223372036854775807, -9223372036854775807} (in this sequence).
#!   The same is true if only  \code{NA}s are submitted with argument \code{na.rm=TRUE}. 
#!  \cr
#!   \code{lim.integer64} returns these limits in proper order \code{-9223372036854775807, +9223372036854775807} and without a \code{\link{warning}}.
#! }
#! \value{
#!   \code{\link{all}} and \code{\link{any}} return a logical scalar\cr
#!   \code{\link{range}} returns a integer64 vector with two elements\cr
#!   \code{\link{min}}, \code{\link{max}}, \code{\link{sum}} and \code{\link{prod}} return a integer64 scalar
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \keyword{ classes }
#! \keyword{ manip }
#! \seealso{ \code{\link{mean.integer64}} \code{\link{cumsum.integer64}} \code{\link{integer64}}  }
#! \examples{
#!   lim.integer64()
#!   range(as.integer64(1:12))
#! }


#! \name{cumsum.integer64}
#! \alias{cummin.integer64}
#! \alias{cummax.integer64}
#! \alias{cumsum.integer64}
#! \alias{cumprod.integer64}
#! \alias{diff.integer64}
#! \title{
#!    Cumulative Sums, Products, Extremes and lagged differences
#! }
#! \description{
#!   Cumulative Sums, Products, Extremes and lagged differences
#! }
#! \usage{
#! \method{cummin}{integer64}(x)
#! \method{cummax}{integer64}(x)
#! \method{cumsum}{integer64}(x)
#! \method{cumprod}{integer64}(x)
#! \method{diff}{integer64}(x, lag = 1L, differences = 1L, \dots)
#! }
#! \arguments{
#!   \item{x}{ an atomic vector of class 'integer64'}
#!   \item{lag}{ see \code{\link{diff}} }
#!   \item{differences}{ see \code{\link{diff}} }
#!   \item{\dots}{ ignored }
#! }
#! \value{
#!   \code{\link{cummin}}, \code{\link{cummax}} , \code{\link{cumsum}} and \code{\link{cumprod}} 
#!      return a integer64 vector of the same length as their input\cr
#!   \code{\link{diff}} returns a integer64 vector shorter by \code{lag*differences} elements \cr
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \keyword{ classes }
#! \keyword{ manip }
#! \seealso{ \code{\link{sum.integer64}} \code{\link{integer64}}  }
#! \examples{
#!   cumsum(rep(as.integer64(1), 12))
#!   diff(as.integer64(c(0,1:12)))
#!   cumsum(as.integer64(c(0, 1:12)))
#!   diff(cumsum(as.integer64(c(0,0,1:12))), differences=2)
#! }


#! \name{c.integer64}
#! \alias{c.integer64}
#! \alias{cbind.integer64}
#! \alias{rbind.integer64}
#! \title{
#!    Concatenating integer64 vectors
#! }
#! \description{
#!   The ususal functions 'c', 'cbind' and 'rbind'
#! }
#! \usage{
#! \method{c}{integer64}(\dots, recursive = FALSE)
#! \method{cbind}{integer64}(\dots)
#! \method{rbind}{integer64}(\dots)
#! }
#! \arguments{
#!   \item{\dots}{ two or more arguments coerced to 'integer64' and passed to \code{\link{NextMethod}} }
#!   \item{recursive}{ logical. If \code{recursive = TRUE}, the function recursively descends through lists (and pairlists) combining all their elements into a vector. }
#! }
#! \value{
#!   \code{\link{c}} returns a integer64 vector of the total length of the input \cr
#!   \code{\link{cbind}} and \code{\link{rbind}} return a integer64 matrix
#! }
#! \note{
#!   R currently only dispatches generic 'c' to method 'c.integer64' if the first argument is 'integer64'
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \keyword{ classes }
#! \keyword{ manip }
#! \seealso{ \code{\link{rep.integer64}} \code{\link{seq.integer64}} 
#!           \code{\link{as.data.frame.integer64}} \code{\link{integer64}}  
#! }
#! \examples{
#!   c(as.integer64(1), 2:6)
#!   cbind(1:6, as.integer(1:6))
#!   rbind(1:6, as.integer(1:6))
#! }


#! \name{rep.integer64}
#! \alias{rep.integer64}
#! \title{
#!    Replicate elements of integer64 vectors
#! }
#! \description{
#!   Replicate elements of integer64 vectors
#! }
#! \usage{
#! \method{rep}{integer64}(x, \dots)
#! }
#! \arguments{
#!   \item{x}{ a vector of 'integer64' to be replicated }
#!   \item{\dots}{ further arguments passed to \code{\link{NextMethod}} }
#! }
#! \value{
#!   \code{\link{rep}} returns a integer64 vector
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \keyword{ classes }
#! \keyword{ manip }
#! \seealso{ \code{\link{c.integer64}} \code{\link{rep.integer64}} 
#!           \code{\link{as.data.frame.integer64}} \code{\link{integer64}}  
#! }
#! \examples{
#!   rep(as.integer64(1:2), 6)
#!   rep(as.integer64(1:2), c(6,6))
#!   rep(as.integer64(1:2), length.out=6)
#! }


#! \name{seq.integer64}
#! \alias{seq.integer64}
#! \title{
#!    integer64: Sequence Generation
#! }
#! \description{
#!   Generating sequence of integer64 values
#! }
#! \usage{
#! \method{seq}{integer64}(from = NULL, to = NULL, by = NULL, length.out = NULL, along.with = NULL, \dots)
#! }
#! \arguments{
#!   \item{from}{ integer64 scalar (in order to dispatch the integer64 method of \code{\link{seq}} }
#!   \item{to}{ scalar }
#!   \item{by}{ scalar }
#!   \item{length.out}{ scalar }
#!   \item{along.with}{ scalar }
#!   \item{\dots}{ ignored }
#! }
#! \details{
#!   \code{seq.integer64} does coerce its arguments 'from', 'to' and 'by' to \code{integer64}.
#!   If not provided, the argument 'by' is automatically determined as \code{+1} or \code{-1},
#!   but the size of 'by' is not calculated as in \code{\link{seq}} (because this might result in a non-integer value).
#! }
#! \value{
#!   an integer64 vector with the generated sequence
#! }
#! \note{
#!   In base R \code{\link{:}} currently is not generic and does not dispatch, see section "Limitations inherited from Base R" in \code{\link{integer64}}
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \keyword{ classes }
#! \keyword{ manip }
#! \seealso{ \code{\link{c.integer64}} \code{\link{rep.integer64}} 
#!           \code{\link{as.data.frame.integer64}} \code{\link{integer64}}  
#! }
#! \examples{
#!   # colon not activated: as.integer64(1):12
#!   seq(as.integer64(1), 12, 2)
#!   seq(as.integer64(1), by=2, length.out=6)
#! }


#! \name{as.data.frame.integer64}
#! \alias{as.data.frame.integer64}
#! \title{
#!    integer64: Coercing to data.frame column
#! }
#! \description{
#!   Coercing integer64 vector to data.frame.
#! }
#! \usage{
#!   \method{as.data.frame}{integer64}(x, \dots)
#! }
#! \arguments{
#!   \item{x}{ an integer64 vector }
#!   \item{\dots}{ passed to NextMethod \code{\link{as.data.frame}} after removing the 'integer64' class attribute }
#! }
#! \value{
#!   a one-column data.frame containing an integer64 vector
#! }
#! \details{
#!   'as.data.frame.integer64' is rather not intended to be called directly,
#!   but it is required to allow integer64 as data.frame columns.
#! }
#! \note{
#!   This is currently very slow -- any ideas for improvement?
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \keyword{ classes }
#! \keyword{ manip }
#! \seealso{ 
#!   \code{\link{cbind.integer64}} \code{\link{integer64}}  %as.vector.integer64 removed as requested by the CRAN maintainer \code{\link{as.vector.integer64}} 
#! }
#! \examples{
#!   as.data.frame.integer64(as.integer64(1:12))
#!   data.frame(a=1:12, b=as.integer64(1:12))
#! }



#! \name{plusclass}
#! \alias{plusclass}
#! \alias{minusclass}
#! \title{
#!    integer64: Maintaining S3 class attribute
#! }
#! \description{
#!   Maintaining integer64 S3 class attribute.
#! }
#! \usage{
#!   plusclass(class, whichclass)
#!   minusclass(class, whichclass)
#! }
#! \arguments{
#!   \item{class}{ NULL or a character vector of class attributes }
#!   \item{whichclass}{ the (single) class name to add or remove from the class vector  }
#! }
#! \value{
#!   NULL or a character vector of class attributes
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \keyword{ classes }
#! \keyword{ manip }
#! \keyword{ internal }
#! \seealso{ 
#!   \code{\link{oldClass}} \code{\link{integer64}}  
#! }
#! \examples{
#!   plusclass("inheritingclass","integer64")
#!   minusclass(c("inheritingclass","integer64"), "integer64")
#! }


# if (!exists(":.default")){
	# ":.default" <- get(":")
	# ":" <- function(from,to)UseMethod(":")
# }

setOldClass("integer64")


identical.integer64 <- function(x, y
, num.eq = FALSE
, single.NA = FALSE
, attrib.as.set = TRUE
, ignore.bytecode = TRUE
)
identical(x=x, y=y
, num.eq = num.eq
, single.NA = single.NA
, attrib.as.set = attrib.as.set
, ignore.bytecode = ignore.bytecode
)


as.integer64 <- function (x, ...) 
UseMethod("as.integer64")

as.bitstring <- function(x, ...)
UseMethod("as.bitstring")



minusclass <- function(class, whichclass){
  if (length(class)){
	  i <- whichclass==class
	  if (any(i))
		class[!i]
	  else
		class
  }else
    class
}

plusclass <- function(class, whichclass){
  if (length(class)){
	  i <- whichclass==class
	  if (any(i))
		class
	  else
		c(class, whichclass)
  }else
    whichclass
}

binattr <- function(e1,e2){
  d1 <- dim(e1)
  d2 <- dim(e2)
  n1 <- length(e1)
  n2 <- length(e2)
  if (length(d1)){
    if (length(d2)){
	  if (!identical(dim(e1),dim(e2)))
		stop("non-conformable arrays")
	}else{
	  if (n2>n1)
	    stop("length(e2) does not match dim(e1)")
	  if (n1%%n2)
		warning("length(e1) not a multiple length(e2)")
	}
	attributes(e1)
  }else{
    if (length(d2)){
	  if (n1>n2)
	    stop("length(e1) does not match dim(n2)")
	  if (n2%%n1)
		warning("length(e2) not a multiple length(e1)")
	  attributes(e2)
	}else{
	  if (n1<n2){
		if (n2%%n1)
			warning("length(e2) not a multiple length(e1)")
	  }else{
		if (n1%%n2)
			warning("length(e1) not a multiple length(e2)")
	  }
	  attributes(e1)
	}
  }
}


integer64 <- function(length=0){
  ret <- double(length)
  oldClass(ret) <- "integer64"
  ret
}

is.integer64 <- function(x)inherits(x, "integer64")

as.integer64.NULL <- function (x, ...){
  ret <- double()
  oldClass(ret) <- "integer64"
  ret
}

as.integer64.integer64 <- function(x, ...)x

as.integer64.double <- function(x, ...){
  ret <- double(length(x))
  .Call("as_integer64_double", x, ret)
  oldClass(ret) <- "integer64"
  ret
}

as.integer64.logical <- as.integer64.integer <- function(x, ...){
  ret <- double(length(x))
  .Call("as_integer64_integer", x, ret)
  oldClass(ret) <- "integer64"
  ret
}

as.integer64.character <- function(x, ...){
  n <- length(x)
  ret <- rep(as.double(NA), n)
  .Call("as_integer64_character", x, ret)
  oldClass(ret) <- "integer64"
  ret
}

as.integer64.factor <- function(x, ...)
as.integer64(unclass(x), ...)

as.double.integer64 <- function(x, ...){
  ret <- double(length(x))
  .Call("as_double_integer64", x, ret)
  ret
}

as.integer.integer64 <- function(x, ...){
  ret <- integer(length(x))
  .Call("as_integer_integer64", x, ret)
  ret
}

as.logical.integer64 <- function(x, ...){
  ret <- logical(length(x))
  .Call("as_logical_integer64", x, ret)
  ret
}

as.character.integer64 <- function(x, ...){
  n <- length(x)
  ret <- rep(as.character(NA), n)
  .Call("as_character_integer64", x, ret)
  ret
}

as.bitstring.integer64 <- function(x, ...){
  n <- length(x)
  ret <- rep(as.character(NA), n)
  .Call("as_bitstring_integer64", x, ret)
  ret
}

# read.table expects S4 as() 
setAs("character","integer64",function(from)as.integer64.character(from))
setAs("integer64","character",function(from)as.character.integer64(from))

# this is a trick to generate NA_integer64_ for namespace export before 
# as.integer64() is available because dll is not loaded
NA_integer64_ <- unserialize(as.raw(c(0x58, 0x0a, 0x00, 0x00, 0x00, 0x02, 0x00, 0x03, 0x03, 
+ 0x00, 0x00, 0x02, 0x03, 0x00, 0x00, 0x00, 0x03, 0x0e, 0x00, 0x00, 
+ 0x00, 0x01, 0x80, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
+ 0x00, 0x04, 0x02, 0x00, 0x00, 0x00, 0x01, 0x00, 0x04, 0x00, 0x09, 
+ 0x00, 0x00, 0x00, 0x05, 0x63, 0x6c, 0x61, 0x73, 0x73, 0x00, 0x00, 
+ 0x00, 0x10, 0x00, 0x00, 0x00, 0x01, 0x00, 0x04, 0x00, 0x09, 0x00, 
+ 0x00, 0x00, 0x09, 0x69, 0x6e, 0x74, 0x65, 0x67, 0x65, 0x72, 0x36, 
+ 0x34, 0x00, 0x00, 0x00, 0xfe)))

"length<-.integer64" <- function(x, value){
  cl <- oldClass(x)
  n <- length(x)
  x <- NextMethod()
  oldClass(x) <- cl
  if (value>n)
    x[(n+1):value] <- 0L
  x
}


format.integer64 <- function(x, justify="right", ...){
  a <- attributes(x)
  x <- as.character(x)
  ret <- format(x, justify=justify, ...)
  a$class <- minusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

print.integer64 <- function(x, quote=FALSE, ...){
  cat("integer64\n")
  a <- attributes(x)
  ret <- as.character(x)
  a$class <- minusclass(a$class, "integer64")
  attributes(ret) <- a
  print(ret, quote=quote, ...)
  invisible(x)
}

"[.integer64" <- function(x,...){
  cl <- oldClass(x)
  ret <- NextMethod()
  oldClass(ret) <- cl
  remcache(ret)
  ret
}

"[<-.integer64" <- function(x,...,value){
  cl <- oldClass(x)
  value <- as.integer64(value)
  ret <- NextMethod()
  oldClass(ret) <- cl
  ret
}

"[[.integer64" <- function(x,...){
  cl <- oldClass(x)
  ret <- NextMethod()
  oldClass(ret) <- cl
  ret
}

"[[<-.integer64" <- function(x,...,value){
  cl <- oldClass(x)
  value <- as.integer64(value)
  ret <- NextMethod()
  oldClass(ret) <- cl
  ret
}

c.integer64 <-
function (..., recursive = FALSE) 
{
	l <- list(...)
	K <- length(l)
	for (k in 1:K){
		if (recursive && is.list(l[[k]])){
			l[[k]] <- do.call("c.integer64", c(l[[k]], list(recursive = TRUE)))
		}else{
			if (!is.integer64(l[[k]])) {
				nam <- names(l[[k]])
				l[[k]] <- as.integer64(l[[k]])
				names(l[[k]]) <- nam
			}
			oldClass(l[[k]]) <- NULL
		}
	}
	ret <- do.call("c", l)
	oldClass(ret) <- "integer64"
	ret
}


cbind.integer64 <- function(...){
  l <- list(...)
	K <- length(l)
  for (k in 1:K){
		if (!is.integer64(l[[k]])){
			nam <- names(l[[k]])
			l[[k]] <- as.integer64(l[[k]])
			names(l[[k]]) <- nam
		}
		oldClass(l[[k]]) <- NULL
  }
  ret <- do.call("cbind", l)
	oldClass(ret) <- "integer64"
  ret
}

rbind.integer64 <- function(...){
  l <- list(...)
	K <- length(l)
  for (k in 1:K){
		if (!is.integer64(l[[k]])){
			nam <- names(l[[k]])
			l[[k]] <- as.integer64(l[[k]])
			names(l[[k]]) <- nam
		}
		oldClass(l[[k]]) <- NULL
  }
  ret <- do.call("rbind", l)
	oldClass(ret) <- "integer64"
  ret
}

# tenfold runtime if using attr() here instead of setattr()
# as.data.frame.integer64 <- function(x, ...){
  # cl <- oldClass(x)
  # oldClass(x) <- minusclass(cl, "integer64")
  # ret <- as.data.frame(x, ...)
  # k <- length(ret)
  # for (i in 1:k)
    # oldClass(ret[[i]]) <- cl
  # ret
# }
as.data.frame.integer64 <- function(x, ...){
  cl <- oldClass(x)
  on.exit(setattr(x, "class", cl))
  setattr(x, "class", minusclass(cl, "integer64"))
  ret <- as.data.frame(x, ...)
  k <- length(ret)
  for (i in 1:k)
   setattr(ret[[i]], "class", cl) 
  ret
}


"rep.integer64" <- function(x, ...){
	cl <- oldClass(x)
	ret <- NextMethod()
	oldClass(ret) <- cl
	ret
}

# FIXME no method dispatch for :
":.integer64" <- function(from, to){
  from <- as.integer64(from)
  to <- as.integer64(to)
  ret <- double(as.integer(to-from+1L))
  .Call("seq_integer64", from, as.integer64(1L), ret)
  oldClass(ret) <- "integer64"
  ret
}

"seq.integer64" <- function(from=NULL, to=NULL, by=NULL, length.out=NULL, along.with=NULL, ...){
    if (is.null(length.out))
      length.out <- length(along.with)
    else 
      length.out <- as.integer(length.out)
	
    if (is.null(by)){
      if (is.null(from) || is.null(to))
	    by <- as.integer64(1L)
	  else
	    by <- as.integer64(sign(to-from))
    }else{
	  by <- as.integer64(by)
	  if ((!is.null(from)) && (!is.null(to)) && sign(by)!=sign(to-from))
        stop("wrong sign of 'by' argument")
    }
  
    if (is.null(from)){
      if (length.out && length(to))
	    from <- to - (length.out-1L)*by
	  else
	    from <- as.integer64(1)
    }else 
	  from <- as.integer64(from)
	  
	if (!length(to)){
	  if (length.out)
	    to <- from + (length.out-1L)*by
      else
		stop("not enough informatoin provided")
	}
	
    if (!length.out){
	  length.out <- (to-from) %/% by + 1L
    }

    if (length.out){
      if (length.out==1L)
        return(from)
      else{
        #return(cumsum(c(from, rep(by, length.out-1L))))
		ret <- double(as.integer(length.out))
		.Call("seq_integer64", from, by, ret)
		oldClass(ret) <- "integer64"
		return(ret)
	  }
    }else
      return(integer64())
}


"+.integer64" <- function(e1, e2){
  if (missing(e2))
    return(e1)
  a <- binattr(e1,e2)
  e1 <- as.integer64(e1)
  e2 <- as.integer64(e2)
  ret <- double(max(length(e1),length(e2)))
  .Call("plus_integer64", e1, e2, ret)
  a$class <- plusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

"-.integer64" <- function(e1, e2){
  if (missing(e2)){
    e2 <- e1
	e1 <- 0L
  }
  a <- binattr(e1,e2)
  e1 <- as.integer64(e1)
  e2 <- as.integer64(e2)
  ret <- double(max(length(e1),length(e2)))
  .Call("minus_integer64", e1, e2, ret)
  a$class <- plusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

"%/%.integer64" <- function(e1, e2){
  a <- binattr(e1,e2)
  e1 <- as.integer64(e1)
  e2 <- as.integer64(e2)
  ret <- double(max(length(e1), length(e2)))
  .Call("intdiv_integer64", e1, e2, ret)
  a$class <- plusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

"%%.integer64" <- function(e1, e2){
  a <- binattr(e1,e2)
  e1 <- as.integer64(e1)
  e2 <- as.integer64(e2)
  ret <- double(max(length(e1), length(e2)))
  .Call("mod_integer64", e1, e2, ret)
  a$class <- plusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}


"*.integer64" <- function(e1, e2){
  a <- binattr(e1,e2)
  ret <- double(max(length(e1),length(e2)))
  if (!is.integer64(e2) && is.double(e2))
    .Call("times_integer64_double", as.integer64(e1), e2, ret)
  else
    .Call("times_integer64_integer64", as.integer64(e1), as.integer64(e2), ret)
  a$class <- plusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

"^.integer64" <- function(e1, e2){
  a <- binattr(e1,e2)
  ret <- double(max(length(e1),length(e2)))
  if (!is.integer64(e2) && is.double(e2))
    .Call("power_integer64_double", as.integer64(e1), e2, ret)
  else
    .Call("power_integer64_integer64", as.integer64(e1), as.integer64(e2), ret)
  a$class <- plusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

"/.integer64" <- function(e1, e2){
  a <- binattr(e1,e2)
  ret <- double(max(length(e1),length(e2)))
  if (!is.integer64(e2) && is.double(e2))
	  .Call("divide_integer64_double", as.integer64(e1), e2, ret)
  else
	  .Call("divide_integer64_integer64", as.integer64(e1), as.integer64(e2), ret)
  a$class <- minusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}


"sign.integer64" <- function(x){
  a <- attributes(x)
  ret <- double(length(x))
  .Call("sign_integer64", x,ret)
  attributes(ret) <- a
  ret
}

"abs.integer64" <- function(x){
  a <- attributes(x)
  ret <- double(length(x))
  .Call("abs_integer64", x,ret)
  attributes(ret) <- a
  ret
}

"sqrt.integer64" <- function(x){
  a <- attributes(x)
  ret <- double(length(x))
  .Call("sqrt_integer64", x,ret)
  a$class <- minusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

"log.integer64" <- function(x, base=NULL){
  a <- attributes(x)
  ret <- double(max(length(x),length(base)))
  if (is.null(base))
	.Call("log_integer64", x, ret)
  else if(length(base)==1){
    .Call("logbase_integer64", x, as.double(base), ret)
  }else{
    .Call("logvect_integer64", x, as.double(base), ret)
  }
  a$class <- minusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

"log10.integer64" <- function(x){
  a <- attributes(x)
  ret <- double(length(x))
  .Call("log10_integer64", x,ret)
  a$class <- minusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

"log2.integer64" <- function(x){
  a <- attributes(x)
  ret <- double(length(x))
  .Call("log2_integer64", x,ret)
  a$class <- minusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

"trunc.integer64" <- function(x, ...)x
"floor.integer64" <- "ceiling.integer64" <- function(x)x

"signif.integer64" <- function(x, digits=6)x

"round.integer64" <- function(x, digits=0){
  if (digits<0){
    a <- attributes(x)
	base <- 10^floor(-digits)
	ret <- (x%/%base) * base
    #a$class <- minusclass(a$class, "integer64")
    attributes(ret) <- a
	ret
  }else
	x
}

"any.integer64" <- function(..., na.rm = FALSE){
  l <- list(...)
  ret <- logical(1)
  if (length(l)==1){
		  .Call("any_integer64", l[[1]], na.rm, ret)
		  ret
  }else{
	  any(sapply(l, function(e){
		if (is.integer64(e)){
		  .Call("any_integer64", e, na.rm, ret)
		  ret
		}else{
		  any(e, na.rm = na.rm)
		}
	  }), na.rm = na.rm)
  }
}

"all.integer64" <- function(..., na.rm = FALSE){
  l <- list(...)
  ret <- logical(1)
  if (length(l)==1){
		  .Call("all_integer64", l[[1]], na.rm, ret)
		  ret
  }else{
	  all(sapply(l, function(e){
		if (is.integer64(e)){
		  .Call("all_integer64", e, na.rm, ret)
		  ret
		}else{
		  all(e, na.rm = na.rm)
		}
	  }), na.rm = na.rm)
  }
}

"sum.integer64" <- function(..., na.rm = FALSE){
  l <- list(...)
  ret <- double(1)
  if (length(l)==1){
		  .Call("sum_integer64", l[[1]], na.rm, ret)
		  oldClass(ret) <- "integer64"
		  ret
  }else{
	  ret <- sapply(l, function(e){
		if (is.integer64(e)){
		  .Call("sum_integer64", e, na.rm, ret)
		  ret
		}else{
		  as.integer64(sum(e, na.rm = na.rm))
		}
	  })
    oldClass(ret) <- "integer64"
	  sum(ret, na.rm = na.rm)
  }
}



"prod.integer64" <- function(..., na.rm = FALSE){
  l <- list(...)
  ret <- double(1)
  if (length(l)==1){
		  .Call("prod_integer64", l[[1]], na.rm, ret)
		  oldClass(ret) <- "integer64"
		  ret
  }else{
      ret <- sapply(l, function(e){
		if (is.integer64(e)){
		  .Call("prod_integer64", e, na.rm, ret)
		  ret
		}else{
		  as.integer64(prod(e, na.rm = na.rm))
		}
	  })
	  oldClass(ret) <- "integer64"
	  prod(ret, na.rm = na.rm)
  }
}

"min.integer64" <- function(..., na.rm = FALSE){
  l <- list(...)
  ret <- double(1)
  noval <- TRUE
  if (length(l)==1){
	if (length(l[[1]]))
	  noval <- FALSE
    .Call("min_integer64", l[[1]], na.rm, ret)
    oldClass(ret) <- "integer64"
  }else{
	  ret <- sapply(l, function(e){
	    if (length(e))
	      noval <<- FALSE
		if (is.integer64(e)){
		  .Call("min_integer64", e, na.rm, ret)
		  ret
		}else{
		  as.integer64(min(e, na.rm = na.rm))
		}
	  })
	  oldClass(ret) <- "integer64"
	  ret <- min(ret, na.rm = na.rm)
  }
  if (noval)
	warning("no non-NA value, returning +9223372036854775807")
  ret
}

"max.integer64" <- function(..., na.rm = FALSE){
  l <- list(...)
  ret <- double(1)
  noval <- TRUE
  if (length(l)==1){
	if (length(l[[1]]))
	  noval <- FALSE
	.Call("max_integer64", l[[1]], na.rm, ret)
	oldClass(ret) <- "integer64"
  }else{
	ret <- sapply(l, function(e){
	    if (length(e))
	      noval <<- FALSE
		if (is.integer64(e)){
		  .Call("max_integer64", e, na.rm, ret)
		  ret
		}else{
		  as.integer64(max(e, na.rm = na.rm))
		}
	})
	oldClass(ret) <- "integer64"
	ret <- max(ret, na.rm = na.rm)
  }
  if (noval)
	warning("no non-NA value, returning -9223372036854775807")
  ret
}


"range.integer64" <- function(..., na.rm = FALSE){
  ret <- double(2)
  l <- list(...)
  noval <- TRUE
  if (length(l)==1){
	if (length(l[[1]]))
	  noval <- FALSE
	.Call("range_integer64", l[[1]], na.rm, ret)
	oldClass(ret) <- "integer64"
  }else{
      ret <- unlist(sapply(l, function(e){
	    if (length(e))
	      noval <<- FALSE
		if (is.integer64(e)){
		  .Call("range_integer64", e, na.rm, ret)
		  ret
		}else{
		  as.integer64(range(e, na.rm = na.rm))
		}
	  }))
	  oldClass(ret) <- "integer64"
	  ret <- range(ret, na.rm = na.rm)
  }
  if (noval)
	warning("no non-NA value, returning c(+9223372036854775807, -9223372036854775807)")
  ret
}

lim.integer64 <- function(){
    ret <- double(2)
	.Call("lim_integer64", ret)
	oldClass(ret) <- "integer64"
	return(ret)
}

"diff.integer64" <- function(x, lag=1L, differences=1L, ...){
  lag <- as.integer(lag)
  n <- length(x)
  d <- differences <- as.integer(differences)
  while(d>0L){
	n <- n - lag
    if (n<=0L){
	  ret <- double()
	  break
	}
	if (d==differences){
	  ret <- double(n)
      .Call("diff_integer64", x, as.integer64(lag), as.integer64(n), ret)
	}else{
	  .Call("diff_integer64", ret, as.integer64(lag), as.integer64(n), ret)
	}
	d <- d - 1L
  }
  length(ret) <- n
  oldClass(ret) <- "integer64"
  ret
}


"cummin.integer64" <- function(x){
  ret <- double(length(x))
  .Call("cummin_integer64", x,ret)
  oldClass(ret) <- "integer64"
  ret
}
"cummax.integer64" <- function(x){

  ret <- double(length(x))
  .Call("cummax_integer64", x,ret)
  oldClass(ret) <- "integer64"
  ret
}

"cumsum.integer64" <- function(x){
  ret <- double(length(x))
  .Call("cumsum_integer64", x,ret)
  oldClass(ret) <- "integer64"
  ret
}

"cumprod.integer64" <- function(x){
  ret <- double(length(x))
  .Call("cumprod_integer64", x,ret)
  oldClass(ret) <- "integer64"
  ret
}



"is.na.integer64" <- function(x){
  a <- attributes(x)
  ret <- logical(length(x))
  .Call("isna_integer64", x, ret)
  a$class <- minusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

"==.integer64" <- function(e1, e2){
  a <- binattr(e1,e2)
  e1 <- as.integer64(e1)
  e2 <- as.integer64(e2)
  ret <- logical(max(length(e1), length(e2)))
  .Call("EQ_integer64", e1, e2, ret)
  a$class <- minusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

"!=.integer64" <- function(e1, e2){
  a <- binattr(e1,e2)
  e1 <- as.integer64(e1)
  e2 <- as.integer64(e2)
  ret <- logical(max(length(e1), length(e2)))
  .Call("NE_integer64", e1, e2, ret)
  a$class <- minusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

"<.integer64" <- function(e1, e2){
  a <- binattr(e1,e2)
  e1 <- as.integer64(e1)
  e2 <- as.integer64(e2)
  ret <- logical(max(length(e1), length(e2)))
  .Call("LT_integer64", e1, e2, ret)
  a$class <- minusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

"<=.integer64" <- function(e1, e2){
  a <- binattr(e1,e2)
  e1 <- as.integer64(e1)
  e2 <- as.integer64(e2)
  ret <- logical(max(length(e1), length(e2)))
  .Call("LE_integer64", e1, e2, ret)
  a$class <- minusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

">.integer64" <- function(e1, e2){
  a <- binattr(e1,e2)
  e1 <- as.integer64(e1)
  e2 <- as.integer64(e2)
  ret <- logical(max(length(e1), length(e2)))
  .Call("GT_integer64", e1, e2, ret)
  a$class <- minusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

">=.integer64" <- function(e1, e2){
  a <- binattr(e1,e2)
  e1 <- as.integer64(e1)
  e2 <- as.integer64(e2)
  ret <- logical(max(length(e1), length(e2)))
  .Call("GE_integer64", e1, e2, ret)
  a$class <- minusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

"&.integer64" <- function(e1, e2){
  a <- binattr(e1,e2)
  ret <- as.logical(e1) & as.logical(e2)
  a$class <- minusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

"|.integer64" <- function(e1, e2){
  a <- binattr(e1,e2)
  ret <- as.logical(e1) | as.logical(e2)
  a$class <- minusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

xor.integer64 <- function(x, y){
  a <- binattr(x,y)
  ret <- as.logical(x) != as.logical(y)
  a$class <- minusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}


"!.integer64" <- function(x){
  a <- attributes(x)
  ret <- !as.logical(x)
  a$class <- minusclass(a$class, "integer64")
  attributes(ret) <- a
  ret
}

# as.vector.integer64 removed as requested by the CRAN maintainer
# as.vector.integer64 <- function(x, mode="any"){
  # ret <- NextMethod()
  # if (mode=="any")
	# oldClass(ret) <- "integer64"
  # ret
# }

# bug in R does not dispatch
is.vector.integer64 <- function(x, mode="any"){
  cl <- minusclass(oldClass(x), "integer64")
  a <- attributes(x)
  a$class <- NULL
  a$names <- NULL
  if (is.na(match(mode, c("any","integer64"))) || length(cl) || length(a) )
    FALSE
  else
    TRUE
}

