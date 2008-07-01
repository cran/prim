\name{prim.box}
\alias{prim.box}
\alias{prim.hdr}
\alias{prim.combine}


\title{PRIM for multivariate data}
\description{
  PRIM for multivariate data.
  %Result is an estimate of the highest density region (HDR).
}
\usage{
prim.box(x, y, box.init=NULL, peel.alpha=0.05, paste.alpha=0.01,
     mass.min=0.05, threshold, pasting=TRUE, verbose=FALSE,
     threshold.type=0)

prim.hdr(prim, threshold, threshold.type)
prim.combine(prim1, prim2)
}

\arguments{
  \item{x}{matrix of data values}
  \item{y}{vector of response values}
  \item{box.init}{initial covering box}
  \item{peel.alpha}{peeling quantile tuning parameter}
  \item{paste.alpha}{pasting quantile tuning parameter}
  \item{mass.min}{minimum mass tuning parameter}
  \item{threshold}{threshold tuning parameter(s)} 
  \item{threshold.type}{threshold direction indicator: 1 = ">=
    threshold", -1 = "<= threshold", 0 = ">= threshold[1] & <= threshold[2]"}
  \item{pasting}{flag for pasting}
  \item{verbose}{flag for printing output during execution}
  \item{prim,prim1,prim2}{objects of type \code{prim}}
}
  


\details{
  The data are \eqn{(\bold{X}_1, Y_1), \dots, (\bold{X}_n, Y_n)}{(X_1, Y_1), \ldots,
    (X_n, Y_n)}  where \eqn{\bold{X}_i}{X_i} is d-dimensional and \eqn{Y_i}{Y_i} is a
  scalar response. PRIM finds modal (and/or anti-modal) regions in the
  conditional expectation \eqn{m(\bold{x}) = \bold{E} (Y | \bold{x}).}{m(x) = E(Y | x).} 
  %These regions are also called the highest density regions (HDR). 

  In general, \eqn{Y_i}{Y_i} can be real-valued. See
  \code{vignette("prim")}.
  Here, we focus on the special case for binary \eqn{Y_i}{Y_i}. Let
  \eqn{Y_i}{Y_i} = 1 when 
  \eqn{\bold{X}_i \sim F^+}{X_i ~ F+}; and \eqn{Y_i}{Y_i} = -1 when
  \eqn{\bold{X}_i \sim 
    F^-}{X_i ~ F-} where \eqn{F^+}{F+} and \eqn{F^-}{F-} are different
  distribution functions. In this set-up, PRIM finds the
  regions where \eqn{F^+}{F+} and \eqn{F^-}{F-} are most different.

  The tuning parameters \code{peel.alpha} and \code{paste.alpha} control
  the `patience' of PRIM. Smaller values involve more patience. Larger
  values less patience. The peeling steps remove data from a box till
  either the box mean is smaller than \code{threshold} or the box mass
  is less than \code{mass.min}.  Pasting is optional, and is used to correct any
  possible over-peeling. The default values for \code{peel.alpha},
  \code{paste.alpha} and \code{mass.min} are taken from Friedman \&
  Fisher (1999).

  The type of PRIM estimate is controlled \code{threshold} and
  \code{threshold.type}:
  \itemize{
    \item{}{For \code{threshold.type=1}, we search for \{\eqn{m(\bold{x}) \geq}{m(x) >=} \code{threshold}\}.}
    
    \item{}{For \code{threshold.type=-1}, we search for \{\eqn{m(\bold{x}) \leq}{m(x) <=} \code{threshold}\}.}
      
    \item{}{For \code{threshold.type=0}, we search for both \{\eqn{m(\bold{x}) \geq}{m(x) >=} \code{threshold[1]}\} and \{\eqn{m(\bold{x}) \leq}{m(x) <=} \code{threshold[2]}\}.}
  }

  There are two ways of using PRIM. One is \code{prim.box} with
  pre-specified threshold(s). This is appropriate when the threshold(s)
  are known to produce good estimates. 

  On the other hand, if the user doesn't provide threshold values then
  \code{prim.box} computes box sequences which cover the data
  range. These can then be pruned at a later stage. \code{prim.hdr}
  allows the user to specify many different threshold values in an
  efficient manner, without having to recomputing the entire PRIM box
  sequence. \code{prim.combine} can be used to join the regions computed
  from \code{prim.hdr}. See the examples below. 
}

\value{
  -- \code{prim.box} produces a PRIM estimate, an object of
  type \code{prim}, which is a list with 8 fields:
  
  \item{x}{list of data matrices}
  \item{y}{list of response variable vectors}
  \item{y.mean}{list of vectors of box mean for y}
  \item{box}{list of matrices of box limits (first row = minima,
    second row = maxima)}
  \item{mass}{vector of box masses (proportion of points inside
    a box)}
  \item{num.class}{total number of PRIM boxes}
  \item{num.hdr.class}{total number of PRIM boxes which form the HDR}
  \item{ind}{threshold direction indicator: 1 = ">= threshold", -1 = "<=threshold"}

  The above lists have \code{num.class} fields, one for each box.

  -- \code{prim.hdr} takes a \code{prim} object and prunes it using 
  different threshold values. Returns another \code{prim} object. This
  is much faster for experimenting with different threshold values than
  calling \code{prim.box} each time. 

  -- \code{prim.combine} combines two \code{prim} objects into a single
  prim object. Usually used in conjunction with \code{prim.hdr}. See examples below.
}

\references{
  Friedman, J.H. \& Fisher, N.I. (1999) Bump-hunting for high
  dimensional data, \emph{Statistics and Computing}, \bold{9}, 123--143. 
} 

%\seealso{}

\examples{
n <- 1000
set.seed(88192)

mus.p <- rbind(c(0,0), c(2,0), c(1, 2), c(2.5, 2))
Sigmas.p <- 0.125*rbind(diag(2), diag(c(0.5, 0.5)),
   diag(c(0.125, 0.25)), diag(c(0.125, 0.25))) 
props.p <- c(0.5, 0.25, 0.125, 0.125)

mus.n <- rbind(c(0,0), c(2,0), c(2.5, 2))
Sigmas.n <- 0.125*rbind(matrix(c(1,-0.6,-0.6,1), nrow=2),
   diag(c(0.5, 0.5)),diag(c(0.125, 0.25))) 
props.n <- c(0.625, 0.25, 0.125)

x.p <- rmvnorm.mixt(n, mus.p, Sigmas.p, props.p)
x.n <- rmvnorm.mixt(n, mus.n, Sigmas.n, props.n)
x <- rbind(x.p, x.n)
y <- c(rep(1, nrow(x.p)), rep(-1, nrow(x.n)))
  ## 1 = positive sample, -1 = negative sample

y.thr <- c(1, -0.35)

## using only one command

x.prim1 <- prim.box(x=x, y=y, threshold=y.thr, threshold.type=0)

## alternative - requires more commands but allows more control
## in intermediate stages

x.prim.hdr.p <- prim.box(x=x, y=y, threshold.type=1,
   threshold=1)

x.prim.n <- prim.box(x=x, y=y, threshold.type=-1)
summary(x.prim.n)
   ## threshold too high, try lower one

x.prim.hdr.n <- prim.hdr(x.prim.n, threshold=-0.35,
   threshold.type=-1)
x.prim2 <- prim.combine(x.prim.hdr.p, x.prim.hdr.n)
       
plot(x.prim2)

summary(x.prim1)
summary(x.prim2) ## should be exactly the same as command above
}
\keyword{multivariate}
