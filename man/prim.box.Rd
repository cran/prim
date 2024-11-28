\name{prim.box}
\alias{prim.box}
\alias{prim.hdr}
\alias{prim.combine}

\title{PRIM for multivariate data}
\description{
  PRIM for multivariate data.
}
\usage{
prim.box(x, y, box.init=NULL, peel.alpha=0.05, paste.alpha=0.01,
     mass.min=0.05, threshold, pasting=TRUE, verbose=FALSE,
     threshold.type=0, y.fun=mean)

prim.hdr(prim, threshold, threshold.type, y.fun=mean)
prim.combine(prim1, prim2, y.fun=mean)
}

\arguments{
  \item{x}{matrix of data values}
  \item{y}{vector of response values}
  \item{y.fun}{function applied to response y. Default is mean.}
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
  \code{paste.alpha} and \code{mass.min} are taken from Friedman &
  Fisher (1999).

  The type of PRIM estimate is controlled \code{threshold} and
  \code{threshold.type}:
  \describe{
    \item{\code{threshold.type=1}}{search for \{\eqn{m(\bold{x}) \geq}{m(x) >=} \code{threshold}\}.}
    
    \item{\code{threshold.type=-1}}{search for \{\eqn{m(\bold{x}) \leq}{m(x) <=} \code{threshold}\}.}
      
    \item{\code{threshold.type=0}}{search for both \{\eqn{m(\bold{x}) \geq}{m(x) >=} \code{threshold[1]}\} and \{\eqn{m(\bold{x}) \leq}{m(x) <=} \code{threshold[2]}\}.}
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


\examples{
data(quasiflow)
qf <- quasiflow[1:1000,1:2]
qf.label <- quasiflow[1:1000,4]

## using only one command
thr <- c(0.25, -0.3)
qf.prim1 <- prim.box(x=qf, y=qf.label, threshold=thr, threshold.type=0)

## alternative - requires more commands but allows more control
## in intermediate stages
qf.primp <- prim.box(x=qf, y=qf.label, threshold.type=1)
   ## default threshold too low, try higher one

qf.primp.hdr <- prim.hdr(prim=qf.primp, threshold=0.25, threshold.type=1)
qf.primn <- prim.box(x=qf, y=qf.label, threshold=-0.3, threshold.type=-1)
qf.prim2 <- prim.combine(qf.primp.hdr, qf.primn)

plot(qf.prim1, alpha=0.2)   ## orange=x1>x2, blue x2<x1
points(qf[qf.label==1,], cex=0.5)
points(qf[qf.label==-1,], cex=0.5, col=2)
}
\keyword{multivariate}
