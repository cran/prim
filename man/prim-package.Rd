\name{prim-package}
\alias{prim} 
\docType{package}

\title{
Patient Rule Induction Method (PRIM)
}

\description{
  PRIM for bump-hunting for high-dimensional regression-type data.   
}

\details{
  The data are
  \eqn{(\bold{X}_1, Y_1), \dots, (\bold{X}_n, Y_n)}{(X_1, Y_1), \ldots,
    (X_n, Y_n)}  where \eqn{\bold{X}_i}{X_i} is d-dimensional and \eqn{Y_i}{Y_i} is a
  scalar response. We wish to find the modal (and/or anti-modal) regions
  in the conditional 
  expectation \eqn{ m(\bold{x}) = \bold{E} (Y | \bold{x}).}{m(x) = E(Y | x).}  
  
  PRIM is a bump-hunting technique introduced by Friedman & Fisher
  (1999), taken from data mining. PRIM estimates are a sequence of
  nested hyper-rectangles (boxes). 
 
  For an overview of this package, see \code{vignette("prim")} for PRIM
  estimation for 2- and 5-dimensional data. 
}

\author{
Tarn Duong <tarn.duong@gmail.com>
}

\references{
  Friedman, J.H. & Fisher, N.I. (1999) Bump-hunting for high
  dimensional data, \emph{Statistics and Computing}, \bold{9}, 123--143.

  Hyndman, R.J.  Computing and graphing highest density
  regions. \emph{American Statistician}, \bold{50}, 120--126. 
}


\keyword{ package }

%\seealso{}

