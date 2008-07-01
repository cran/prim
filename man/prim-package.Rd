\name{prim}
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
  %These regions are also known as the highest density regions (HDR).
  
  PRIM is a bump-hunting technique introduced by Friedman \& Fisher
  (1999), taken from data mining.
  PRIM estimates are a sequence of nested hyper-rectangles (boxes).
  %The
  %boxes which exceed a threshold comprise the HDR estimate. 

  For an overview of this package, see \code{vignette("prim")} for PRIM
  estimation for 2- and 5-dimensional data. 
}

\author{
Tarn Duong <tduong@maths.unsw.edu.au>
}

\references{
  %Duong, T., Koch, I., \& Wand, M.P. (2007) A generalised chi-squared
  %test for comparing sample from data rich sources: an example from flow
  %cytometry. \emph{In preparation.}
  Friedman, J.H. \& Fisher, N.I. (1999) Bump-hunting for high
  dimensional data, \emph{Statistics and Computing}, \bold{9}, 123--143.

  Hyndman, R.J.  Computing and graphing highest density
  regions. \emph{The American Statistician}, \bold{50}, 120--126. 
}


\keyword{ package }

%\seealso{}

