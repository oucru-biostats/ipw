\name{healthdat}
\Rdversion{1.1}
\alias{healthdat}
\docType{data}
\title{IQ, Income and Health}

\description{A simulated dataset containing IQ, income and health score measurements in 1000 individuals.}

\usage{data(healthdat)}

\format{
A data frame with 1000 rows, with each row corresponding to a separate individual. The following variables are included:
\describe{
\item{\code{id}}{individual ID.}
\item{\code{iq}}{IQ score.}
\item{\code{income}}{gross monthly income (EUR).}
\item{\code{health}}{health score (0-100).}
}}

\details{
In these simulated data, IQ is a confounder for the effect of income on health. 
}

\references{
Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse  Probability Weighting. \emph{Journal of Statistical Software}, \bold{43}(13), 1-23. \url{http://www.jstatsoft.org/v43/i13/}.
}

\author{Willem M. van der Wal \email{w.m.vanderwal@amc.uva.nl}}

\seealso{\code{\link{basdat}}, \code{\link{haartdat}}, \code{\link{healthdat}}, \code{\link{ipwplot}}, \code{\link{ipwpoint}}, \code{\link{ipwtm}}, \code{\link{timedat}}, \code{\link{tstartfun}}.}

\examples{
#see ?ipwpoint for example
}

\keyword{datasets}
