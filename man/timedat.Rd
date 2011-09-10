\name{timedat}
\Rdversion{1.1}
\alias{timedat}
\docType{data}
\title{HIV: TB and Survival (Longitudinal Measurements)}

\description{Simulated dataset. Time varying CD4 measurements of 386 HIV positive individuals. Time of first active tuberculosis, time of death and individual end time of the patients are included in dataset \code{\link{basdat}}.}

\usage{data(timedat)}
\format{
	A data frame with 6291 observations on the following 3 variables.
	\describe{
		\item{\code{id}}{patient ID.}
		\item{\code{fuptime}}{follow-up time (days since HIV seroconversion).}
		\item{\code{cd4count}}{CD4 count measured at fuptime.}
	}
}
\details{
These simulated data are used together with data in \code{\link{basdat}} in a detailed causal modelling example using inverse probability weighting (IPW). See \code{\link{ipwtm}} for the example. Data were simulated using the algorithm described in Van der Wal e.a. (2009).
}

\references{
Cole, S.R. & Hernán, M.A. (2008). Constructing inverse probability weights for marginal structural models. \emph{American Journal of Epidemiology}, \bold{168}(6), 656-664.

Robins, J.M., Hernán, M.A. & Brumback, B.A. (2000). Marginal structural models and causal inference in epidemiology. \emph{Epidemiology}, \bold{11}, 550-560.

Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse  Probability Weighting. \emph{Journal of Statistical Software}, \bold{43}(13), 1-23. \url{http://www.jstatsoft.org/v43/i13/}.

Van der Wal W.M., Prins M., Lumbreras B. & Geskus R.B. (2009). A simple G-computation algorithm to quantify the causal effect of a secondary illness on the progression of a chronic disease. \emph{Statistics in Medicine}, \bold{28}(18), 2325-2337.
}

\author{Willem M. van der Wal \email{w.m.vanderwal@amc.uva.nl}}

\seealso{\code{\link{basdat}}, \code{\link{haartdat}}, \code{\link{ipwplot}}, \code{\link{ipwpoint}}, \code{\link{ipwtm}}, \code{\link{timedat}}, \code{\link{tstartfun}}.}

\examples{
#See ?ipwtm for example
}
\keyword{datasets}
