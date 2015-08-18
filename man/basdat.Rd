\name{basdat}
\Rdversion{1.1}
\alias{basdat}
\docType{data}
\title{HIV: TB and Survival (Baseline Data)}

\description{Simulated dataset. Baseline data of 386 HIV positive individuals, including time of first active tuberculosis, time of death, individual end time. Time varying CD4 measurements of these patients are included in dataset \code{\link{timedat}}.}

\usage{data(basdat)}
\format{
A data frame with 386 observations on the following 4 variables.
\describe{
\item{\code{id}}{patient ID.}
\item{\code{Ttb}}{time of first active tuberculosis, measured in days since HIV seroconversion.}
\item{\code{Tdeath}}{time of death, measured in days since HIV seroconversion.}
\item{\code{Tend}}{individual end time (either death or censoring), measured in days since HIV seroconversion.}
}
}

\details{
These simulated data are used together with data in \code{\link{timedat}} in a detailed causal modelling example using inverse probability weighting (IPW). See \code{\link{ipwtm}} for the example. Data were simulated using the algorithm described in Van der Wal e.a. (2009).
}

\references{
Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse  Probability Weighting. \emph{Journal of Statistical Software}, \bold{43}(13), 1-23. \url{http://www.jstatsoft.org/v43/i13/}.

Van der Wal W.M., Prins M., Lumbreras B. & Geskus R.B. (2009). A simple G-computation algorithm to quantify the causal effect of a secondary illness on the progression of a chronic disease. \emph{Statistics in Medicine}, \bold{28}(18), 2325-2337.
}

\author{Willem M. van der Wal \email{w.m.vanderwal@amc.uva.nl}}

\seealso{\code{\link{basdat}}, \code{\link{haartdat}}, \code{\link{ipwplot}}, \code{\link{ipwpoint}}, \code{\link{ipwtm}}, \code{\link{timedat}}, \code{\link{tstartfun}}.}

\examples{
#see ?ipwtm for example
}
\keyword{datasets}
