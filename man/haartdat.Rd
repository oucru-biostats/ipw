\name{haartdat}
\Rdversion{1.1}
\alias{haartdat}
\docType{data}

\title{HAART and Survival in HIV Patients}

\description{
Survival data measured in 1200 HIV positive patients. Start of follow-up is HIV seroconversion. Each row corresponds to a 100 day interval of follow-up time, using the counting process notation. Patients can initiate HAART therapy. CD4 count is a confounder for the effect of HAART on mortality.
}

\usage{data(haartdat)}
\format{
\describe{
\item{\code{patient}}{patient ID}
\item{\code{tstart}}{starting time for each interval of follow-up, measured in days since HIV seroconversion}
\item{\code{fuptime}}{end time for each interval of follow-up, measured in days since HIV seroconversion}
\item{\code{haartind}}{indicator for the initiation of HAART therapy at the end of the interval (0=HAART not initiated/1=HAART initiated).}
\item{\code{event}}{indicator for death at the end of the interval (0=alive/1=died)}
\item{\code{sex}}{sex (0=male/1=female)}
\item{\code{age}}{age at the start of follow-up (years)}
\item{\code{cd4.sqrt}}{square root of CD4 count, measured at fuptime, before \code{haartind}}
}
}
\details{These data were simulated.

Patients can initiate HAART at \code{fuptime=0}. Therefore, to allow the fitting of a model predicting initiation of HAART, starting time for the first interval within each patient is negative (-100).
}

\references{
Van der Wal W.M. & Geskus R.B. (2011). ipw: An R Package for Inverse  Probability Weighting. \emph{Journal of Statistical Software}, \bold{43}(13), 1-23. \url{http://www.jstatsoft.org/v43/i13/}.
}

\author{Willem M. van der Wal \email{w.m.vanderwal@amc.uva.nl}}

\seealso{\code{\link{basdat}}, \code{\link{haartdat}}, \code{\link{ipwplot}}, \code{\link{ipwpoint}}, \code{\link{ipwtm}}, \code{\link{timedat}}, \code{\link{tstartfun}}.}

\examples{
#see ?ipwtm for example
}

\keyword{datasets}
