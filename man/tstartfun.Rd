\name{tstartfun}
\Rdversion{1.1}
\alias{tstartfun}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Compute Starting Time For Counting Process Notation}

\description{Function to compute starting time for intervals of follow-up, when using the counting process notation. Within each unit under observation (usually individuals), computes starting time equal to:
\itemize{
\item time of previous record when there is a previous record.
\item 0 for first record, when time of first record is after 0.
\item -1 for first record, when time of first record is 0.
}}


\usage{
tstartfun(id, timevar, data)
}

\arguments{
\item{id}{numerical vector, uniquely identifying the units under observation, within which the longitudinal measurements are taken.}
\item{timevar}{numerical vector, representing follow-up time, starting at 0.}
\item{data}{dataframe containing \code{id} and \code{timevar}.}
}


\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{Numerical vector containing starting time for each record. In the same order as the records in \code{data}, to facilitate merging.}

\section{Missing values}{\code{id} and \code{timevar} should not contain missing values.}

\references{
%% ~put references to the literature/web site here ~
}

\author{Willem M. van der Wal \email{w.m.vanderwal@amc.uva.nl}}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{basdat}}, \code{\link{haartdat}}, \code{\link{healthdat}}, \code{\link{ipwplot}}, \code{\link{ipwpoint}}, \code{\link{ipwtm}}, \code{\link{timedat}}, \code{\link{tstartfun}}.}

\examples{
#data
mydata1 <- data.frame(
   patient = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2),
   time.days = c(14, 34, 41, 56, 72, 98, 0, 11, 28, 35))

#compute starting time for each interval
mydata1$tstart <- tstartfun(patient, time.days, mydata1)

#result
mydata1

#see also ?ipwtm for example
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{survival}

