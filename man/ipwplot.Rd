\name{ipwplot}
\Rdversion{1.1}
\alias{ipwplot}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Plot Inverse Probability Weights
}
\description{
For time varying weights: display boxplots within strata of follow-up time.
For point treatment weights: display density plot.
}

\usage{
ipwplot(weights, timevar = NULL, binwidth = NULL, logscale = TRUE,
        xlab = NULL, ylab = NULL, main = "", ref = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{weights}{numerical vector of inverse probability weights to plot.}
\item{timevar}{numerical vector representing follow-up time. When specified, boxplots within strata of follow-up time are displayed. When left unspecified, a density plot is displayed.}
\item{binwidth}{numerical value indicating the width of the intervals of follow-up time; for each interval a boxplot is made. Ignored when \code{timevar} is not specified.}
\item{logscale}{logical value. If \code{TRUE}, weights are plotted on a logarithmic scale.}
\item{xlab}{label for the horizontal axis.}
\item{ylab}{label for the vertical axis.}
\item{main}{main title for the plot.}
\item{ref}{logical value. If \code{TRUE}, a reference line is plotted at \code{y=1}.}
\item{\dots}{additional arguments passed to \code{\link{boxplot}} (when \code{timevar} is specified) or \code{\link{plot}} (when \code{timevar} is not specified).}
}

\value{
A plot is displayed.
}


\author{Willem M. van der Wal \email{w.m.vanderwal@amc.uva.nl}}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{basdat}}, \code{\link{haartdat}}, \code{\link{ipwplot}}, \code{\link{ipwpoint}}, \code{\link{ipwtm}}, \code{\link{timedat}}, \code{\link{tstartfun}}.}

\examples{
#see ?ipwpoint and ?ipwtm for examples
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
