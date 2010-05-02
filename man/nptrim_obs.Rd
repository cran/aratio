\name{nptrim_obs}
\alias{nptrim_obs}
\title{Nonparametric outlier trimming procedure
}
\description{
A nonparametric trimming procedure for identifying outlier observations in a variable.
}
\usage{
nptrim_obs(x, k = 3, type=6, print = TRUE)

}
\arguments{
  \item{x}{Variable to be trimmed of outliers.}
  \item{k}{Factor determining range for trimming. Default: k=3}
  \item{type}{An integer between 1 and 9 selecting one of the nine quantile algorithms used by R's \emph{quantile} command. 
The algorithms are listed in the R's documentation for \emph{quantile}.  Default:  \emph{type=6}.}
  \item{print}{If print=T, presents summary information on x and trimming results.}
}
\value{
An n-vector indicating which observations are identified as outliers.
}
\references{
International Association of Assessing Officers, "Standard on Ratio Studies", \emph{Assessment Journal} (September/October 1999), 12-64.  Section 6.6 and Table 1.
}
\author{
Daniel McMillen
}

\examples{
data(dupage99)
dupage99$ratio <- dupage99$av/dupage99$price
dropobs <- nptrim_obs(dupage99$ratio,k=3)
dupage99 <- dupage99[dropobs==FALSE,]
}

\details{ 	

Interquartile range = quantile(x, .75) - quantile(x, .25) = q75 - q25 \cr
	Lower bound for trimming:  lo = q25 - k*(q75 - q25) \cr
 	Upper bound for trimming:  hi = q75 + k*(q75 - q25) \cr
	Output variable equals one if x<lo or x>hi. }

