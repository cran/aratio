\name{lwr1_minh}
\alias{lwr1_minh}
\title{Locally Weighted Regression bandwidth or window selection}
\description{
Called by lwr1 to determine which value of potential bandwidths or window sizes minimizes a criterion function.  
}
\usage{
lwr1_minh(form, window = 0, bandwidth = 0, kern = "tcub", 
  method = "gcv", print = TRUE, alldata = FALSE, data=NULL)

}
\arguments{
  \item{form}{
A formula of the type y~x.
}
  \item{window}{Vector of possible bandwidths.  Default:  none.
}
  \item{bandwidth}{
Vector of possible bandwidths.  Default:  none.
}
  \item{kern}{
Kernel weighting functions.  Default is the tri-cube. Options include  "rect", "tria", "epan", "bisq", "tcub", "trwt", and "gauss".
}
  \item{method}{ The selection criterion when searching for the optimal value for the bandwidth or window.  Must be in quotes.  
The default is the generalized cross-validation criterion, or "gcv".  The alternative is the cross-validation criterion, or "cv".  
See the description of the lwr1 function for more information.}

  \item{print}{If print=T, prints summary of estimation results for each potential window or bandwidth size.  Default:  print=T.
}

  \item{alldata}{If alldata=T, each observation is used as a target value for x.  
When alldata=F, the function is estimated at a set of points chosen by the locfit program using an adaptive decision tree approach, 
and the akima program is used to interpolate to the full set of observations.  Specifying alldata=T can lead to long estimation times.
}
  \item{data }{A data frame containing the data. Default:  use data in the current working directory.}
}

\value{
Returns the value of the chosen bandwidth or window
}

\author{
Daniel McMillen
}

\examples{
	set.seed(23849103)
	n = 1000
	x <- runif(n,0,2*pi)
	x <- sort(x)
	ybase <- x - .1*(x^2) + sin(x) - cos(x) -.5*sin(2*x) + .5*cos(2*x)
	sig = sd(ybase)/2
	y <- ybase + rnorm(n,0,sig)
	bandwidth <- lwr1_minh(y~x, bandwidth =  seq(from=.10, to=.90, by=.10))
	window <- lwr1_minh(y~x, window=c(.10, .15, .20, .25), method="cv", kern="epan")

}
