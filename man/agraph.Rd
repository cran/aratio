\name{agraph}
\alias{agraph}
\title{Histograms with density plots for assessment ratio data}

\description{agraph creates a histogram of assessment ratios created from assessed values and sales prices.  
Options allow kernel density plots and a normal distribution to be added to the diagrams.  
Lines can be added to the diagrams for statutory assessment rates and user-inputted 95\% confidence intervals for statistics such as the mean or median.}

\usage{
agraph(av, price, cftitle = "95\% CI", title = "Assessment Ratios", 
  width = 0.01, legloc = "none", freq = FALSE, normdens = FALSE, 
  kdens = FALSE, cfint = c(FALSE, FALSE), statute = FALSE, file = "none")
}

\arguments{
  \item{av}{Variable representing assessed value.}
  \item{price}{Variable representing price.}
  \item{cftitle}{Title for cfint in legend if legloc is not equal to "none" and cfint is not equal to c(F,F).  
Defaults to cfint="95\% CI" if cfint option is specified.
}
  \item{title}{Main title for histogram plot.  Default is title="Assessment Ratios"}
  \item{width}{Bin width for histogram.  Default is width=.01.  Smaller values show more detail.}
  \item{legloc}{Legend location.  Default is "none".  Other standard R options are available:  
"bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", and "center."
}
  \item{freq}{Indicates whether y-axis should be frequency (freq=T) or percent (freq=F).  Default is freq=F.}
  \item{normdens}{If normdens=T, normal density function is added to the plot.  Default is normdens=F.}
  \item{kdens}{If kerndens=T, kernel density function is added to the plot.  Default is kerndens=F.}
  \item{cfint}{Designed to add lines for a confidence interval to the plot. Values must be entered in a vector with two entries, e.g., cfint=c(.40, .60).  
Can be used to add a line for a single value by entering the same value twice, e.g., cfint=c(.5, .5).  Useful in combination with rstats function.  
Default is cfint=c(F,F), which means that no lines are added to the plot.
}
  \item{statute}{Adds a single vertical line to the plot.  Default is no line:  statute=F.  
Designed to add a line representing the statutory assessment rate, e.g, statute=.33
}
  \item{file}{Results are written to a file rather than the screen if file is specified.  Default is file="none".  
Example of a possible file location:  file="/dir1/dir2/graph.png"
}
}

\value{
agraph  is designed to produce a plot that is written to a file or produced on the screen.  The plot can be stored for later use, e.g., aplot <- agraph(av,price).}
\author{
Daniel McMillen
}



\examples{
data(dupage99)
dupage99$ratio <- dupage99$av/dupage99$price
dropobs1 <- nptrim_obs(dupage99$ratio,k=3)
dropobs2 <- (dupage99$price<quantile(dupage99$price,.01))|
  	(dupage99$price>quantile(dupage99$price,.99))
dupage99 <- dupage99[dropobs1==FALSE & dropobs2==FALSE,]
o <- order(dupage99$price)
dupage99 <- dupage99[o,]
attach(dupage99)

# add kernel density and normal density plots; smaller bin width for histotgram
agraph(av, price,normdens=TRUE, kdens=TRUE, width=.001)

# combine with ratio_stats to estimated 95% confidence interval assessment ratio mean.  
# vary title of the plot.  Add legend.  Frequency rather than percents for histogram.
# indicate statutory assessment ratio
rstat <- ratio_stats(av, price, nboot=1000)
lo = rstat[2,5]
hi = rstat[2,6]
agraph(av, price, normdens=TRUE, kdens=TRUE, legloc= "topleft", cfint=c(lo,hi), 
  freq=TRUE, title= "Histogram of Assessment Ratios", statute=.33)

# plot mean only; vary title for cfint
amean = rstat[1,1]
agraph(av, price, cfint=c(amean,amean), legloc="topright", 
	cftitle= "Mean", kdens=TRUE, statute=.33)
}


