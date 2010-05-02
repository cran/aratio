\name{ratio_stats}
\alias{ratio_stats}
\title{Assessment ratio study statistics and confidence intervals.}
\description{
Produces a table of standard assessment ratio study statistics with bootstrap standard errors and confidence intervals.
}
\usage{ratio_stats(av, price, nboot = 0, type=6, print = TRUE)}
\arguments{
  \item{av}{Assessed value.}
  \item{price}{Price.}
  \item{nboot}{Number of replications of bootstrap resampling procedure.  Default:  nboot = 0.}
  \item{type}{An integer between 1 and 9 selecting one of the nine quantile algorithms used by R's \emph{quantile} command. 
The algorithms are listed in the R's documentation for \emph{quantile}.  Default:  \emph{type=6}.}
  \item{print}{If print = T, prints the table of results.  Default print=T.}
}
\value{
If nboot=0, produces a vector holding the sample mean, median, value weighted mean, 
coefficient of dispersion (cod), value-weighted cod, and price-related differential (prd).

If nboot>0, produces a 6x6 matrix with the following form, with the cell entries indicating the row, column combinations:

\tabular{lcccccc}{
\tab Sample Statistc \tab Bootstrap Std err \tab lo1 \tab up1 \tab lo2 \tab up2 \cr
Mean \tab 1,1 \tab 1,2 \tab 1,3 \tab 1,4 \tab 1,5 \tab 1,6 \cr
Median \tab 2,1 \tab 2,2 \tab 2,3 \tab 2,4 \tab 2,5 \tab 2,6 \cr
Value-weighted mean \tab 3,1 \tab 3,2 \tab 3,3 \tab 3,4 \tab 3,5 \tab 3,6 \cr
cod \tab 4,1 \tab 4,2 \tab 4,3 \tab 4,4 \tab 4,5 \tab 4,6 \cr
Value-weighted cod  \tab 5,1 \tab 5,2 \tab 5,3 \tab 5,4 \tab 5,5 \tab 5,6 \cr
prd  \tab 6,1 \tab 6,2 \tab 6,3 \tab 6,4 \tab 6,5 \tab 6,6 \cr }

}
\references{
Gloudemans, Robert.  \emph{Mass Appraisal of Real Property.}  Chicago:  International Association of Assessing Officers, 1999.  Chapter 5.

International Association of Assessing Officers, "Standard on Ratio Studies", \emph{Assessment Journal} (September/October 1999), 12-64.

}

\author{
Daniel McMillen
}



\examples{
data(dupage99)
dupage99$ratio <- dupage99$av/dupage99$price
dropobs1 <- nptrim_obs(dupage99$ratio,k=3)
dropobs2 <- (dupage99$price<quantile(dupage99$price,.01))|
  	(dupage99$price>quantile(dupage99$price,.99))
dupage99 <- dupage99[dropobs1==FALSE&dropobs2==FALSE,]
o <- order(dupage99$price)
dupage99 <- dupage99[o,]
attach(dupage99)

set.seed(484717)
rstat <- ratio_stats(av, price, nboot=1000)
lo = rstat[2,5]
hi = rstat[2,6]
agraph(av, price, normdens=TRUE, kdens=TRUE, legloc= "topleft", 
	cfint=c(lo,hi), freq=TRUE, title= "Histogram of Assessment Ratios", 
	statute=.33)
}


\details{
Statistics for ratio <- av/price. Denote av by A and price by P in the following. \cr
1. Mean \cr
2. Median \cr
3. Value-weighted mean: \eqn{\sum_i \frac{(P_i*Ratio_i)}{\sum_i (P_i)} = \frac{mean(A)}{mean(P)}}{sum_i (P_i*Ratio_i)/sum_i (P_i) = mean(A)/mean(P)}\cr
4. Coefficient of dispersion: \eqn{cod=100*\sum_i \frac{|ratio_i - median(ratio)|}{n*median(ratio)}}{cod=100*sum(|ratio_i - median(ratio)|)/n*median(ratio)} \cr
5. Value-weighted coefficient of dispersion: \eqn{100*\sum_i\frac{P_i*|ratio_i-median(ratio)|}{median(ratio)*\sum_i P_i}}{100*sum(P_i*|ratio_i-median(ratio)|)/(median(ratio)*sum(P_i)} \cr
6. Price-related differential: mean/value weighted mean. \cr 

The package boot is used to calculate standard errors for all six statistics.  
The procedure involves drawing with replacement nboot new samples of size n from the original set of av, 
price pairs and re-calculating the statistics for each of the nboot samples.  
The standard error that is reported in the second column of the table is simply the standard deviation of the 
nboot values of the statistic in question, and the bootstrap 95 percent confidence interval is the original sample statistic 
plus or minus 1.96 times the standard error.  The lower and upper bounds of this version of the 95 percent confidence interval 
are reported in columns 3 and 4 of the table and are labeled "lo1" and "hi1".  
This confidence interval is most appropriate if the statistic follows a normal distribution.  
If the statistic does not follow a symmetric distribution, the "percentile method" will produce a more accurate confidence interval.  
Let B denote the vector containing the nboot re-calculations of the test statistics; the percentile method 95 percent confidence interval 
is simply lo2 = quantile(b, .025) and hi2 = quantile(b, .975).

}
