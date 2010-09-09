nptrim_obs <-
function(x,k=3,print=TRUE,data=NULL) {
  attach(data,warn.conflicts=FALSE)

  missobs <- is.na(x)

  n = length(x)
  nmiss = sum(missobs)
  qmin = min(x,na.rm=TRUE)
  q25 = quantile(x,.25,na.rm=TRUE)
  q75 = quantile(x,.75,na.rm=TRUE)
  qmax = max(x,na.rm=TRUE)
  lo = q25 - k*(q75-q25)
  hi = q75 + k*(q75-q25)
  trimobs <- ifelse(!is.na(x)&(x<lo|x>hi),1,0)
  ntrim = sum(trimobs)
  dropobs <- is.na(x)|trimobs==1

  trim_vect <- c(n,nmiss,n-nmiss,ntrim,sum(dropobs),n-nmiss-ntrim,qmin,q25,q75,qmax,lo,hi)
  names(trim_vect) <- c("Total number of observations","Number of missing observations","Number of non-missing observations",
     "Number of non-missing observatons trimmed","Total number of observations dropped","Number of non-missing observations after trimming",
     "minimum","25th percentile","75th percentile","maximum","lower bound for trimming","upper bound for trimming")
  x.out <- data.frame(trim_vect)
  colnames(x.out) <- " "
  if (print==TRUE) {print(x.out)}
  names(trim_vect) <- NULL
  return(dropobs)
}

