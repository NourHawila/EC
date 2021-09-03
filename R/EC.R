#' Restricted MLE of Variance
#'
#' This function calculates the maximum likelihood estimator (MLE) of the variance
#' restricted to P.T-P.C=-Delta.
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param Delta  numeric between -1 and 1 representing the constraint for P.C-P.T
#' @return list of length 2 where the first element is MR.T representing the MLE for the treatment group and the second element is MR.C representing the MLE for the control group
#' @examples
#' #calculates MLE of variance of 5/10 and 7/10 under the constraint P.T=P.C-0.2
#' restricted_ml(x.T=5,x.C=7,N.T=10,N.C=10,Delta=0.2)
#' @references
#' \insertRef{Miettinen:85}{EC}
#' @export
restricted_ml = function(x.T, x.C, N.T, N.C, Delta) {
  R.T = x.T/N.T
  R.C = x.C/N.C
  N = N.C+N.T
  C = x.C+x.T
  L3 = N
  L2 = (N.C+2*N.T) * Delta - N - C
  L1 = (N.T*Delta-N-2*x.T)*Delta+C
  L0 = x.T*Delta*(1-Delta)
  Q = (L2^3)/((3*L3)^3) - (L1*L2) / (6*(L3^2))+ L0/(2*L3)

  if (Q > -(1e-21) & Q < 0)	{Q =  Q - 1e-12}
  if (Q < (1e-21) & Q >= 0)	{Q = Q + 1e-12}
  P = (Q > 0)*(sqrt((L2^2) / ((3*L3)^2) - L1/(3*L3))) - (Q < 0) * (sqrt((L2^2)/((3*L3)^2) - L1/(3*L3)))

  TEMP = Q/(P^3)
  A = (1/3) * (pi+acos(pmin(pmax(TEMP, -1.0), 1.0)))
  MR.T = 2*P*cos(A) - L2/(3*L3)
  MR.C = MR.T + Delta

  my_list = list(MR.T = MR.T, MR.C = MR.C)
  return(my_list)
}



#' Delta-projected Test Statistic for Non-inferiority Hypothesis
#'
#' This function  calculate statistics for testing non-inferiority based on Santner & Snell (1980), Blackwelder (1982),
#' Miettinen & Nurminen (1985) and Farrington & Manning (1990)
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param Delta0  numeric between 0 and 1 representing the non-inferiority margin
#' @param method character representing the method for ordering criterion("MN","FM","SS","Blackwelder")
#' @return numeric representing the value of the test statistic
#' @examples
#' #First example calculates the Miettinen & Nurminen test statistic
#' #for the Rodary et al. (1989) study with proportion of success
#' #in the treatment group being 83/88 and 69/76 for the control
#' #with a 10% noninferiority margin.
#' #Second example calculates the Farrington & Manning test statistic
#' #for the Fries et al. (1993) study with proportion of success
#' #in the treatment group being 8/15 and 3/15 for the control
#' #with a 20% noninferiority margin.
#' orderfun(x.T=83,x.C=69,N.T=88,N.C=76,Delta0=0.1,method="MN")
#' orderfun(x.T=8,x.C=3,N.T=15,N.C=15,Delta0=0.2,method="FM")
#' @references
#' \insertRef{Santner:80}{EC}
#'
#' \insertRef{Blackwelder:82}{EC}
#'
#' \insertRef{Miettinen:85}{EC}
#'
#' \insertRef{Farrington:90}{EC}
#' @export
orderfun = function(x.T, x.C, N.T, N.C, Delta0, method) {
  rml = restricted_ml(
    x.T = x.T,
    x.C = x.C,
    N.T = N.T,
    N.C = N.C,
    Delta = Delta0
  )

  N = N.T+N.C
  R.T = x.T/N.T
  R.C = x.C/N.C

  if (method == "MN") {
    MR.T = rml$MR.T
    MR.C = rml$MR.C
    DEN = (MR.T*(1-MR.T)/N.T + MR.C*(1-MR.C)/N.C)

    mn.res = ifelse(DEN<=0, 0, (R.T-R.C+Delta0) / sqrt(DEN))
    out = mn.res
  }

  if (method == "FM") {
    MR.T = rml$MR.T
    MR.C = rml$MR.C
    DEN <- (MR.T*(1-MR.T)/N.T + MR.C*(1-MR.C)/N.C)

    fm.res = ifelse(DEN<=0, 0, sqrt(N/(N-1))*((R.T-R.C+Delta0) / sqrt(DEN)))
    out = fm.res
  }

  if (method == "SS") {
    ss.res = R.T-R.C
    out = ss.res
  }


  if (method == "Blackwelder") {
    DEN = (R.T*(1-R.T)/N.T + R.C*(1-R.C)/N.C)
    blackwelder.res = (R.T-R.C)/sqrt(DEN)
    out = blackwelder.res
  }

  out
}



#' Order Statistic for all 2x2 Tables
#'
#' This function calculate the order statistics for all 2x2 tables
#'
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param Delta0  numeric between 0 and 1 representing the non-inferiority margin
#' @param method character representing the method for ordering criterion("MN","FM","SS","Blackwelder")
#' @return array of dimensions (N.T+1)x(N.C+1) where the (i,j) element is the order statistic for x.T=i and x.C=j
#' @examples
#' #16x16 array of the the Santner & Snell test statistic for all
#' #possible 2x2 tables arising from a sample size of 15 in the
#' #treatment and control group and 30% noninferiority margin.
#' order_mat(N.T=15,N.C=15,Delta0=0.3,method="SS")
#' @references
#' \insertRef{Santner:80}{EC}
#' @export
order_mat = function(N.T, N.C, Delta0, method) {
  res = array(NA, dim = c(N.T + 1, N.C + 1))

  for (i1 in 0:(N.T)) {
    for (i2 in 0:(N.C)) {
      res[i1 + 1, i2 + 1] = orderfun(
        x.T = i1,
        x.C = i2,
        N.T = N.T,
        N.C = N.C,
        Delta0 = Delta0,
        method = method
      )
    }
  }
  res
}

#' Verifying Barnard's Criterion
#'
#' This function checks whether Barnard's Criterion (Biometrika, 1947) is satisfied
#'
#' @param mat array of dimensions (N.T+1)x(N.C+1) where the (i,j) element is the order statistic for x.T=i and x.C=j
#' @return logical. TRUE if Barnard's Criterion is satisfied, FALSE otherwise
#' @examples
#' #Checks if the Barnard criterion is satisfied for the
#' #Miettenin & Nurminen ordering statistic with N.T=15,
#' #N.C=15 and Delta0=30%
#' barnard_check(order_mat(N.T=15,N.C=15,Delta0=0.3,method="MN"))
#' @export
#' @references
#' \insertRef{Barnard:47}{EC}
barnard_check=function(mat){
  N.T=dim(mat)[1]-1
  N.C=dim(mat)[2]-1

  out=TRUE
  for (i in 1:(N.T+1)) {
    for (j in 1:(N.C+1)) {

      if((i+1)<=(N.T+1)) {out=out*(mat[i,j]<=mat[i+1,j])}
      if((j+1)<=(N.C+1)) {out=out*(mat[i,j]>=mat[i,j+1])}
    }}
  if(out==1) return(TRUE)
  if(out==0) return(FALSE)
}

#' Chan's Exact P-value
#'
#' This function calculates the exact p-value based on Chan (1998)
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param Delta0  numeric between 0 and 1 representing the non-inferiority margin
#' @param method character representing the method for ordering criterion("MN","FM","SS","Blackwelder")
#' @param lower logical. TRUE for the null hypothesis P.T-P.C<=-Delta0.FALSE for the null hypothesis P.T-P.C>-Delta0
#' @param tol positive numeric representing the increment size for domain of Delta. Default is set to 0.001.
#' @return numeric representing Chan's exact p-value
#' @examples
#' #The first example is taken from Rodary et al. (1989) which was analyzed by Chan (1998)
#' #The second example is taken from Hawila (2021)
#' chan_pval(x.T=83,x.C=69,N.T=88,N.C=76,Delta0=0.1,method="MN")
#' chan_pval(x.T=2,x.C=0,N.T=15,N.C=10,Delta0=0.12,method="MN")
#' @references
#' \insertRef{Chan:98}{EC}
#' @export
chan_pval <- function(x.T, x.C, N.T, N.C, Delta0, method, lower = TRUE,tol=1e-3) {
  myorder_mat = order_mat(N.T, N.C, Delta0 = Delta0, method = method)
  obs = myorder_mat[x.T + 1, x.C + 1]

  P.T.vec = seq(max(0,-Delta0), min(1, 1 - Delta0), by = tol)
  pvals = rep(NA, length(P.T.vec))
  for (ptvar in 1:length(P.T.vec)) {
    P.T = P.T.vec[ptvar]
    P.C = P.T + Delta0

    if(lower==TRUE){
      pvals[ptvar] = sum((myorder_mat>=obs) * dbinom(0:N.T, N.T, P.T) %o%dbinom(0:N.C, N.C, P.C))
    }
    if(lower==FALSE){
      pvals[ptvar] = sum((myorder_mat<=obs) * dbinom(0:N.T, N.T, P.T) %o%dbinom(0:N.C, N.C, P.C))

    }
  }
  max(pvals)

}


#' CZ p-value
#'
#' This function computes the Chan & Zhang p-value given in Hawila (2021)
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param Delta0  numeric between 0 and 1 representing the non-inferiority margin
#' @param method character representing the method for ordering criterion("MN","FM","SS","Blackwelder")
#' @param tol positive numeric representing the tolerance for calculation
#' @return numeric representing the Chan & Zhang p-value
#' @examples
#' #This example is taken from Hawila (2021)
#' cz_pval(x.T=2,x.C=0,N.T=15,N.C=10,Delta0=0.12,method="MN")
#' @references
#' \insertRef{Hawila:21}{EC}
#' @export
cz_pval <- function(x.T, x.C, N.T, N.C, Delta0, method,tol=1e-3){
  deltas=seq(-0.999,-Delta0,by=tol)
  res=array(dim=length(deltas))
  for(d in 1:length(deltas)){
    res[d] = chan_pval(x.T,x.C,N.T,N.C,-deltas[d],method)
  }
  max(res)
}

#' Chan p-value for All 2x2 Tables
#'
#' This function calculates the Chan p-values for all 2x2 tables given N.T, N.C and an ordering criterion
#'
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param Delta0  numeric between 0 and 1 representing the non-inferiority margin
#' @param method character representing the method for ordering criterion("MN","FM","SS","Blackwelder")
#' @return array of dimension (N.T+1)x(N.C+1)  where the (i,j) element is the exact p-value for x.T=i and x.C=j
#' @examples
#'#8x6 array of the the Chan p-values based on Farrington & Manning
#'#statistic for all possible 2x2 tables arising from a sample size
#'#of 7 in the treatment group and 5 in the control group and 10%
#'#noninferiority margin.
#'chan_mat(N.T=7,N.C=5,Delta0=0.1,method="FM")
#'@references
#' \insertRef{Chan:98}{EC}
#' @export
chan_mat = function(N.T, N.C, Delta0, method) {
  mat = array(NA, dim = c(N.T + 1, N.C + 1))
  for (i1 in 0:(N.T)) {
    for (i2 in 0:(N.C)) {
      mat[i1 + 1, i2 + 1] = chan_pval(
        x.T = i1,
        x.C = i2,
        N.T = N.T,
        N.C = N.C,
        Delta0 = Delta0,
        method = method
      )
    }
  }
  mat
}

#' Exact-Corrected Test Statistic
#'
#' This function calculates the Exact-Corrected test statistic
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param Delta  numeric between -1 and 1 representing the constraint for P.C-P.T
#' @param Delta0  numeric between 0 and 1 representing the non-inferiority margin
#' @return numeric representing the value of the exact-corrected test statistic
#' @examples
#' #First example calculates the exact-corrected test statistic for
#' #the Rodary et al. (1989) study with proportion of success in the treatment
#' #group being 83/88 and 69/76 for the control with a 10%
#' #noninferiority margin.
#' #Second example calculates the exact-corrected test statistic for
#' #the Fries et al. (1993) study with proportion of success in the treatment
#' #group being 8/15 and 3/15 for the control with a 20% noninferiority margin.
#' orderfunEC(x.T=83,x.C=69,N.T=88,N.C=76,Delta=0.1,Delta0=0.1)
#' orderfunEC(x.T=8,x.C=3,N.T=15,N.C=15,Delta=0.2,Delta0=0.2)
#' @references
#' \insertRef{Hawila:21}{EC}
#' @export
orderfunEC = function(x.T, x.C, N.T, N.C, Delta, Delta0) {
  ECval=0
  rml = restricted_ml(x.T = x.T, x.C = x.C, N.T = N.T, N.C = N.C, Delta = Delta0)
  MR.T = rml$MR.T
  MR.C = rml$MR.C
  DEN_obs <- (MR.T * (1 - MR.T) / N.T + MR.C * (1 - MR.C) / N.C)
  pval=chan_pval(x.T=x.T,x.C=x.C,N.T=N.T,N.C=N.C,Delta0=Delta0,method="MN",lower=TRUE)
  if(pval>1e-10 & pval<(1-1e-10)) {
    ECval=sqrt(DEN_obs)*(orderfun(x.T=x.T, x.C=x.C, N.T=N.T, N.C=N.C,Delta0=Delta0,method="MN")-qnorm(1-pval))
  }

  rml = restricted_ml(x.T = x.T, x.C = x.C, N.T = N.T, N.C = N.C, Delta = Delta)
  MR.T = rml$MR.T
  MR.C = rml$MR.C
  DEN <- (MR.T * (1 - MR.T) / N.T + MR.C * (1 - MR.C) / N.C)
  ECterm=ECval/sqrt(DEN)

  list(Z=(orderfun(x.T=x.T, x.C=x.C, N.T=N.T, N.C=N.C, Delta0=Delta,method="MN"))-ECterm, ECterm=ECterm,ratio=DEN_obs/DEN)
}

#'Null Likelihood
#'
#' This function evaluates the Binomial likelihood under P.T-P.C=-Delta0
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param P.T  numeric between 0 and 1 representing the proportion of responders in the treatment group
#' @param Delta0  numeric between 0 and 1 representing the non-inferiority margin
#' @return numeric representing the probability of getting the observed outcome under the null hypothesis
#' @examples
#' #The probability of getting 10/20 successes in the treatment group
#' #and 8/20 in the placebo group when the proportion of responders
#' #is 30% and the noninferiority margin is 10% can be calculated by
#' likelihood_null(x.T=10, x.C=8, N.T=20, N.C=20, P.T=0.3, Delta0=0.1)
#' @export
likelihood_null = function(x.T, x.C, N.T, N.C, P.T, Delta0) {
  P.C = P.T + Delta0
  dbinom(x.T, N.T, P.T) * dbinom(x.C, N.C, P.C)
}

#' Level of Chan's Exact p-value
#'
#' This function calculates the level of the exact p-value based on Chan (1998)
#'
#' @param alpha  numeric between 0 and 1 representing the significance level
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param Delta0  numeric between 0 and 1 representing the non-inferiority margin
#' @param method character representing the method for ordering criterion("MN","FM","SS","Blackwelder")
#' @return numeric representing the level using Chan's p-value method
#' @examples
#' #Level of Chan's p-value method for the Fries et al. study (1993)
#' chan_level(alpha=0.05,N.T=15,N.C=15,Delta0=0.2,method="MN")
#' @references
#' \insertRef{Chan:03}{EC}
#' @export
chan_level = function(alpha, N.T, N.C, Delta0, method) {
  mat = chan_mat(
    N.T = N.T,
    N.C = N.C,
    Delta0 = Delta0,
    method = method
  )
  ind = which(mat <= alpha, arr.ind = T)
  P.T.vec = seq(0, min(1, 1 - Delta0), by = .1)
  res = rep(NA, length(P.T.vec))
  for (ptvar in 1:length(P.T.vec)) {
    tot = 0
    for (indvar in 1:dim(ind)[1]) {
      tot = tot + likelihood_null(
        x.T = ind[indvar, 1]-1,
        x.C = ind[indvar, 2]-1,
        N.T = N.T,
        N.C = N.C,
        P.T = P.T.vec[ptvar],
        Delta0 = Delta0
      )
    }

    res[ptvar] = tot
  }
  list(mat = mat, dimind = dim(ind), res = res, level = max(res))
}

#' Non-inferiority Confidence Intervals
#'
#' This function computes the confidence interval limits given an ordering criterion with or
#' without the exact-correction
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param Delta0  numeric between 0 and 1 representing the non-inferiority margin
#' @param method character representing the method for ordering criterion("MN","FM","SS","Blackwelder")
#' @param EC logical. TRUE for the exact-corrected confidence limits. FALSE for default method without exact-correction
#' @param alpha  numeric between 0 and 1 representing the significance level
#' @param tol positive numeric representing the tolerance for convergence
#' @return list of length 2 (D.lower, D.upper) representing the lower and upper confidence limits
#' @examples
#' #These two examples demonstrate the confidence intervals for the
#' #Rodary et al. study with and without the exact-correction.
#' confintZ(x.T=83,x.C=69,N.T=88,N.C=76,Delta0=0.1,method="MN", EC=TRUE)
#' confintZ(x.T=83,x.C=69,N.T=88,N.C=76,Delta0=0.1,method="MN", EC=FALSE)
#' @references
#' \insertRef{Hawila:21}{EC}
#'
#' \insertRef{Miettinen:85}{EC}
#'
#' \insertRef{Farrington:90}{EC}
#' @export
confintZ = function(x.T, x.C, N.T, N.C, Delta0, method, EC, alpha=.05, tol=1e-10) {
  ECval=0
  count=0
  if(EC==TRUE){
    rml = restricted_ml(x.T = x.T, x.C = x.C, N.T = N.T, N.C = N.C, Delta = Delta0)
    MR.T = rml$MR.T
    MR.C = rml$MR.C
    DEN_obs <- (MR.T * (1 - MR.T) / N.T + MR.C * (1 - MR.C) / N.C)
    pval=chan_pval(x.T=x.T,x.C=x.C,N.T=N.T,N.C=N.C,Delta0=Delta0,method=method,lower = TRUE)
    if(pval<0.999999999999999) {
      ECval=sqrt(DEN_obs)*(orderfun(x.T=x.T, x.C=x.C, N.T=N.T, N.C=N.C, Delta0=Delta0, method)-qnorm(1-pval))
    }
    if(pval>0.999999999999999) {
      count=100
    }
  }

  chi2=qchisq(1-alpha,1)

  if(count!=100) {
    D1=x.T/N.T-x.C/N.C
    D2=-1
    repeat{
      count=count+1
      Dmid=(D1+D2)/2
      rml = restricted_ml(x.T = x.T, x.C = x.C, N.T = N.T, N.C = N.C, Delta = -Dmid)
      MR.T = rml$MR.T
      MR.C = rml$MR.C
      DEN <- (MR.T * (1 - MR.T) / N.T + MR.C * (1 - MR.C) / N.C)
      chi2.stat=ifelse(DEN>0, (orderfun(x.T=x.T,x.C=x.C,N.T=N.T,N.C=N.C,Delta0=-Dmid,method=method)-ECval/sqrt(DEN))^2, 0)

      if(chi2.stat < chi2){D1=Dmid}
      if(chi2.stat > chi2){D2=Dmid}

      if(abs(chi2.stat-chi2)<tol | count>=70 | chi2.stat>1e8){break}
    }

    D.lower=Dmid

    D1=x.T/N.T-x.C/N.C
    D2=1
    count=0
    repeat{
      count=count+1
      Dmid=(D1+D2)/2
      rml = restricted_ml(x.T = x.T, x.C = x.C, N.T = N.T, N.C = N.C, Delta = -Dmid)
      MR.T = rml$MR.T
      MR.C = rml$MR.C
      DEN <- (MR.T * (1 - MR.T) / N.T + MR.C * (1 - MR.C) / N.C)
      chi2.stat=ifelse(DEN>0, (orderfun(x.T=x.T,x.C=x.C,N.T=N.T,N.C=N.C,Delta0=-Dmid,method=method)-ECval/sqrt(DEN))^2, 0)

      if(chi2.stat < chi2){D1=Dmid}
      if(chi2.stat > chi2){D2=Dmid}

      if(abs(chi2.stat-chi2)<tol | count>=70 | chi2.stat>1e8){break}
    }

    D.upper=Dmid
  }

  if(count==100) {
    D.lower = -1
    D.upper = 1
  }

  list(D.lower=D.lower,D.upper=D.upper,count=count)
}


#' Chan and Zhang Confidence Interval
#'
#' This function calculates the confidence interval based on Chan and Zhang (1999)
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param method character representing the method for ordering criterion("MN","FM","SS","Blackwelder")
#' @param alpha  numeric between 0 and 1 representing the significance level
#' @param tol positive numeric representing the increment size for domain of Delta. Default is set to 0.001.
#' @param width positive numeric representing the range from starting values based on Miettinen & Nurminen confidence limits.
#' @return list of length 2 (D.lower, D.upper) representing the lower and upper Chan and Zhang confidence limits
#' @examples
#'#Chan & Zhang confidence interval for the Rodary et al. study. (1989)
#'#chan_zhang(x.T=83,x.C=69,N.T=88,N.C=76, method="MN")
#'@references
#' \insertRef{Chan:99}{EC}
#' @export
chan_zhang <- function(x.T, x.C, N.T, N.C, method, alpha=.05, tol=1e-3, width=0.3) {
  d_LL = confintZ(x.T, x.C, N.T, N.C, Delta0=0, method,EC=F)$D.lower
  d_UL = confintZ(x.T, x.C, N.T, N.C, Delta0=0, method,EC=F)$D.upper

  deltaL=seq(max(-d_LL-width,-0.9999),min(-d_LL+width,0.9999),tol)
  deltaU=seq(max(-d_UL-width,-0.9999),min(-d_UL+width,0.9999),tol)
  probsL=rep(NA, length(deltaL))
  probsU=rep(NA, length(deltaL))

    for (i in 1:length(deltaL)) probsL[i]=chan_pval(x.T, x.C, N.T, N.C, deltaL[i], method ="MN", lower=T)
    for (i in 1:length(deltaU)) probsU[i]=chan_pval(x.T, x.C, N.T, N.C, deltaU[i],method = "MN", lower=F)

    if(length(deltaL[probsL>alpha/2])>0) LB=max(deltaL[probsL>alpha/2],na.rm=T)
    if(length(deltaL[probsL>alpha/2])==0) LB=1
    if(length(deltaU[probsU>alpha/2])>0) UB=min(deltaU[probsU>alpha/2],na.rm=T)
    if(length(deltaU[probsU>alpha/2])==0) UB=-1

  D.lower=-LB
  D.upper=-UB

  list(D.lower=D.lower, D.upper=D.upper)
}

#' Level of Confidence Interval
#'
#' This function evaluates true level of decision rule based on confidence interval method
#'
#' @param alpha  numeric between 0 and 1 representing the significance level
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param Delta0  numeric between 0 and 1 representing the non-inferiority margin
#' @param method character representing the method for ordering criterion("MN","FM","SS","Blackwelder")
#' @param EC logical. TRUE for the exact-corrected confidence limits. FALSE for default method without exact-correction. Only relevant if CZ=F
#' @param tolEC positive numeric representing the tolerance for confidence interval convergence
#' @param CZ logical. TRUE for Chan and Zhang confidence limits. FALSE for Exact-corrected (EC=T) or default method without correction (EC=F)
#' @param tolCZ positive numeric representing the increment size for domain of Delta. Default is set to 0.001.
#' @param width positive numeric representing the range from starting values based on Miettinen & Nurminen confidence limits.
#' @return numeric representing the level of test based on the specified confidence interval method
#' @examples
#' #The three examples calculate the level of the exact-corrected,
#' #Chan & Zhang and Asymptotic confidence interval based on the
#' #Miettenen and Nurminen test statistic where each group has a
#' #sample size of 10, alpha is 0.1 and the noninferiority margin
#' #is 20%.
#' ci_level(alpha=0.1,N.T=10,N.C=10,Delta0=0.2,method="MN",EC=TRUE,tolEC=1e-4,CZ=FALSE,tolCZ=1e-3,width=1e-3)
#' #ci_level(alpha=0.1,N.T=10,N.C=10,Delta0=0.2,method="MN",EC=FALSE,tolEC=1e-4,CZ=TRUE,tolCZ=1e-3,width=1e-3)
#' ci_level(alpha=0.1,N.T=10,N.C=10,Delta0=0.2,method="MN",EC=FALSE,tolEC=1e-4,CZ=FALSE,tolCZ=1e-3,width=1e-3)
#' @references
#' \insertRef{Hawila:21}{EC}
#' \insertRef{Chan:99}{EC}
#' \insertRef{Miettinen:85}{EC}

#' @export
ci_level <- function(alpha, N.T, N.C, Delta0, method, EC, tolEC, CZ, tolCZ, width) {
  M=matrix(0,N.T+1,N.C+1)
  for(i in 0:N.T) {
    for(j in 0:N.C) {
      if(CZ) {
        ci = chan_zhang(x.T=i,x.C=j,N.T=N.T,N.C=N.C,method=method,tol=tolCZ,width=width,alpha=alpha)
        lb=ci$D.lower
        ub=ci$D.upper
      }
      if(!CZ) {
        ci = confintZ(x.T=i,x.C=j,N.T=N.T,N.C=N.C,Delta0=Delta0,method=method,EC=EC,alpha=alpha,tol=tolEC)
        lb=ci$D.lower
        ub=ci$D.upper
        count=ci$count
      }
      ci.rej = (lb>=-Delta0)

      if(!CZ) {
        if(lb>(-Delta0) & count<70) {
          M[i+1,j+1]=1
        }}
      if(CZ) {
        if(lb>(-Delta0)) {
          M[i+1,j+1]=1
        }}

    }

  }
  ind = which(M==1, arr.ind = T)
  P.T.vec = seq(0, min(1, 1 - Delta0), by = .1)
  res = rep(NA, length(P.T.vec))

  if(dim(ind)[1]==0) {
    level=0
  }
  if(dim(ind)[1]>0) {
    for (ptvar in 1:length(P.T.vec)) {
      tot = 0
      for (indvar in 1:dim(ind)[1]) {
        tot = tot + likelihood_null(
          x.T = ind[indvar, 1]-1,
          x.C = ind[indvar, 2]-1,
          N.T = N.T,
          N.C = N.C,
          P.T = P.T.vec[ptvar],
          Delta0 = Delta0
        )
      }
      res[ptvar] = tot
      level=max(res)
    }
  }
  list(level=level)
}

#' @importFrom Rdpack reprompt
