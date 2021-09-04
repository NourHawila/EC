#' Restricted MLE of variance
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
#' restricted_mle(x.T=5,x.C=7,N.T=10,N.C=10,Delta=0.2)
#' @references
#' \insertRef{Miettinen:85}{EC}
#' @export
restricted_mle = function(x.T, x.C, N.T, N.C, Delta) {
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



#' Test statistics for noninferiority hypothesis
#'
#' This function  calculate statistics for testing noninferiority based on Santner & Snell (1980), Blackwelder (1982),
#' Miettinen & Nurminen (1985) and Farrington & Manning (1990).
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param delta0  numeric between 0 and 1 representing the noninferiority margin
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
#' stat_general(x.T=83,x.C=69,N.T=88,N.C=76,delta0=0.1,method="MN")
#' stat_general(x.T=8,x.C=3,N.T=15,N.C=15,delta0=0.2,method="FM")
#' @references
#' \insertRef{Santner:80}{EC}
#'
#' \insertRef{Blackwelder:82}{EC}
#'
#' \insertRef{Miettinen:85}{EC}
#'
#' \insertRef{Farrington:90}{EC}
#' @export
stat_general = function(x.T, x.C, N.T, N.C, delta0, method) {
  rml = restricted_mle(
    x.T = x.T,
    x.C = x.C,
    N.T = N.T,
    N.C = N.C,
    Delta = delta0
  )

  N = N.T+N.C
  R.T = x.T/N.T
  R.C = x.C/N.C

  if (method == "MN") {
    MR.T = rml$MR.T
    MR.C = rml$MR.C
    DEN = (MR.T*(1-MR.T)/N.T + MR.C*(1-MR.C)/N.C)

    mn.res = ifelse(DEN<=0, 0, (R.T-R.C+delta0) / sqrt(DEN))
    out = mn.res
  }

  if (method == "FM") {
    MR.T = rml$MR.T
    MR.C = rml$MR.C
    DEN <- (MR.T*(1-MR.T)/N.T + MR.C*(1-MR.C)/N.C)

    fm.res = ifelse(DEN<=0, 0, sqrt(N/(N-1))*((R.T-R.C+delta0) / sqrt(DEN)))
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



#' Order statistic for all 2x2 tables
#'
#' This function calculate the order statistics for all 2x2 tables
#'
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param delta0  numeric between 0 and 1 representing the noninferiority margin
#' @param method character representing the method for ordering criterion("MN","FM","SS","Blackwelder")
#' @return array of dimensions (N.T+1)x(N.C+1) where the (i,j) element is the order statistic for x.T=i and x.C=j
#' @examples
#' #16x16 array of the the Santner & Snell test statistic for all
#' #possible 2x2 tables arising from a sample size of 15 in the
#' #treatment and control group and 30% noninferiority margin.
#' order_mat(N.T=15,N.C=15,delta0=0.3,method="SS")
#' @references
#' \insertRef{Santner:80}{EC}
#' @export
order_mat = function(N.T, N.C, delta0, method) {
  res = array(NA, dim = c(N.T + 1, N.C + 1))

  for (i1 in 0:(N.T)) {
    for (i2 in 0:(N.C)) {
      res[i1 + 1, i2 + 1] = stat_general(
        x.T = i1,
        x.C = i2,
        N.T = N.T,
        N.C = N.C,
        delta0 = delta0,
        method = method
      )
    }
  }
  res
}

#' Verifying Barnard's criterion
#'
#' This function checks whether Barnard's criterion (Biometrika, 1947) is satisfied
#'
#' @param mat array of dimensions (N.T+1)x(N.C+1) where the (i,j) element is the order statistic for x.T=i and x.C=j
#' @return logical. TRUE if Barnard's criterion is satisfied, FALSE otherwise
#' @examples
#' #Checks if the Barnard criterion is satisfied for the
#' #Miettenin & Nurminen ordering statistic with N.T=15,
#' #N.C=15 and delta0=30%
#' Barnard_check(order_mat(N.T=15,N.C=15,delta0=0.3,method="MN"))
#' @export
#' @references
#' \insertRef{Barnard:47}{EC}
Barnard_check=function(mat){
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

#' p-value corresponding to the Miettinen & Nurminen (1985) confidence interval
#'
#' This function calculates the asymptotic p-value described in Hawila & Berg (2021)
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param delta0  numeric between 0 and 1 representing the noninferiority margin
#' @return numeric representing the asymptotic p-value associated with the Miettinen & Nurminen statistic
#' @examples
#' #The example is taken from Rodary et al. (1989) which was analyzed by Chan (1998)
#' pval_MN(x.T=83,x.C=69,N.T=88,N.C=76,delta0=0.1)
#' @references
#' \insertRef{Miettinen:85}{EC}
#' \insertRef{Hawila:21}{EC}
#' @export
pval_MN <- function(x.T, x.C, N.T, N.C, delta0) {
  1-pnorm(stat_general(x.T, x.C, N.T, N.C, delta0,method="MN"))
}

#' p-value corresponding to the Wald confidence interval
#'
#' This function calculates the Wald p-value
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param delta0  numeric between 0 and 1 representing the noninferiority margin
#' @return numeric representing the asymptotic p-value associated with the Miettinen & Nurminen statistic
#' @examples
#' #The example is taken from Rodary et al. (1989) which was analyzed by Chan (1998)
#' pval_MN(x.T=83,x.C=69,N.T=88,N.C=76,delta0=0.1)
#' @export
pval_Wald <- function(x.T, x.C, N.T, N.C, delta0) {
  p.T = x.T/N.T
  p.C = x.C/N.C
  se = sqrt(p.T*(1-p.T)/N.T + p.C*(1-p.C)/N.C)
  z = (p.T-p.C+delta0)/se
  1-pnorm(z)
}

#' Exact p-value of Chan (1998)
#'
#' This function calculates the exact p-value based on Chan (1998)
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param delta0  numeric between 0 and 1 representing the noninferiority margin
#' @param method character representing the method for ordering criterion("MN","FM","SS","Blackwelder")
#' @param lower logical. TRUE for the null hypothesis P.T-P.C<=-delta0.FALSE for the null hypothesis P.T-P.C>-delta0
#' @param tol positive numeric representing the increment size for domain of Delta. Default is set to 0.001.
#' @return numeric representing Chan's exact p-value
#' @examples
#' #The first example is taken from Rodary et al. (1989) which was analyzed by Chan (1998)
#' #The second example is taken from Hawila & Berg (2021)
#' pval_Chan(x.T=83,x.C=69,N.T=88,N.C=76,delta0=0.1,method="MN")
#' pval_Chan(x.T=2,x.C=0,N.T=15,N.C=10,delta0=0.12,method="MN")
#' @references
#' \insertRef{Chan:98}{EC}
#' @export
pval_Chan <- function(x.T, x.C, N.T, N.C, delta0, method="MN", lower = TRUE,tol=1e-3) {
  myorder_mat = order_mat(N.T, N.C, delta0 = delta0, method = method)
  obs = myorder_mat[x.T + 1, x.C + 1]

  P.T.vec = seq(max(0,-delta0), min(1, 1 - delta0), by = tol)
  pvals = rep(NA, length(P.T.vec))
  for (ptvar in 1:length(P.T.vec)) {
    P.T = P.T.vec[ptvar]
    P.C = P.T + delta0

    if(lower==TRUE){
      pvals[ptvar] = sum((myorder_mat>=obs) * dbinom(0:N.T, N.T, P.T) %o%dbinom(0:N.C, N.C, P.C))
    }
    if(lower==FALSE){
      pvals[ptvar] = sum((myorder_mat<=obs) * dbinom(0:N.T, N.T, P.T) %o%dbinom(0:N.C, N.C, P.C))

    }
  }
  max(pvals)

}


#' p-value corresponding to the Chan & Zhang confidence interval as described in Hawila & Berg (2021)
#'
#' This function computes the Chan & Zhang p-value given in Hawila & Berg (2021)
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param delta0  numeric between 0 and 1 representing the noninferiority margin
#' @param method character representing the method for ordering criterion("MN","FM","SS","Blackwelder")
#' @param tol positive numeric representing the tolerance for calculation
#' @return numeric representing the Chan & Zhang p-value
#' @examples
#' #This example is taken from Hawila & Berg (2021)
#' pval_CZ(x.T=2,x.C=0,N.T=15,N.C=10,delta0=0.12,method="MN")
#' @references
#' \insertRef{Hawila:21}{EC}
#' @export
pval_CZ <- function(x.T, x.C, N.T, N.C, delta0, method="MN",tol=1e-3){
  deltas=seq(-0.999,-delta0,length=1/tol)
  res=array(dim=length(deltas))
  for(d in 1:length(deltas)){
    res[d] = pval_Chan(x.T,x.C,N.T,N.C,-deltas[d],method)
  }
  max(res)
}

#' Chan (1998) exact p-value for all 2x2 tables
#'
#' This function calculates the Chan p-values for all 2x2 tables given N.T, N.C and an ordering criterion;
#' This function is primarily used in the \code{\link[EC]{size_Chan}} function.
#'
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param delta0  numeric between 0 and 1 representing the noninferiority margin
#' @param method character representing the method for ordering criterion("MN","FM","SS","Blackwelder")
#' @return array of dimension (N.T+1)x(N.C+1)  where the (i,j) element is the exact p-value for x.T=i and x.C=j
#' @examples
#'#8x6 array of the the Chan p-values based on Farrington & Manning
#'#statistic for all possible 2x2 tables arising from a sample size
#'#of 7 in the treatment group and 5 in the control group and 10%
#'#noninferiority margin.
#'mat_Chan(N.T=7,N.C=5,delta0=0.1,method="FM")
#'@references
#' \insertRef{Chan:98}{EC}
#' @seealso [EC::size_Chan]
#' @export
mat_Chan = function(N.T, N.C, delta0, method="MN") {
  mat = array(NA, dim = c(N.T + 1, N.C + 1))
  for (i1 in 0:(N.T)) {
    for (i2 in 0:(N.C)) {
      mat[i1 + 1, i2 + 1] = pval_Chan(
        x.T = i1,
        x.C = i2,
        N.T = N.T,
        N.C = N.C,
        delta0 = delta0,
        method = method
      )
    }
  }
  mat
}

#' Exact-corrected test statistic
#'
#' This function calculates the exact-corrected test statistic
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param Delta  numeric between -1 and 1 representing the constraint for P.C-P.T
#' @param delta0  numeric between 0 and 1 representing the noninferiority margin
#' @return numeric representing the value of the exact-corrected test statistic
#' @examples
#' #First example calculates the exact-corrected test statistic for
#' #the Rodary et al. (1989) study with proportion of success in the treatment
#' #group being 83/88 and 69/76 for the control with a 10%
#' #noninferiority margin.
#' #Second example calculates the exact-corrected test statistic for
#' #the Fries et al. (1993) study with proportion of success in the treatment
#' #group being 8/15 and 3/15 for the control with a 20% noninferiority margin.
#' stat_EC(x.T=83,x.C=69,N.T=88,N.C=76,Delta=0.1,delta0=0.1)
#' stat_EC(x.T=8,x.C=3,N.T=15,N.C=15,Delta=0.2,delta0=0.2)
#' @references
#' \insertRef{Hawila:21}{EC}
#' @export
stat_EC = function(x.T, x.C, N.T, N.C, Delta, delta0) {
  ECval=0
  rml = restricted_mle(x.T = x.T, x.C = x.C, N.T = N.T, N.C = N.C, Delta = delta0)
  MR.T = rml$MR.T
  MR.C = rml$MR.C
  DEN_obs <- (MR.T * (1 - MR.T) / N.T + MR.C * (1 - MR.C) / N.C)
  pval=pval_Chan(x.T=x.T,x.C=x.C,N.T=N.T,N.C=N.C,delta0=delta0,method="MN",lower=TRUE)
  if(pval>1e-10 & pval<(1-1e-10)) {
    ECval=sqrt(DEN_obs)*(stat_general(x.T=x.T, x.C=x.C, N.T=N.T, N.C=N.C,delta0=delta0,method="MN")-qnorm(1-pval))
  }

  rml = restricted_mle(x.T = x.T, x.C = x.C, N.T = N.T, N.C = N.C, Delta = Delta)
  MR.T = rml$MR.T
  MR.C = rml$MR.C
  DEN <- (MR.T * (1 - MR.T) / N.T + MR.C * (1 - MR.C) / N.C)
  ECterm=ECval/sqrt(DEN)

  list(Z=(stat_general(x.T=x.T, x.C=x.C, N.T=N.T, N.C=N.C, delta0=Delta,method="MN"))-ECterm, ECterm=ECterm,ratio=DEN_obs/DEN)
}

#' Binomial likelihood under the null hypothesis
#'
#' This function evaluates the binomial likelihood under P.T-P.C=-delta0
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param P.T  numeric between 0 and 1 representing the proportion of responders in the treatment group
#' @param delta0  numeric between 0 and 1 representing the noninferiority margin
#' @return numeric representing the probability of getting the observed outcome under the null hypothesis
#' @examples
#' #The probability of getting 10/20 successes in the treatment group
#' #and 8/20 in the placebo group when the proportion of responders
#' #is 30% and the noninferiority margin is 10% can be calculated by
#' likelihood_null(x.T=10, x.C=8, N.T=20, N.C=20, P.T=0.3, delta0=0.1)
#' @export
likelihood_null = function(x.T, x.C, N.T, N.C, P.T, delta0) {
  P.C = P.T + delta0
  dbinom(x.T, N.T, P.T) * dbinom(x.C, N.C, P.C)
}

#' Maximal size of Chan's exact p-value
#'
#' This function calculates the maximal size of Chan's (1998) exact p-value; see Hawila & Berg (2021)
#'
#' @param alpha  numeric between 0 and 1 representing the significance level
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param delta0  numeric between 0 and 1 representing the noninferiority margin
#' @param method character representing the method for ordering criterion("MN","FM","SS","Blackwelder")
#' @return numeric representing the level using Chan's p-value method
#' @examples
#' #Level of Chan's p-value method for the Fries et al. study (1993)
#' size_Chan(alpha=0.05,N.T=15,N.C=15,delta0=0.2,method="MN")
#' @references
#' \insertRef{Hawila:21}{EC}
#' \insertRef{Chan:03}{EC}
#' @export
size_Chan = function(alpha, N.T, N.C, delta0, method="MN") {
  mat = mat_Chan(
    N.T = N.T,
    N.C = N.C,
    delta0 = delta0,
    method = method
  )
  ind = which(mat <= alpha, arr.ind = T)
  P.T.vec = seq(0, min(1, 1 - delta0), by = .1)
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
        delta0 = delta0
      )
    }

    res[ptvar] = tot
  }
  list(mat = mat, dimind = dim(ind), res = res, level = max(res))
}

#' Noninferiority confidence intervals
#'
#' This function computes the confidence interval limits given an ordering criterion with or
#' without the exact-correction
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param delta0  numeric between 0 and 1 representing the noninferiority margin
#' @param method character representing the method for ordering criterion("MN","FM","SS","Blackwelder")
#' @param EC logical. TRUE for the exact-corrected confidence limits. FALSE for default method without exact-correction
#' @param alpha  numeric between 0 and 1 representing the significance level
#' @param tol positive numeric representing the tolerance for convergence
#' @return list of length 2 (ci.lower, ci.upper) representing the lower and upper confidence limits
#' @examples
#' #These two examples demonstrate the confidence intervals for the
#' #Rodary et al. study with and without the exact-correction.
#' ci_general(x.T=83,x.C=69,N.T=88,N.C=76,delta0=0.1,method="MN", EC=TRUE)
#' ci_general(x.T=83,x.C=69,N.T=88,N.C=76,delta0=0.1,method="MN", EC=FALSE)
#' @references
#' \insertRef{Hawila:21}{EC}
#'
#' \insertRef{Miettinen:85}{EC}
#'
#' \insertRef{Farrington:90}{EC}
#' @export
ci_general = function(x.T, x.C, N.T, N.C, delta0, method="MN", EC, alpha=.05, tol=1e-10) {
  ECval=0
  count=0
  if(EC==TRUE){
    rml = restricted_mle(x.T = x.T, x.C = x.C, N.T = N.T, N.C = N.C, Delta = delta0)
    MR.T = rml$MR.T
    MR.C = rml$MR.C
    DEN_obs <- (MR.T * (1 - MR.T) / N.T + MR.C * (1 - MR.C) / N.C)
    pval=pval_Chan(x.T=x.T,x.C=x.C,N.T=N.T,N.C=N.C,delta0=delta0,method=method,lower = TRUE)
    if(pval<0.999999999999999) {
      ECval=sqrt(DEN_obs)*(stat_general(x.T=x.T, x.C=x.C, N.T=N.T, N.C=N.C, delta0=delta0, method)-qnorm(1-pval))
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
      rml = restricted_mle(x.T = x.T, x.C = x.C, N.T = N.T, N.C = N.C, Delta = -Dmid)
      MR.T = rml$MR.T
      MR.C = rml$MR.C
      DEN <- (MR.T * (1 - MR.T) / N.T + MR.C * (1 - MR.C) / N.C)
      chi2.stat=ifelse(DEN>0, (stat_general(x.T=x.T,x.C=x.C,N.T=N.T,N.C=N.C,delta0=-Dmid,method=method)-ECval/sqrt(DEN))^2, 0)

      if(chi2.stat < chi2){D1=Dmid}
      if(chi2.stat > chi2){D2=Dmid}

      if(abs(chi2.stat-chi2)<tol | count>=70 | chi2.stat>1e8){break}
    }

    ci.lower=Dmid

    D1=x.T/N.T-x.C/N.C
    D2=1
    count=0
    repeat{
      count=count+1
      Dmid=(D1+D2)/2
      rml = restricted_mle(x.T = x.T, x.C = x.C, N.T = N.T, N.C = N.C, Delta = -Dmid)
      MR.T = rml$MR.T
      MR.C = rml$MR.C
      DEN <- (MR.T * (1 - MR.T) / N.T + MR.C * (1 - MR.C) / N.C)
      chi2.stat=ifelse(DEN>0, (stat_general(x.T=x.T,x.C=x.C,N.T=N.T,N.C=N.C,delta0=-Dmid,method=method)-ECval/sqrt(DEN))^2, 0)

      if(chi2.stat < chi2){D1=Dmid}
      if(chi2.stat > chi2){D2=Dmid}

      if(abs(chi2.stat-chi2)<tol | count>=70 | chi2.stat>1e8){break}
    }

    ci.upper=Dmid
  }

  if(count==100) {
    ci.lower = -1
    ci.upper = 1
  }

  list(ci.lower=ci.lower,ci.upper=ci.upper,count=count)
}


#' Confidence interval using Chan & Zhang (1999) methodology
#'
#' This function calculates the confidence interval based on Chan & Zhang (1999) given an ordering
#' criterion specified in the method argument. The default ordering criterion (method="MN") is
#' the ordering criterion used in Chan & Zhang (1999).
#'
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param method character representing the method for ordering criterion("MN","FM","SS","Blackwelder")
#' @param alpha  numeric between 0 and 1 representing the significance level
#' @param tol positive numeric representing the increment size for domain of Delta. Default is set to 0.001.
#' @param eps positive numeric representing the range from starting values based on Miettinen & Nurminen confidence limits.
#' @return list of length 2 (ci.lower, ci.upper) representing the lower and upper Chan and Zhang confidence limits
#' @examples
#'#Chan & Zhang confidence interval for the Rodary et al. study. (1989)
#'#ci_CZ(x.T=83,x.C=69,N.T=88,N.C=76, method="MN",alpha=0.05)
#'@references
#' \insertRef{Chan:99}{EC}
#' @export
ci_CZ <- function(x.T, x.C, N.T, N.C, method="MN", alpha=.05, tol=1e-3, eps=0.3) {
  d_LL = ci_general(x.T, x.C, N.T, N.C, delta0=0, method,EC=F)$ci.lower
  d_UL = ci_general(x.T, x.C, N.T, N.C, delta0=0, method,EC=F)$ci.upper

  deltaL=seq(max(-d_LL-eps,-0.9999),min(-d_LL+eps,0.9999),tol)
  deltaU=seq(max(-d_UL-eps,-0.9999),min(-d_UL+eps,0.9999),tol)
  probsL=rep(NA, length(deltaL))
  probsU=rep(NA, length(deltaL))

    for (i in 1:length(deltaL)) probsL[i]=pval_Chan(x.T, x.C, N.T, N.C, deltaL[i], method ="MN", lower=T)
    for (i in 1:length(deltaU)) probsU[i]=pval_Chan(x.T, x.C, N.T, N.C, deltaU[i],method = "MN", lower=F)

    if(length(deltaL[probsL>alpha/2])>0) LB=max(deltaL[probsL>alpha/2],na.rm=T)
    if(length(deltaL[probsL>alpha/2])==0) LB=1
    if(length(deltaU[probsU>alpha/2])>0) UB=min(deltaU[probsU>alpha/2],na.rm=T)
    if(length(deltaU[probsU>alpha/2])==0) UB=-1

  ci.lower=-LB
  ci.upper=-UB

  list(ci.lower=ci.lower, ci.upper=ci.upper)
}

#' Maximal size of tests based on various confidence interval methods
#'
#' This function evaluates maximum size of the decision rule based on a given confidence interval method; see Hawila & Berg (2021)
#'
#' @param alpha  numeric between 0 and 1 representing the significance level
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param delta0  numeric between 0 and 1 representing the noninferiority margin
#' @param method character representing the method for ordering criterion("MN","FM","SS","Blackwelder")
#' @param EC logical. TRUE for the exact-corrected confidence limits. FALSE for default method without exact-correction. Only relevant if CZ=F
#' @param tolEC positive numeric representing the tolerance for confidence interval convergence
#' @param CZ logical. TRUE for Chan and Zhang confidence limits. FALSE for Exact-corrected (EC=T) or default method without correction (EC=F)
#' @param tolCZ positive numeric representing the increment size for domain of Delta. Default is set to 0.001.
#' @param eps positive numeric representing the range from starting values based on Miettinen & Nurminen confidence limits.
#' @return numeric representing the level of test based on the specified confidence interval method
#' @examples
#' #The three examples calculate the level of the exact-corrected,
#' #Chan & Zhang and Asymptotic confidence interval based on the
#' #Miettenen and Nurminen test statistic where each group has a
#' #sample size of 10, alpha is 0.1 and the noninferiority margin
#' #is 20%.
#' size_general(alpha=0.1,N.T=10,N.C=10,delta0=0.2,method="MN",EC=TRUE,tolEC=1e-4,CZ=FALSE,tolCZ=1e-3,eps=1e-3)
#' #size_general(alpha=0.1,N.T=10,N.C=10,delta0=0.2,method="MN",EC=FALSE,tolEC=1e-4,CZ=TRUE,tolCZ=1e-3,eps=1e-3)
#' size_general(alpha=0.1,N.T=10,N.C=10,delta0=0.2,method="MN",EC=FALSE,tolEC=1e-4,CZ=FALSE,tolCZ=1e-3,eps=1e-3)
#' @references
#' \insertRef{Hawila:21}{EC}
#' \insertRef{Chan:99}{EC}
#' \insertRef{Miettinen:85}{EC}
#' @export
size_general <- function(alpha, N.T, N.C, delta0, method="MN", EC, tolEC=1e-6, CZ, tolCZ=1e-3, eps=0.3) {
  M=matrix(0,N.T+1,N.C+1)
  count=0
  for(i in 0:N.T) {
    for(j in 0:N.C) {
      if(CZ) {
        ci = ci_CZ(x.T=i,x.C=j,N.T=N.T,N.C=N.C,method=method,tol=tolCZ,eps=eps,alpha=alpha)
        lb=ci$ci.lower
        ub=ci$ci.upper
      }
      if(!CZ) {
        if(method !="Wald") ci = ci_general(x.T=i,x.C=j,N.T=N.T,N.C=N.C,delta0=delta0,method=method,EC=EC,alpha=alpha,tol=tolEC)
        if(method == "Wald") ci = ci_Wald(x.T=i,x.C=j,N.T=N.T,N.C=N.C,alpha=alpha)
        lb=ci$ci.lower
        ub=ci$ci.upper
        if(method != "Wald") count=ci$count
      }
      ci.rej = (lb>=-delta0)

      if(!CZ) {
        if(lb>(-delta0) & count<70) {
          M[i+1,j+1]=1
        }}
      if(CZ) {
        if(lb>(-delta0)) {
          M[i+1,j+1]=1
        }}

    }

  }
  ind = which(M==1, arr.ind = T)
  P.T.vec = seq(0, min(1, 1 - delta0), by = .1)
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
          delta0 = delta0
        )
      }
      res[ptvar] = tot
      level=max(res)
    }
  }
  list(level=level)
}


#' Confidence interval using the exact-corrected methodology (Hawila & Berg, 2021)
#'
#' This is a wrapper function of \code{\link[EC]{ci_general}}
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param delta0  numeric between 0 and 1 representing the noninferiority margin
#' @param alpha  numeric between 0 and 1 representing the significance level
#' @param tol positive numeric representing the tolerance for convergence
#' @return list of length 2 (ci.lower, ci.upper) representing the lower and upper confidence limits
#' @examples
#' #This examples demonstrates the confidence intervals for the
#' #Rodary et al. study using the exact-corrected method
#' ci_EC(x.T=83,x.C=69,N.T=88,N.C=76,delta0=0.1,alpha=0.05)
#' @references
#' \insertRef{Hawila:21}{EC}
#' @seealso [EC::ci_general]
#' @export
ci_EC <- function(x.T, x.C, N.T, N.C, delta0,alpha = 0.05, tol = 1e-10){
  ci_general(x.T=x.T, x.C=x.C, N.T=N.T, N.C=N.C, delta0=delta0, method="MN", EC=T, alpha=alpha,tol=tol)
}

#' Confidence interval using the Miettinen & Nurminen (1985) methodology
#'
#' This is a wrapper function of \code{\link[EC]{ci_general}}
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param alpha  numeric between 0 and 1 representing the significance level
#' @param tol positive numeric representing the tolerance for convergence
#' @return list of length 2 (ci.lower, ci.upper) representing the lower and upper confidence limits
#' @examples
#' #This examples demonstrates the confidence intervals for the
#' #Rodary et al. study using the exact-corrected method
#' ci_MN(x.T=83,x.C=69,N.T=88,N.C=76,alpha=0.05)
#' @references
#' \insertRef{Miettinen:85}{EC}
#' @seealso [EC::ci_general]
#' @export
ci_MN <- function(x.T, x.C, N.T, N.C,alpha = 0.05, tol = 1e-10){
  ci_general(x.T=x.T, x.C=x.C, N.T=N.T, N.C=N.C, delta0=0, method="MN", EC=F, alpha=alpha,tol=tol)
}

#' Confidence interval using the Wald methodology
#'
#'
#' @param x.T  positive integer representing the observed number of responders in the treatment group
#' @param x.C  positive integer representing the observed number of responders in the control group
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param alpha  numeric between 0 and 1 representing the significance level
#' @return list of length 2 (ci.lower, ci.upper) representing the lower and upper confidence limits
#' @examples
#' #This examples demonstrates the confidence intervals for the
#' #Kim et al. study using Wald's method
#' ci_Wald(x.T=173,x.C=174,N.T=181,N.C=181,alpha=0.05)
#' @references
#' \insertRef{Altman:13}{EC}
#' \insertRef{Fagerland:15}{EC}
#' @export
ci_Wald <- function(x.T, x.C, N.T, N.C,alpha = 0.05){
  p.T = x.T/N.T
  p.C = x.C/N.C
  se = sqrt(p.T*(1-p.T)/N.T + p.C*(1-p.C)/N.C)
  ci.lower = p.T-p.C - qnorm(1-alpha/2)*se
  ci.upper = p.T-p.C + qnorm(1-alpha/2)*se
  return(list(ci.lower=ci.lower, ci.upper=ci.upper))
}



#' Maximal size of the exact-corrected confidence interval
#'
#' This function evaluates maximum size of the decision rule based on the exact-corrected confidence interval
#'
#' @param alpha  numeric between 0 and 1 representing the significance level
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param delta0  numeric between 0 and 1 representing the noninferiority margin
#' @param tolEC positive numeric representing the tolerance for confidence interval convergence
#' @return numeric representing the level of test based on the specified confidence interval method
#' @examples
#' #The three examples calculate the level of the exact-corrected,
#' #Chan & Zhang and Asymptotic confidence interval based on the
#' #Miettenen and Nurminen test statistic where each group has a
#' #sample size of 10, alpha is 0.1 and the noninferiority margin
#' #is 20%.
#' size_EC(alpha=0.1,N.T=10,N.C=10,delta0=0.2,tolEC=1e-4)
#' @references
#' \insertRef{Hawila:21}{EC}
#' @export
size_EC <- function(alpha, N.T, N.C, delta0,tolEC=1e-4){
size_general(alpha, N.T, N.C, delta0, method="MN", EC=T, tolEC=tolEC,CZ=F)
}



#' Maximal size of the Chan & Zhang (1999) confidence interval
#'
#' This function evaluates maximum size of the decision rule based on Chan & Zhang's confidence interval
#'
#' @param alpha  numeric between 0 and 1 representing the significance level
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param delta0  numeric between 0 and 1 representing the noninferiority margin
#' @param tolCZ positive numeric representing the increment size for domain of Delta. Default is set to 0.001.
#' @param eps positive numeric representing the range from starting values based on Miettinen & Nurminen confidence limits.
#' @return numeric representing the level of test based on the specified confidence interval method
#' @examples
#' #The three examples calculate the level of the exact-corrected,
#' #Chan & Zhang and Asymptotic confidence interval based on the
#' #Miettenen and Nurminen test statistic where each group has a
#' #sample size of 10, alpha is 0.1 and the noninferiority margin
#' #is 20%.
#' size_CZ(alpha=0.1,N.T=10,N.C=10,delta0=0.2,tolCZ=1e-3)
#' @references
#' \insertRef{Hawila:21}{EC}
#' @export
size_CZ <- function(alpha, N.T, N.C, delta0,tolCZ=1e-3, eps=1e-3){
  size_general(alpha, N.T, N.C, delta0, method="MN", EC=F, tolCZ=tolCZ,CZ=T,eps=eps)
}

#' Maximal size of the Miettinen & Nurminen (1985) confidence interval
#'
#' This function evaluates maximum size of the decision rule based on the Miettinen & Nurminen confidence interval
#'
#' @param alpha  numeric between 0 and 1 representing the significance level
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param delta0  numeric between 0 and 1 representing the noninferiority margin
#' @param tolEC positive numeric representing the tolerance for confidence interval convergence
#' @return numeric representing the level of test based on the specified confidence interval method
#' @examples
#' #The example calculates the size of the confidence interval based on the
#' #Miettenen and Nurminen test statistic where each group has a
#' #sample size of 10, alpha is 0.1 and the noninferiority margin
#' #is 20%.
#' size_MN(alpha=0.1,N.T=10,N.C=10,delta0=0.2,tolEC=1e-4)
#' @references
#' \insertRef{Hawila:21}{EC}
#' @export
size_MN <- function(alpha, N.T, N.C, delta0,tolEC=1e-4){
  size_general(alpha, N.T, N.C, delta0, method="MN", EC=F, tolEC=tolEC,CZ=F)
}


#' Maximal size of the Wald confidence interval
#'
#' This function evaluates maximum size of the decision rule based on the Wald confidence interval
#'
#' @param alpha  numeric between 0 and 1 representing the significance level
#' @param N.T  positive integer representing the sample size in the treatment group
#' @param N.C  positive integer representing the sample size in the control group
#' @param delta0  numeric between 0 and 1 representing the noninferiority margin
#' @return numeric representing the level of test based on the specified confidence interval method
#' @examples
#' #The example calculates the size of the confidence interval based on the
#' #Wald's method where each group has a
#' #sample size of 10, alpha is 0.1 and the noninferiority margin
#' #is 20%.
#' size_Wald(alpha=0.1,N.T=10,N.C=10,delta0=0.1)
#' @references
#' \insertRef{Hawila:21}{EC}
#' @export
size_Wald <- function(alpha, N.T, N.C,delta0){
  size_general(alpha, N.T, N.C, method="Wald",EC=F,CZ=F,delta0=delta0)
}



#' @importFrom Rdpack reprompt
