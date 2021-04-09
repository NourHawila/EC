#' Restricted MLE of Variance
#'
#' This function calculate the maximum likelihood estimator (MLE) of the variance
#' restricted to P.T-P.C=-Delta
#'
#' @param x.T  observed number of responders in Treatment group
#' @param x.C  observed number of responders in Control group
#' @param N.T  sample size in Treatment group
#' @param N.C  sample size in Control group
#' @param Delta  constraint for P.C-P.T
#' @return a vector (MR.T, MR.C) where MR.T and MR.C are the restricted MLE for the Treatment and Control group respectively
#' @export
restricted.ml = function(x.T, x.C, N.T, N.C, Delta) {
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



#' Test Statistic for Non-inferiority Hypothesis
#'
#' This function  calculate statistics for testing non-inferiority based on (Santner & Snell (1980), Blackwelder (1982),
#' Miettinen & Nurminen (1985) and Farrington & Manning (1990))
#'
#' @param x.T  observed number of responders in Treatment group
#' @param x.C  observed number of responders in Control group
#' @param N.T  sample size in Treatment group
#' @param N.C  sample size in Control group
#' @param Delta0  non-inferiority margin
#' @param method method for ordering criterion
#' @return statistic value of ordering criterion
#' @export
orderfun = function(x.T, x.C, N.T, N.C, Delta0, method) {
  rml = restricted.ml(
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

  if (method == "all") {
    out = list("MN" = mn.res, "FM" = fm.res, "SS" = ss.res, "Blackwelder" = blackwelder.res)
  }
  out
}



#' Order Statistic for all 2x2 Tables
#'
#' This function calculate order statistic for all 2x2 tables given N.T, N.C and an ordering criterion
#'
#' @param N.T  sample size in Treatment group
#' @param N.C  sample size in Control group
#' @param Delta0 non-inferiority margin
#' @param method method for ordering criterion
#' @return (N.T+1)x(N.C+1) array where the (i,j) element is the order statistic for x.T=i and x.C=j
#' @export
order.mat = function(N.T, N.C, Delta0, method) {
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
#' This function checks whether Barnard's Criterion is verified
#'
#' @param mat (N.T+1)x(N.C+1) array where the (i,j) element is the order statistic for x.T=i and x.C=j
#' @return TRUE if Barnard's Criterion is satisfied, FALSE otherwise
#' @export
barnard.check=function(mat){
  N.T=dim(mat)[1]-1
  N.C=dim(mat)[2]-1

  out=TRUE
  for (i in 1:(N.T+1)) {
    for (j in 1:(N.C+1)) {

      if((i+1)<=(N.T+1)) {out=out*(mat[i,j]<=mat[i+1,j])}
      if((j+1)<=(N.C+1)) {out=out*(mat[i,j]>=mat[i,j+1])}
    }}
  out
}

#' Chan's Exact P-value
#'
#' This function calculates the exact p-value based on Chan (1998)
#'
#' @param x.T  observed number of responders in Treatment group
#' @param x.C  observed number of responders in Control group
#' @param N.T  sample size in Treatment group
#' @param N.C  sample size in Control group
#' @param Delta0  non-inferiority margin
#' @param method method for ordering criterion
#' @param lower TRUE for the null hypothesis P.T-P.C<=-Delta0, FALSE for the null hypothesis P.T-P.C>-Delta0
#' @param tol increment size for domain of Delta, default set to 0.001
#' @return Exact p-value
#' @export
chan.pval <- function(x.T, x.C, N.T, N.C, Delta0, method, lower = TRUE,tol=1e-3) {
  myorder.mat = order.mat(N.T, N.C, Delta0 = Delta0, method = method)
  obs = myorder.mat[x.T + 1, x.C + 1]

  P.T.vec = seq(max(0,-Delta0), min(1, 1 - Delta0), by = tol)
  pvals = rep(NA, length(P.T.vec))
  for (ptvar in 1:length(P.T.vec)) {
    P.T = P.T.vec[ptvar]
    P.C = P.T + Delta0

    if(lower==TRUE){
      pvals[ptvar] = sum((myorder.mat>=obs) * dbinom(0:N.T, N.T, P.T) %o%dbinom(0:N.C, N.C, P.C))
    }
    if(lower==FALSE){
      pvals[ptvar] = sum((myorder.mat<=obs) * dbinom(0:N.T, N.T, P.T) %o%dbinom(0:N.C, N.C, P.C))

    }
  }
  max(pvals)

}

#' Chan p-value for All 2x2 Tables
#'
#' This function calculates the Chan p-values for all 2x2 tables given N.T, N.C and an ordering criterion
#'
#' @param N.T  Sample size in Treatment group
#' @param N.C  Sample size in Control group
#' @param Delta0  Non-inferiority margin
#' @param method method for ordering criterion
#' @return (N.T+1)x(N.C+1) array where the (i,j) element is the exact p-value for x.T=i and x.C=j
#' @export
#'
chan.mat = function(N.T, N.C, Delta0, method) {
  mat = array(NA, dim = c(N.T + 1, N.C + 1))
  for (i1 in 0:(N.T)) {
    for (i2 in 0:(N.C)) {
      mat[i1 + 1, i2 + 1] = chan.pval(
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
#' @param x.T  observed number of responders in Treatment group
#' @param x.C  observed number of responders in Control group
#' @param N.T  sample size in Treatment group
#' @param N.C  sample size in Control group
#' @param Delta  constraint for P.C-P.T
#' @param Delta0  non-inferiority margin
#' @return Exact-Corrected statistic
#' @export
orderfunEC = function(x.T, x.C, N.T, N.C, Delta, Delta0) {
  ECval=0
  rml = restricted.ml(x.T = x.T, x.C = x.C, N.T = N.T, N.C = N.C, Delta = Delta0)
  MR.T = rml$MR.T
  MR.C = rml$MR.C
  DEN_obs <- (MR.T * (1 - MR.T) / N.T + MR.C * (1 - MR.C) / N.C)
  pval=chan.pval(x.T=x.T,x.C=x.C,N.T=N.T,N.C=N.C,Delta0=Delta0,method="MN",lower=TRUE)
  if(pval>1e-10 & pval<(1-1e-10)) {
    ECval=sqrt(DEN_obs)*(orderfun(x.T=x.T, x.C=x.C, N.T=N.T, N.C=N.C,Delta0=Delta0,method="MN")-qnorm(1-pval))
  }

  rml = restricted.ml(x.T = x.T, x.C = x.C, N.T = N.T, N.C = N.C, Delta = Delta)
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
#' @param x.T  observed number of responders in Treatment group
#' @param x.C  observed number of responders in Control group
#' @param N.T  sample size in Treatment group
#' @param N.C  sample size in Control group
#' @param P.T  proportion of responders in Treatment group
#' @param Delta0  non-inferiority margin
#' @return Exact p-value
#' @export
likelihood.null = function(x.T, x.C, N.T, N.C, P.T, Delta0) {
  P.C = P.T + Delta0
  dbinom(x.T, N.T, P.T) * dbinom(x.C, N.C, P.C)
}

#' Level of Chan's Exact p-value
#'
#' This function calculates the level of the exact p-value based on Chan (1998)
#'
#' @param alpha  significance level
#' @param N.T  sample size in Treatment group
#' @param N.C  sample size in Control group
#' @param Delta0  non-inferiority margin
#' @param method method for ordering criterion
#' @return True level of Chan's p-value
#' @export
chan.level = function(alpha, N.T, N.C, Delta0, method) {
  mat = chan.mat(
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
      tot = tot + likelihood.null(
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
#' @param x.T  observed number of responders in Treatment group
#' @param x.C  observed number of responders in Control group
#' @param N.T  sample size in Treatment group
#' @param N.C  sample size in Control group
#' @param Delta0  non-inferiority margin
#' @param method method for ordering criterion
#' @param EC TRUE if confidence limits are exact-corrected to match Chan's exact p-value
#' @param alpha significance level
#' @param tol tolerance for convergence
#' @return List of lower and upper confidence limits (D.lower, D.upper)
#' @export
confint.z = function(x.T, x.C, N.T, N.C, Delta0, method, EC, alpha=.05, tol=1e-10) {
  ECval=0
  count=0
  if(EC==TRUE){
    rml = restricted.ml(x.T = x.T, x.C = x.C, N.T = N.T, N.C = N.C, Delta = Delta0)
    MR.T = rml$MR.T
    MR.C = rml$MR.C
    DEN_obs <- (MR.T * (1 - MR.T) / N.T + MR.C * (1 - MR.C) / N.C)
    pval=chan.pval(x.T=x.T,x.C=x.C,N.T=N.T,N.C=N.C,Delta0=Delta0,method=method,lower = TRUE)
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
      rml = restricted.ml(x.T = x.T, x.C = x.C, N.T = N.T, N.C = N.C, Delta = -Dmid)
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
      rml = restricted.ml(x.T = x.T, x.C = x.C, N.T = N.T, N.C = N.C, Delta = -Dmid)
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

  list(D.lower=D.lower,D.upper=D.upper)
}


#' Chan and Zhang Confidence Interval
#'
#' This function calculates the confidence interval based on Chan and Zhang (1999)
#'
#' @param x.T  observed number of responders in Treatment group
#' @param x.C  observed number of responders in Control group
#' @param N.T  sample size in Treatment group
#' @param N.C  sample size in Control group
#' @param Delta0  non-inferiority margin
#' @param method method for ordering criterion
#' @param alpha significance level
#' @param tol increment size for domain of Delta, default set to 0.001
#' @param width range from starting values based on Miettinen & Nurminen confidence limits
#' @return Chan and Zhang Confidence Interval
#' @export
chan.zhang <- function(x.T, x.C, N.T, N.C, method, alpha=.05, tol=1e-3, width=0.3) {
  d_LL = confint.z(x.T, x.C, N.T, N.C, Delta0=0, method,EC=F)$D.lower
  d_UL = confint.z(x.T, x.C, N.T, N.C, Delta0=0, method,EC=F)$D.upper

  deltaL=seq(max(-d_LL-width,-0.9999),min(-d_LL+width,0.9999),tol)
  deltaU=seq(max(-d_UL-width,-0.9999),min(-d_UL+width,0.9999),tol)
  probsL=rep(NA, length(deltaL))
  probsU=rep(NA, length(deltaL))

    for (i in 1:length(deltaL)) probsL[i]=chan.pval(x.T, x.C, N.T, N.C, deltaL[i], method ="MN", lower=T)
    for (i in 1:length(deltaU)) probsU[i]=chan.pval(x.T, x.C, N.T, N.C, deltaU[i],method = "MN", lower=F)

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
#' @param alpha  significance level
#' @param N.T  sample size in Treatment group
#' @param N.C  sample size in Control group
#' @param Delta0  non-inferiority margin
#' @param method method for ordering criterion
#' @param EC TRUE if confidence limits are exact-corrected to match Chan's exact p-value, only relevant if CZ=F
#' @param tolEC tolerance for convergence
#' @param CZ TRUE if Chan and Zhang method used, FALSE otherwise
#' @param tolCZ increment size for domain of Delta, default set to 0.001
#' @param width range from starting values based on Miettinen & Nurminen confidence limits
#' @return Level of test based on the specified confidence interval
#' @export
ci.level <- function(alpha, N.T, N.C, Delta0, method, EC, tolEC, CZ, tolCZ, width) {
  M=matrix(0,N.T+1,N.C+1)
  for(i in 0:N.T) {
    for(j in 0:N.C) {
      if(CZ) {
        ci = chan.zhang(x.T=i,x.C=j,N.T=N.T,N.C=N.C,method=method,tol=tolCZ,width=width,alpha=alpha)
        lb=ci$D.lower
        ub=ci$D.upper
        count=ci$count
      }
      if(!CZ) {
        ci = confint.z(x.T=i,x.C=j,N.T=N.T,N.C=N.C,Delta0=Delta0,method=method,EC=EC,alpha=alpha,tol=tolEC)
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
        tot = tot + likelihood.null(
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

