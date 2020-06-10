#' You can use this r code function urt to test for the existence of unit roots
#' @title ur.test
#' @details
#' @param Z input one dimensional time series
#' @param K0_vec lag K0 vector, if missing, the default value we choose K0_vec=c(0,1,2,3,4)
#' @param con_vec constant vector for ck, if missing, the default value we use 0.55
#' @param alpha significance level
#' @return A dataframe, "1" means we reject the null hypothesis and "0" means we do not reject the null hypothesis
#' @export
#' @import tseries
#' @import sandwich
#' @examples
#' N=100
#' Z=arima.sim(list(ar=c(0.9)), n = 2*N, sd=sqrt(1))
#' con_vec=c(0.45,0.55,0.65)
#' K0_vec=c(0,1,2)
#' ur.test(Z,K0_vec=K0_vec, con_vec=con_vec,alpha=0.05)
#' ur.test(Z,alpha=0.05)


ur.test <- function(Z, K0_vec=K0_vec, con_vec=con_vec, alpha=alpha) {

     args = as.list(match.call())
     if(is.null(args$K0_vec)){
      K0_vec=c(0,1,2,3,4);
     }
     args = as.list(match.call())
     if(is.null(args$con_vec)){
       con_vec=0.55;
     }
     require(tseries)                      ## kpss package
     require(sandwich)                     ## long-run variance package

#  if(con_vec[1]!=0){
#  	print("Untruncated needs constant equals to 0");  break ;             ## untruncated needs con_vec equals to 0
#  }
  Tnvec=NULL; nm=NULL;
for(kk in 1:length(K0_vec)){

  K0=K0_vec[kk]+1                       #eg. K0=1, gamma(0)
  nm=c(nm,paste("K0=", K0-1, sep=""))

  n=length(Z)                           ## sample size
  N=floor(n/2)
  N1=2*N-K0

  sgn_matrix=matrix(0,N1,K0)            ### sign matrix
  for(t in 1:N1){
    for(k in 1:K0){
      sgn_matrix[t,k]=sign(k+t-N-1-0.5) ##  eg. K0=1, gamma(0)
    }
  }

  Y=Z;  DY=diff(Y) ## diffential Z

  au_Y =drop(acf(Y,lag.max =K0+1, type = c("covariance"),plot = FALSE)$acf)   ## gamma(Y)
  au_DY=drop(acf(DY,lag.max=K0+1, type = c("covariance"),plot = FALSE)$acf)   ## gamma(X)

  short_Var=var(DY)                                         ## shortrun variance
  long_Var=n*lrvar(DY, type = "andrews", prewhite = FALSE)  ## longrun variance
  ratio_Var=short_Var/long_Var                              ## variance ratio


  ## esttimate rho
  Z2=Z[1:(n-1)];  Z1=Z[2:n]
  DZ2=diff(Z2);   DZ1=diff(Z1)
  rho_hat=lm(DZ2~DZ1)$coefficients[2]; bb=1+rho_hat         ## rho_hat

  au_Ratio=(au_Y[1]+au_Y[2])/(au_DY[1]+au_DY[2])            ## ratio

  Y1=Y[1:N]                                   ##data spliting
  Y2=Y[(N+1):(2*N)]

  ## auto covariance
  auto_cov=drop(acf(Y,lag.max =K0+1,type = c("covariance"),plot = FALSE)$acf)
  auto_cov1=drop(acf(Y1,lag.max =K0+1,type = c("covariance"),plot = FALSE)$acf)
  auto_cov2=drop(acf(Y2,lag.max =K0+1,type = c("covariance"),plot = FALSE)$acf)

  T1=sum((auto_cov1[1:(K0)])^2)
  T2=sum((auto_cov2[1:(K0)])^2)             ## test statistics

  ## construct Qt
  N1=2*N-K0; Y=Y-mean(Y)
  ft=Y[1:N1]%*%t(rep(1,K0));
  for(t in 1:N1){
    ft[t,]=ft[t,]*Y[(t):(t+K0-1)]                       ## gamma(0), .... gamma(K0-1) data
  }
  gamma_hat=t(auto_cov[1:(K0)]%*%t(rep(1,N1)))         ## estimate gamma(0),.... gamma(K0-1)
  ytk=2*(ft-gamma_hat)*sgn_matrix                       ## construct ytk
  sgn_Auto=t(sign(auto_cov1[1:(K0)])%*%t(rep(1,N1)))    ## sign auto covariance
  xitk=2*ytk*gamma_hat; Qt=apply(xitk,1,sum)               ## Qt
  lr_Qt=lrvar(Qt, type = "andrews", prewhite = FALSE)   ## long-run variance Qt

  ## test procedure

  kappa=2/(ratio_Var*bb);                          ## kappa
  if(lr_Qt>0){                                     ## variance >0
     #cv=qnorm(1-alpha)*sqrt(lr_Qt)+T1
     #Tnvec=c(Tnvec,T2>cv)                         ## no truncated
    for(tt in 1:length(con_vec)){
      ck=con_vec[tt]
      if (au_Ratio<=(ck*kappa*N^{3/5})){
        th_d=10^5                                  ## truncated belongs to H0
      }  else{th_d=0.1*log(N)}                     ## truncated belongs to H1

      cv=min(1.645*sqrt(lr_Qt)+T1, th_d)
      Tnvec=c(Tnvec, T2>cv)
    }
  }
}
    res.table=matrix(as.numeric(Tnvec), length(K0_vec), byrow=T)
    rownames(res.table)=nm        #rownames ("K0=1", "K0=2")
    return(list(result=res.table))

}

