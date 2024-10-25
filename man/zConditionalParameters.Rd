\name{zConditionalParameters}
\alias{zConditionalParameters}

\title{
Conditional parameters for z, given x
}
\description{
This function calculates and returns conditional parameters for z, given x are being already considered in the model
}
\usage{
zConditionalParameters(mean0, R0, z, x.var, z.var)
}

\arguments{
  \item{mean0}{
Mean vector for multivariate random vector under the null hypothesis. Dimensions: kx1
}
  \item{R0}{
 Correlations matrix for multivariate random vector under the null hypothesis. Dimensions kxk

}
  \item{z}{
vector of random observation. Dimensions kx1
}
  \item{x.var}{
Elements of z that are already considered in the model
}
  \item{z.var}{
element of z whose contribution to C_k|C_k-1,C_k-2,...,C_1 is going to be calculated
}
}

\value{
A list containing
\item{muC }{conditional mean for z}
\item{RC}{Conditional variance for z}
%% ...
}
\references{
Paper
}
\author{
Dr. Burcu Aytaçoğlu (burcuaytacoglu@gmail.com)
Dr. Diana Barraza-Barraza (diana.barraza@ujed.mx),
Dr. Víctor G. Tercero-Gómez (victor.tercero@tec.mx),
Dr. A. Eduardo Cordero-Franco (lalo.cordero@gmail.com),

}

\examples{
        k<-3
      sigma0 = matrix(diag(rep(1,k)),ncol = k)
      mu0 = matrix(c(0,0,0), ncol = 1)
      Weights = diag(c(0.5, 0.25,0.25))

      library(mvtnorm)
      set.seed(1000)
      X = matrix(ncol= 1, data = rmvnorm(n = 1, mean = mu0, sigma = sigma0))
      Z = (X - mu0)/sqrt(as.numeric(diag(sigma0)))
      Corr<-get.R(Sigma0 = sigma0)

      zConditionalParameters(mean0 = mu0, R0 = Corr, z = Z, x.var = 1:2, z.var = 3)

}
