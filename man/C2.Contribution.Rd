\name{C2.Contribution}
\alias{C2.Contribution}

\title{
Contribution of variable z.var to C^2
}
\description{
Returns contribution of variable z.var to C^2, even if there are no previous variables in the model
}
\usage{
C2.Contribution(z, mean0, W, R, x.var, z.var = NULL)
}

\arguments{
  \item{z}{
observation vector, kx1, where z[x.var, ] correspond to variables already in the model
}
  \item{mean0}{
Mean vector for multivariate random vector under the null hypothesis. Dimensions: kx1
}
  \item{W}{
matrix of variables weigths, kxk
}
  \item{R}{
correlation matrix, kxk
}
  \item{x.var}{
vector indicating variables already present in the model. length: k-1
}
  \item{z.var}{
scalar indicating variables to be included. Defaults to NULL, indicating there are no previous variables in the model
}
}

\value{
C2.k.extra, scalar containing the contribution of variable z to C_k
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

C2.Contribution(z = Z, W = Weights, R = Corr, x.var = 1:2, z.var = 3)

}
