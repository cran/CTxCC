\name{SimulatedDistributionC2}
\alias{SimulatedDistributionC2}
\title{
Distribution for C2, through simulation of its values
}
\description{
Simulates s instances of C^2_k given 1 to k-1 variables are already in the model. Obtains the quantile indicated by alpha
}
\usage{
SimulatedDistributionC2(z, R.C, mu.C, W, R, A, x.var, z.var, alpha, s)
}

\arguments{
  \item{z}{
observation vector, kx1
}
  \item{R.C}{
scalar, conditional covariance for z given x
}
  \item{mu.C}{
scalar, conditional mean for z given x
}
  \item{W}{
matrix of variables weigths, kxk
}
  \item{R}{
correlation matrix, kxk
}
  \item{A}{
list containing matrix decomposition of A, preferably, obtained from function decomposeA
}
  \item{x.var}{
vector indicating variables already present in the model. length: k-1
}
  \item{z.var}{
scalar indicating variable to be included
}
  \item{alpha}{
quantile(s) of the distribution
}
  \item{s}{
scalar indicating amount of simulations
}
}

\value{
Quantile(s) of the simulated distribution
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

A<-decomposeA(W = Weights, R = Corr, x.var = 1:2, z.var = 3)

Par<-zConditionalParameters(mean0 = mu0, R0 = Corr, z = Z, x.var = 1:2, z.var = 3)
SimulatedDistributionC2(z = Z, R.C = Par$RC, mu.C = Par$muC, W = Weights, R = Corr,
                        A = A, x.var = 1:2, z.var = Z, alpha = 0.95, s = 1000 )
}
