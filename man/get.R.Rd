\name{get.R}
\alias{get.R}

\title{
Get Correlation matrix from a Covariance matrix
}
\description{
Returns a correlation matrix from a variance-covariance matrix
}
\usage{
get.R(Sigma0)
}

\arguments{
  \item{Sigma0}{
variance-covariance matrix of dimensions kxk
}
}

\value{
\item{R }{correlation matrix correspondig to Sigma0}
%%  \item{comp2 }{Description of 'comp2'}
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
k<-6 # variables
B<-matrix(runif(n = k*k),ncol= k)### creating random matrix for sigma
sigma = B\%*\%t(B)
get.R(Sigma0=sigma)
}
