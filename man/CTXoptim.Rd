\name{CTXoptim}
\alias{CTXoptim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimum Design of a CTX Control Chart
}
\description{
Computes the optimal weight vector and control limit for a CTX control chart that
minimizes the overall out-of-control Average Run Length (ARL) for a given set of
shift vectors. The optimization can be performed using either a single-start or a
multistart approach, depending on the number of variables in the dataset.
For the single-start case, the BFGS method is used, while in the multistart case,
the optimization is performed using either BFGS or Nelder–Mead.}

\usage{
CTXoptim(delta, R, alpha,
n_starts, method_primary,
method_secondary, scale_sum,jitter_sd, seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{Numeric. Shift size or shift vector considered in the optimization.}
  \item{R}{Numeric. Correlation matrix or covariance structure used in the CTX control chart.}
  \item{alpha}{Numeric. Significance level (Type I error rate) for the control chart.}
  \item{n_starts}{Integer. Number of starting points for the optimization. Default is \code{30}.}
  \item{method_primary}{Character. Optimization method for the single-start case. Default is \code{"BFGS"}.}
  \item{method_secondary}{Character. Optimization method for the multistart case (alternative to primary). Default is \code{"Nelder-Mead"}.}
  \item{scale_sum}{Numeric. Scaling constant applied to the weight vector. Default is \code{1}.}
  \item{jitter_sd}{Numeric. Standard deviation of the jitter added when generating starting points in multistart optimization. Default is \code{0.3}.}
  \item{seed}{Integer. Random seed for reproducibility in multistart optimization. Default is \code{123}.}
}

\details{
The optimization is performed over proxy weights (log-transformed) to ensure positivity, and the weights are normalized to sum to one. The user can choose between single-start and multistart optimization, depending on the number of variables (and out-of-control ARLs) to be monitored. A quasi-Newton method (\code{BFGS}) is used to minimize the sum of out-of-control ARLs across the specified shift scenarios. In the multistart case, if BFGS fails, the Nelder–Mead algorithm is applied instead. The control limit is computed to satisfy the desired in-control ARL, with \code{alpha = 1/ARL0}.
}
\value{
    \item{weights}{The optimized weight vector. Sum up to scale_sum}
    \item{control_limit}{The computed control limit.}
    \item{objective}{The optimal value of the objective function.}
    \item{method}{The optimization method used.}
    \item{n_starts}{The number of starting points used in the optimization.}

}

\references{
Paper (pending complete bibliographic reference when available).

}
\author{
Dr. Burcu Aytaçoğlu (\email{burcuaytacoglu@gmail.com}),
Dr. Diana Barraza-Barraza (\email{diana.barraza@ujed.mx}),
Dr. Víctor G. Tercero-Gómez (\email{victor.tercero@tec.mx}),
Dr. A. Eduardo Cordero-Franco (\email{lalo.cordero@gmail.com})
}

\examples{
\dontrun{
# Define shift scenarios
delta1 <- c(1, 0.5)
delta <- list(delta1)
rho <- 0.75

 # Correlation matrix and false alarm rate
 R <- diag(2); R[1,2] <- R[2,1] <- rho
alpha <- 1/200

 # Compute optimal weights and control limit
 CTXoptim(delta = delta, R = R, alpha = alpha,
                       n_starts = 30, scale_sum = 1, seed = 123)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
