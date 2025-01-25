\name{p2}
\alias{p2}

\title{
  MCMC estimates for the p2 model
}
\description{
  Estimates p2 model parameters with the adaptive random walk algorithm as described in Zijlstra, Van Duijn and Snijders (2009) <doi: 10.1348/000711007X255336>.
}
\usage{
p2(net, sender = NULL, receiver = NULL, density = NULL, reciprocity = NULL, 
burnin = NULL, sample = NULL, adapt = NULL, seed = NULL)
}

\arguments{
  \item{net}{Directed dichotomous n*n network (digraph).}
  \item{sender}{Optional sender covariates of lenght n.}
  \item{receiver}{Optinal receiver covariates of length n.}
  \item{density}{Optional density covariates of dimensions n*n.}
  \item{reciprocity}{Optional symmetric reciprocity covariates of dimensions n*n.}
  \item{burnin}{Optional specification of number of burn-in iterations (default is 10000).}
  \item{sample}{Optional specification of number of MCMC samples (default is 40000).}
  \item{adapt}{Optional number of adaptive sequenses (default is 100).}
  \item{seed}{Optional specification of random seed (delfault is 1).}
}

\value{
Returns a matrix with MCMC means, standard deviations, quantiles and estimated effective sample sizes for p2 parameters.
}
\references{
Zijlstra, B.J.H., Duijn, M.A.J. van, and  Snijders, T.A.B. (2009). MCMC estimation for the $p_2$ network regression model with crossed random effects. \emph{British Journal of Mathematical and Statistical Psychology}, 62, 143-166.
}
\author{
Bonne J.H. Zijlstra \email{b.j.h.zijlstra@uva.nl}
}

\examples{

# create a very small network with covariates for illustrative purposes
S <- c(1,0,1,0,1,1,0,1,0,1)
REC <- (S*-1)+1
D1 <- matrix(c(0,1,0,1,0,1,0,1,0,1,
              0,0,0,1,0,1,0,1,0,1,
              1,1,0,0,1,0,0,0,0,0,
              1,1,1,0,1,0,0,0,0,1,
              1,0,1,0,0,1,1,0,1,0,
              0,0,0,0,0,0,1,1,1,1,
              0,0,0,0,0,1,0,1,0,1,
              1,0,0,0,0,1,1,0,1,1,
              0,1,0,1,0,1,0,1,0,0,
              1,0,1,1,1,0,0,0,0,0), ncol=10)
D2 <- abs(matrix(rep(S,10), byrow = FALSE, ncol= 10) -
            matrix(rep(REC,10), byrow = TRUE, ncol= 10))
R <- D1*t(D1)
Y <- matrix(c(0,1,1,1,1,1,0,0,1,1,
              0,0,0,1,1,1,0,0,1,0,
              1,1,0,1,1,1,0,0,1,1,
              1,1,1,0,1,1,0,1,1,0,
              1,1,1,1,0,1,1,0,1,1,
              0,1,1,1,1,0,1,1,1,0,
              1,0,1,0,1,1,0,1,0,1,
              0,1,1,1,0,1,1,0,1,1,
              1,0,1,0,1,0,1,1,0,1,
              1,1,1,0,0,1,1,1,1,0), ncol=10) 

# estimate p2 model
p2(Y,sender= ~ S, receiver =  ~ REC, density = ~ D1 + D2, reciprocity= ~ R,
   burnin = 100, sample = 400, adapt = 10)
# Notice: burn-in, sample size and number of adaptive sequenses are 
# much smaller than recommended to keep computation time low.
# recommended code: 
\dontrun{
p2(Y,sender= ~ S, receiver =  ~ REC, density = ~ D1+ D2, reciprocity= ~ R)
}
}
\keyword{p2 model}

