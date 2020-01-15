# SampleID...TimeSeries...Value
# 52-001-201...01...23
# 52-001-201...02...NA
# 52-001-201...03...30
# the test for github desktop

LPRP = function(Z, x, p, h){
  #Local polynomial regression estimator with Gaussian kernel
  #the first column of Z is covariates; second column is response
  #x is the point to predict
  X = Z[, 1]; Y = Z[, 2]
  n = length(X)
  K = c()
  for (i in 1 : n){
    K[i] = dnorm((X[i] - x) / h)
  }
  W = diag(K / h)
  design = matrix(0, n, p + 1)
  for (j in 0 : p){
    design[, j + 1] = (X - x)^j
  }
  D = t(design) %*% W %*% design
  if (abs(det(D)) <= 10^(-9)) cat("Error: try smaller p or larger h", "\n")
  if (abs(det(D)) > 10^(-9)){
    beta = solve(t(design) %*% W %*% design) %*% t(design) %*% W %*% matrix(Y, n, 1)
    return(beta)
  }
}

Fcurve <- function(Z, x, h) {
  # h is a vector of bandwidth for prediction of response and its derivatives
  # y is the predicted value, d1 is the first derivative, d2 is the second derivative, d3 is the third derivative
  y = LPRP(Z, x, 1, h[1])[1]
  d1 = LPRP(Z, x, 2, h[2])[2]
  d2 = 2 * LPRP(Z, x, 3, h[3])[3]
  d3 = 3 * 2 * LPRP(Z, x, 4, h[4])[4]
  return(c(y, d1, d2, d3))
}
