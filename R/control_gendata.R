control_gendata <-
  function(n,p,k.lin,k.nonlin,deg,sigma, fixed_X, fixed_beta){
    alpha0 <- 0.05
    alpha <- c(matrix(rnorm(k.lin), ncol=1), rep(0,p-k.lin))
    if (missing(fixed_beta)){
      beta <- matrix(c(rep(0,k.lin*deg),5*rnorm(deg*k.nonlin), rep(0,deg*(p-k.nonlin-k.lin))),nrow=deg*p, ncol=1)
    } else{
      beta <- matrix(fixed_beta, ncol = 1)
    }
    if (missing(fixed_X)){
      X = scale(matrix(rnorm(n*p), n, p))
    } else{
      X = fixed_X
    }
    U <- NULL
    for(j in 1:p) {
      Uj <- matrix(poly(X[,j],raw = T, deg), ncol=deg)
      U <- cbind(U, Uj)
    }
    y_mat <- NULL
    for(k in 1:p) {
      yk = U[,((k-1)*deg+1):(k*deg)] %*% beta[((k-1)*deg+1):(k*deg)] + X[,k] * alpha[k]
      y_mat <- cbind(y_mat, yk)
    }
    y_mat = scale(y_mat, scale = F)
    truef = alpha0 + rowSums(y_mat)
    y <- truef + sigma*rnorm(n)
    probs <- 1/(1 + exp(-truef))
    yb <- rbinom(n,1,probs)
    list(X=X,y=y,yb=yb, y_mat = y_mat, alpha=alpha, beta=beta, U=U, truef=truef)
  }

