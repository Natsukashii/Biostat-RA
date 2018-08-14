my_gendata <-
  function(n,p,k.lin,k.nonlin,deg,sigma, fixed_X, fixed = FALSE){
    alpha0 <- 0.05
    if (!fixed){
      alpha <- c(matrix(10 + 3*rnorm(k.lin), ncol=1), rep(0,p-k.lin))
      beta <- matrix(c(rep(0,k.lin*deg),20*rnorm(deg*k.nonlin), rep(0,deg*(p-k.nonlin-k.lin))),nrow=deg*p, ncol=1)
    } else{
      alpha <- c(matrix(rep(13, k.lin), ncol=1), rep(0,p-k.lin))
      beta <- matrix(c(rep(0,k.lin*deg),rep(beta_nonlinear, k.nonlin), rep(0,deg*(p-k.nonlin-k.lin))),nrow=deg*p, ncol=1)
    }
    if (missing(fixed_X)){
      X = scale(matrix(rnorm(n*p), n, p))
    } else{
      X = fixed_X
    }
    U <- NULL
    for(j in 1:p) {
      Uj <- matrix(poly(X[,j],deg), ncol=deg)
      U <- cbind(U, Uj)
      X[,j] <- Uj[,1]
    }
    truef = alpha0 + X%*%alpha + U%*%beta
    y <- truef + sigma*rnorm(n)
    probs <- 1/(1 + exp(-truef))
    yb <- rbinom(n,1,probs)
    list(X=X,y=y,yb=yb, alpha=alpha, beta=beta, U=U, truef=truef)
  }