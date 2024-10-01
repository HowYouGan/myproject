library(profvis)

compute_Xty = function (X, y) {
  cm = attr(X,"scaled:center")
  csd = attr(X,"scaled:scale")
  ytX = crossprod(y,X)

  # Scale Xty.
  if (!is.null(attr(X,"matrix.type")))

    # When X is a trend filtering matrix.
    scaled.Xty = compute_tf_Xty(attr(X,"order"),y)/csd
  else

    # When X is an ordinary sparse/dense matrix.
    scaled.Xty = t(ytX/csd)

  # Center Xty.
  centered.scaled.Xty = scaled.Xty - cm/csd * sum(y)
  return(as.numeric(centered.scaled.Xty))
}

single_effect_regression =
  function (y, X, V, residual_variance = 1, prior_weights = NULL,
            optimize_V = c("none", "optim", "uniroot", "EM", "simple"),
            check_null_threshold = 0) {
    optimize_V = match.arg(optimize_V)
    Xty = compute_Xty(X,y)
    betahat = (1/attr(X,"d")) * Xty
    shat2 = residual_variance/attr(X,"d")
    if (is.null(prior_weights))
      prior_weights = rep(1/ncol(X),ncol(X))
    if (optimize_V != "EM" && optimize_V != "none")
      V = optimize_prior_variance(optimize_V,betahat,shat2,prior_weights,
                                  alpha = NULL,post_mean2 = NULL,V_init = V,
                                  check_null_threshold = check_null_threshold)

    # log(bf) for each SNP
    lbf = dnorm(betahat,0,sqrt(V + shat2),log = TRUE) -
      dnorm(betahat,0,sqrt(shat2),log = TRUE)

    # Deal with special case of infinite shat2 (e.g., happens if X does
    # not vary).
    lbf[is.infinite(shat2)] = 0
    maxlbf = max(lbf)

    # w is proportional to BF, but subtract max for numerical stability.
    w = exp(lbf - maxlbf)

    # Posterior prob for each SNP.
    w_weighted = w * prior_weights
    weighted_sum_w = sum(w_weighted)
    alpha = w_weighted / weighted_sum_w
    post_var = (1/V + attr(X,"d")/residual_variance)^(-1) # Posterior variance.
    post_mean = (1/residual_variance) * post_var * Xty
    post_mean2 = post_var + post_mean^2 # Second moment.

    # BF for single effect model.
    lbf_model = maxlbf + log(weighted_sum_w)
    loglik = lbf_model + sum(dnorm(y,0,sqrt(residual_variance),log = TRUE))

    if(optimize_V == "EM")
      V = optimize_prior_variance(optimize_V,betahat,shat2,prior_weights,
                                  alpha,post_mean2,
                                  check_null_threshold = check_null_threshold)

    return(list(alpha = alpha,mu = post_mean,mu2 = post_mean2,lbf = lbf,
                lbf_model = lbf_model,V = V,loglik = loglik))
  }
