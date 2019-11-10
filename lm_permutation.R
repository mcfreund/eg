## permuted p values for linear model ---

lm.permute <- function(df, yname, nresamples = 1E4, add.intercept = TRUE, alternative = "two.sided") {
  ## NB: NOT CONFIGURED FOR USE WHEN is.factor(df$column).
  ## categorical IVs can still be used, but the levels just need to be specified as separate integer columns
  ## (e.g., dummy variables).
  
  if (any(sapply(df, is.factor))) stop("NOT CONFIGURED FOR USE WITH FACTOR CLASS")
  
  ## get observed estimates ----
  
  xnames <- names(df)[-grep(yname, names(df))]
  lmform <- as.formula(paste0(yname, " ~ ", paste0(xnames, collapse = " + ")))
  
  observed <- as.data.frame(coef(summary(lm(lmform, df))))
  names(observed) <- c("b", "se", "t", "p")
  
  
  ## get permuted estimates ---- 
  
  ## convert to matrices (for speed)
  if (add.intercept) {
    X <- cbind(int = 1, as.matrix(df[xnames]))
  } else {
    X <- cbind(as.matrix(df[xnames]))
  }
  
  y <- as.matrix(df[yname])
  B <- matrix(NA, nrow = ncol(X), ncol = nresamples)  ## to hold all resamples
  
  resamples <- replicate(round(nresamples * 1.2), sample.int(nrow(y), nrow(y)))  ## oversample to remove duplicates
  resamples <- unique(resamples, MARGIN = 2)  ## remove duplicate permutations
  resamples <- resamples[, 1:nresamples]  ## downsample to nresamples
  
  ## fit models, but do it fast:
  for (ii in seq_len(nresamples)) B[, ii] <- c(solve(t(X) %*% X, t(X) %*% y[resamples[, ii]]))
  
  
  ## calculate p values ---- 
  
  ## proportion of permuted betas with...
  if (alternative == "two.sided") {
    ## ...more extreme estimates:
    observed$p.permuted <- rowSums(abs(B) >= abs(observed$b)) / nresamples
  } else if (alternative == "less") {
    ## ...smaller estimates:
    observed$p.permuted <- rowSums(B < observed$b) / nresamples
  } else if (alternative == "greater") {
    ## ...larger estimates:
    observed$p.permuted <- rowSums(B > observed$b) / nresamples
  }

  observed
  
}

df <- as.data.frame(matrix(rnorm(800), ncol = 8))
lm.permute(df, "V1")
