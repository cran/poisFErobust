pois.fe.robust <-
  function(outcome, xvars, group.name, data, 
           qcmle.coefs = NULL, allow.set.key = FALSE, index.name = NULL){
    if (!is.null(index.name)){
      warning("index.name is deprecated. It is not necessary.")
    }
    # If coefficients are not supplied, then estimate them.
    if (is.null(qcmle.coefs)){
      qcmle.fit <- glmmML::glmmboot(reformulate(xvars, response = outcome),
                                    family = "poisson", data = data,
                                    cluster = get(group.name), boot = 0)
      qcmle.coefs <- qcmle.fit$coefficients
    }
    if (!is.data.table(data)){
      stop(paste("data must be a data.table. Use as.data.table() to coerce.",
                 "For example, data <- as.data.table(data)."))
    }
    if (!is.factor(data[[group.name]])){
      stop("group.name must be a factor.")
    }
    if (allow.set.key == FALSE){
      data <- copy(data)
      warning("allow.set.key is FALSE. TRUE is recommended to avoid copying ",
              "the data; however, TRUE will sort the data in-place by ",
              "group.name.")
    }
    setkeyv(data, c(group.name))
    
    # Estimate robust standard errors. 
    # Naming convention aligns with Wooldridge (1999).
    p <- length(qcmle.coefs)
    id.rle <- rle(as.integer(data[[group.name]]))
    bigN <- length(id.rle$values)
    outlist <- vector("list", length = bigN)
    
    bigK <- matrix(0, p, p)
    bigA <- matrix(0, p, p)
    bigB <- matrix(0, p, p)
    
    xvarcols <- match(xvars, colnames(data))
    ycol <- match(outcome, colnames(data))
    
    if (anyNA(data[, c(ycol, xvarcols, match(group.name, colnames(data))),
                   with = FALSE])){
      stop("NA detected. Please remove NA observations.")
    }
    
    starti <- 0L
    endi <- 0L
    
    for(i.index in 1:bigN) {
      starti <- endi + 1
      endi <- endi + id.rle$lengths[i.index]
      indicies <- starti:endi
      y <- data[indicies, ycol, with = FALSE]
      n_i <- sum(y)
      xvarcolmatrix <- as.matrix(data[indicies, xvarcols, with = FALSE]) # TxP
      this.expxb.vec <- exp(xvarcolmatrix %*% qcmle.coefs) # Tx1
      sumfor_it <- sum(this.expxb.vec)
      p_i <- this.expxb.vec / sumfor_it
      outlam <- (1/sumfor_it^2)*(sumfor_it * as.vector(this.expxb.vec) * 
                 xvarcolmatrix - 
                 this.expxb.vec %*% colSums(as.vector(this.expxb.vec) *
                                              xvarcolmatrix))
      W_i <- diag(as.vector(1/p_i))
      u_i <- as.matrix(y) - n_i * p_i
      outlist[[i.index]] <- list("L_i" = outlam, "W_i" = W_i, "u_i" = u_i)
      
      bigK <- bigK + 1/bigN * (t(outlam)) %*% (n_i * outlam)
      bigA <- bigA + 1/bigN * n_i * (t(outlam) %*% W_i %*% outlam)
      bigB <- bigB + 1/bigN * t(outlam) %*% W_i %*% u_i %*% t(u_i) %*% W_i %*% 
        outlam
    }
    
    bigAinv <- solve(bigA)
    se.rob <- sqrt(diag((bigAinv %*% bigB %*% bigAinv)/ bigN))
    names(se.rob) <- xvars
    
    rmat <- matrix(0, bigN, p)
    for (i in 1:bigN){
      rmat[i,] <- t(outlist[[i]]$u_i) %*% 
        (outlist[[i]]$L_i - outlist[[i]]$W_i %*% 
           outlist[[i]]$L_i %*% bigAinv %*% t(bigK))
    }
    
    ssr <- sum(.lm.fit(x = rmat, y = rep(1, nrow(rmat)))$residuals^2)
    p.val.cond.mean <- 1 - pchisq(q = bigN - ssr, df = p)
    
    return(list(coefficients = qcmle.coefs, se.robust = se.rob, 
                p.value = p.val.cond.mean))
  }
