###### FUNCTIONS ########

get.R<-function(Sigma0){
  # Returns R matrix from variance-covariance matrix Sigma0
  #     Reference: Eq. (16) from Paper Criticality Assessment for
  #                Enhanced Multivariate Process Monitoring
  #     Author: Dr. Víctor G. Tercero-Gómez, Dr. Diana Barraza-Barraza,
  #             Dr. A. Eduardo Cordero-Franco, Dr. Burcu Aytaçoğlu
  #
  #     Date: October 6, 2023
  #     Versión: 1.0

  if(matrixcalc::is.positive.definite(Sigma0)){
    sigma0.dim = dim(Sigma0)
    dim.rows = sigma0.dim[1]
    dim.cols = sigma0.dim[2]

    R = matrix(rep(NA,dim.rows*dim.cols), ncol = dim.cols)
    for(ro in 1:dim.rows){
      for(co in 1:dim.cols){
        R[ro,co] = Sigma0[ro,co] / sqrt(as.numeric(diag(Sigma0))[ro]) / sqrt(as.numeric(diag(Sigma0))[co])
      }
    }
    return(R)
  }else{

    stop(paste0("Matrix Sigma0 is not positive definite."))
    return(NULL)
  }
}

decomposeA <-function(W,R,x.var, z.var){
  # Returns decomposition of matrix A
  #     Reference: Proposition Distribution of a C^2 contribution
  #                Eq. (41) from Paper Criticality Assessment for
  #                Enhanced Multivariate Process Monitoring
  #
  #     Author: Dr. Víctor G. Tercero-Gómez, Dr. Diana Barraza-Barraza,
  #             Dr. A. Eduardo Cordero-Franco, Dr. Burcu Aytaçoğlu
  #
  #     Date: October 6, 2023
  #     Versión: 1.0
  #   Parameteres
  #       W : matrix of variables weigths, kxk
  #       R : correlation matrix, kxk
  #       x.var: vector indicating variables already present in the model. length: k-1
  #       z.var: scalar indicating variables to be included.
  #       length(z.var) + length(x.var) = k

  k.prov<-c(x.var, z.var)
  dimW <- dim(W[k.prov,k.prov]) ## in order to check if R and W have the same dimmensions
  dimR <- dim(R[k.prov,k.prov])

  if(matrixcalc::is.square.matrix(W)&&matrixcalc::is.square.matrix(R)){ ### checks if matrices are square
    if(dimW[1]==dimR[1]){###checks if are of the same dimensions and
      if((length(x.var)+length(z.var))==dimW[1]){ ### checks that k = length(x.var)+length(z.var)
        ######## Proposition 4.3
        ### Eq. 41
        a.I<- W%*%solve(R)%*%t(W)#W[k.prov,k.prov]%*%solve(R[k.prov,k.prov])%*%t(W[k.prov,k.prov]);a.I #change "A" to "a" to avoid environment
        a<-solve(a.I);a
        axx<-a[x.var,x.var];axx
        azx<-matrix(a[x.var,z.var], ncol = 1);azx
        az2<-a[z.var,z.var];az2

        ### Eq 44
        a2zx <- az2-t(azx)%*%solve(axx)%*%azx
        a2zx
        return(list(Axx=axx, Azx = azx, Az2 = az2, A2zx=a2zx))
      }else{
        stop("W and R dimensions do not match with z.var and x.var lengths")
      }
    }else{
      stop("W and R dimensions do not match")
    }
  }else{
    stop("Please, check dimensions for W and R")
  }
}

zConditionalParameters<-function(mean0, R0, z, x.var, z.var){
  # Returns conditional parameters for z, given x already in model
  #     Reference: Proposition Distribution of a C^2 contribution
  #                Eq. (48) from Paper Criticality Assessment for
  #                Enhanced Multivariate Process Monitoring
  #
  #     Author: Dr. Víctor G. Tercero-Gómez, Dr. Diana Barraza-Barraza,
  #             Dr. A. Eduardo Cordero-Franco, Dr. Burcu Aytaçoğlu
  #
  #
  #     Date: October 6, 2023
  #     Versión: 1.0
  # R0 :  original correlation matrix for data distribution, kxk
  # z :   observation vector, kx1, where
  #        z[x.var, ] correspond to variables already in the model
  #        z[z.var, ] corresponds to new variable in the model
  # x.var: vector indicating variables already present in the model. length: k-1
  # z.var: scalar indicating variables to be included.
  # length(z.var) + length(x.var) = k

  k = length(z.var) + length(x.var)
  x.length<-length(x.var)
  dimR0<-dim(R0)


  # if(dimR0[1]==k){ ## checks dimensions for R0 and k
  #  if(NROW(z)==k){## checks dimensions for z and k, z should be a kx1 matrix
  z<-matrix(z, ncol=1)#nrow=k)

  # Calculation of z's conditional mean given x
  mu.C= 0+R0[x.var,z.var]%*%solve(R0[x.var,x.var])%*%(z[x.var,])

  # Calculation of z's conditional variance given x
  R.C = R0[z.var,z.var]-
    matrix(R0[x.var,z.var],nrow = 1)%*%
    solve(matrix(R0[x.var,x.var],ncol = x.length))%*%
    matrix(R0[x.var,z.var],ncol = 1)

  return(list(muC = mu.C, RC=R.C))
  #    }else{
  #     print("z does not have k elements")
  #  }
  #  }else{
  #   print("R0 dimensions do not match with z.var and x.var lengths")  }
}

C2.DecisionLimit<-function(z,mu.C, R.C, A, x.var, alpha){
  # Returns conditional decision limit for z, given x already in model
  #     Reference: Proposition Distribution of a C^2 contribution
  #                from Paper Criticality Assessment for
  #                Enhanced Multivariate Process Monitoring
  #
  #     Author: Dr. Víctor G. Tercero-Gómez, Dr. Diana Barraza-Barraza,
  #             Dr. A. Eduardo Cordero-Franco, Dr. Burcu Aytaçoğlu
  #
  #
  #     Date: October 6, 2023
  #     Versión: 1.0
  #
  # z :     observation vector, kx1, where z[x.var, ] correspond to variables already in the model
  # A:      list containing matrix decomposition of A, preferably, obtained from function
  #         decomposeA
  # R.C:    scalar, conditional covariance for z given x,
  # mu.C :  scalar, conditional mean for z given x
  # x.var:  vector indicating variables already present in the model. length: k-1
  # alpha : confidence level for decision limit

  # change of variables
  x = z[x.var,]
  lambda = (R.C/(A$A2zx)) # Eq. 49
  m<-(mu.C/sqrt(R.C)-((t(A$Azx) %*% solve(A$Axx)) %*% x)/sqrt(R.C)) # Eq. 49
  nc = m^2

  conditionalCL = lambda*qchisq(p = alpha, df = 1, ncp = as.numeric(nc))
  return(conditionalCL)
}

C2.Contribution<-function(z, mean0, W, R, x.var, z.var=NULL){
  # Returns contribution of variable z.var to C^2, even if there are no previous
  # variables in the model
  #     Reference: Proposition Distribution of a C^2 contribution
  #                from Paper Criticality Assessment for
  #                Enhanced Multivariate Process Monitoring
  #
  #     Author: Dr. Víctor G. Tercero-Gómez, Dr. Diana Barraza-Barraza,
  #             Dr. A. Eduardo Cordero-Franco, Dr. Burcu Aytaçoğlu
  #
  #
  #     Date: October 6, 2023
  #     Versión: 1.0
  #
  # x.var: vector indicating variables already present in the model. length: k-1
  # z.var: scalar indicating variables to be included. Defaults to NULL, indicating there are no
  #       previous variables in the model
  # z :     observation vector, kx1, where
  #         z[x.var, ] correspond to variables already in the model
  # W :   matrix of variables weigths, kxk
  # R : correlation matrix, kxk


  if(is.null(z.var)){ ### If there is no previous variable in model
    C2.k = t(z[x.var,]) %*% W[x.var,x.var] %*% solve(R[x.var,x.var]) %*% W[x.var,x.var] %*% z[x.var,] # Eq. 24
    p.value<-pchisq(q = C2.k/(W[x.var,x.var]^2), df = 1, ncp = 0, lower.tail = F)
    return(list(C2.k=C2.k, p.value = p.value))
  }else{ ### if there are x.var previous variable in model
    k.z<-c(x.var,z.var) # full model up to variable z.var

    C2.k = t(z[k.z,]) %*% W[k.z,k.z] %*% solve(R[k.z,k.z]) %*% W[k.z,k.z] %*% z[k.z,] # Eq. 24

    C2.k.1 =  t(z[x.var,]) %*% W[x.var,x.var] %*% solve(R[x.var,x.var]) %*% W[x.var,x.var] %*% z[x.var,]
    C2.k.extra <- C2.k-C2.k.1

    #Calculation for p.value
    A<-decomposeA(W = W, R = R, x.var = x.var, z.var = z.var)
    x = z[x.var,]

    Par<-zConditionalParameters(mean0 = mean0, R0 = R, z = z, x.var = x.var, z.var = z.var)

    lambda = (Par$RC/(A$A2zx)) # Eq. 49
    m<-(Par$muC/sqrt(Par$RC)-((t(A$Azx) %*% solve(A$Axx)) %*% x)/sqrt(Par$RC)) # Eq. 49
    nc = m^2

    p.value<-pchisq(q = C2.k.extra/lambda, df = 1, ncp = nc, lower.tail = F)
    return(list(C2.k=C2.k.extra, p.value = p.value))
    #return(C2.k.extra)
  }
}

C2.allPerms<-function(z, W, R){
  # Returns a matrix with values for C^2_1 and C^2_k|C^2_k-1,C^2_k-2, ..., C^2_1, k=2,3, 4...
  #
  #     Reference: Proposition Distribution of a C^2 contribution
  #                from Paper Criticality Assessment for
  #                Enhanced Multivariate Process Monitoring
  #
  #     Author: Dr. Víctor G. Tercero-Gómez, Dr. Diana Barraza-Barraza,
  #             Dr. A. Eduardo Cordero-Franco, Dr. Burcu Aytaçoğlu
  #
  #
  #     Date: October 6, 2023
  #     Versión: 1.0
  #
  # z :     observation vector, kx1,
  # W :   matrix of variables weigths, kxk
  # R : correlation matrix, kxk

  k<-nrow(z)

  if(dim(W)[1]==k&dim(R)[1]==k){ ## checks for W and R dimensions


    permutations<-matrix(unlist(combinat::permn(x = 1:k)), ncol = k, byrow = T)
    #order of variables addition to model are as they appear in column 1, 2, ...

    c2.values <-matrix(NA, ncol = 2*k, nrow = nrow(permutations))
    p.values<-matrix(NA, ncol = k, nrow = nrow(permutations))
    for(ro in 1:nrow(permutations)){

      per<-permutations[ro, ] #it is easier to reference to these values in a shor-named variable
      c2.values[ro,1:k]<-per # to save the order in which variables were added

      for(co in 1:ncol(permutations)){

        if(co == 1){ ### if there are no previous variables
          Aux1<-C2.Contribution(z = z, W = W, R = R, x.var = per[co], z.var = NULL)
          c2.values[ro,k+co]<-Aux1$C2.k ## Aux1 for cleaner variable management
          p.values[ro,co]<-Aux1$p.value
          rm(Aux1) # remove it from environment to avoid contamination
        }else{
          Aux2<-C2.Contribution(z = z, W = W, R = R,
                                x.var = per[1:(co-1)], ### these variables are already in model
                                z.var = per[co]) ### variable whose conditional contribution is going to be
          ## calculated
          #c2.values[ro,k+co]<-Aux2
          c2.values[ro,k+co]<-Aux2$C2.k ## Aux1 for cleaner variable management
          p.values[ro,co]<-Aux2$p.value
          rm(Aux2)

        }
      }
    }
    c2.values<-as.data.frame(c2.values)
    colnames(c2.values)[1:k]<-paste0("Enter.",1:k) ### first three columns are order in which variables entered the model
    colnames(c2.values)[(k+1):(2*k)]<-paste0("C2.1") ## Value of C^2(1)
    colnames(c2.values)[(k+2):(2*k)]<-paste0("C2.k.",(1):(k-1)) ### Value of C^2_k|C^2_k-1,C^2_k-2,...,C^2_1 are already in model

    c2.values$Sum<-rowSums(c2.values[,(k+1):(2*k)]) # Checking if sum of C^2 in all rows is the same (it should be)
    ### These lines work for calculating C^2_k and C^2_k|C^2_k-1,C^2_k-2, ..., C^2_1 are
    p.values<-as.data.frame(p.values)
    colnames(p.values)<-paste0("Enter.",1:k)
    return(list(C2.Contr=c2.values, p.values=p.values))

  }else{
    stop(message = "Matrix W and/or R dimensions do not match dimensions for z")
    return(NULL)
  }


}


SimulatedDistributionC2<-function(z, R.C, mu.C, W, R, A, x.var, z.var, alpha, s){
  # Obtains distribution for C2, through simulation of its values
  #
  #     Reference: Proposition Distribution of a C^2 contribution
  #                from Paper Criticality Assessment for
  #                Enhanced Multivariate Process Monitoring
  #
  #     Author: Dr. Víctor G. Tercero-Gómez, Dr. Diana Barraza-Barraza,
  #             Dr. A. Eduardo Cordero-Franco, Dr. Burcu Aytaçoğlu
  #
  #
  #     Date: October 8, 2023
  #     Versión: 1.0
  #
  # z :     observation vector, kx1,
  # R.C:    scalar, conditional covariance for z given x,
  # mu.C :  scalar, conditional mean for z given x
  # W :     matrix of variables weigths, kxk
  # R :     correlation matrix, kxk
  # A:      list containing matrix decomposition of A, preferably, obtained from function
  #         decomposeA
  # x.var: vector indicating variables already present in the model. length: k-1
  # z.var: scalar indicating variables to be included.
  # s:     scalar indicating amount of simulations
  # alppha: quantile of the distribution

  C2.k.extra<-vector()

  for(i in 1:s){
    #  print(i)
    t<-rnorm(1)
    z[3,1]<-t*sqrt(R.C)+mu.C

    C2.3.s = t(z) %*% W %*% solve(R) %*% W %*% z
    C2.k.1 =  t(z[x.var,]) %*% W[x.var,x.var] %*% solve(R[x.var,x.var]) %*% W[x.var,x.var] %*% z[x.var,]
    C2.k.extra[i] <- C2.3.s-C2.k.1

    C2.k.Eq49s<-(R.C/(A$A2zx))*(t+mu.C/sqrt(R.C)-(t(A$Azx)%*%solve(A$Axx))%*%z[x.var,]/sqrt(R.C))^2
  }
  return(quantile(x = C2.k.extra, probs = alpha))
}

wChisq.arl <- function(delta, R, h, w){

  W=diag(w)

  dimR <- dim(R)
  dimW <- dim(W)

  if (any(dimR > 6)) stop("If number of variables is greater than 6 than use larger weights, e.g., weights that sum up to 100 or more.")

  delta <- as.matrix(delta)
  A <- W%*%solve(R)%*%W
  mu.0 <- matrix(0,ncol=1,nrow=dimR[1])
  mu <- mu.0+delta


  B <- expm::sqrtm(R)%*%A%*%expm::sqrtm(R)
  ev <- eigen(B)
  E_val <- ev$values
  E_vec <- ev$vectors
  V <- t(E_vec)%*%solve(expm::sqrtm(R))%*%mu

  nc <- V^2

  p.value <- CompQuadForm::davies(h,E_val,h=rep(1, length(E_val)),nc)$Qq

  arl <- 1/p.value

  return(list(arl=arl))
}

wChisq.CLim <- function(w,R,alpha){

  W=diag(w)

  dimR<-dim(R)
  dimW<-dim(W)

  if (any(dimR > 6)) stop("If number of variables is greater than 6 than use larger weights, e.g., weights that sum up to 100 or more.")

  if(matrixcalc::is.square.matrix(W)&&matrixcalc::is.square.matrix(R)){ ### checks if matrices are square
    if(dimW[1]==dimR[1]){###checks if matrices are of the same dimensions

      A <- W%*%solve(R)%*%W
      mu.0 <- matrix(0,ncol=dimR[1],nrow=1)
      B <- expm::sqrtm(R)%*%A%*%expm::sqrtm(R)
      ev <- eigen(B)
      E_val <- ev$values
      E_vec <- ev$vectors
      V <- t(E_vec)%*%solve(expm::sqrtm(R))%*%t(mu.0)

      # Finding the Control limit that gives false alarm rate=alpha by bisection method
      f <- function(x) CompQuadForm::davies(x,E_val,h=rep(1, length(E_val)),V^2)$Qq-alpha
      a<-0.5
      b<-1000000
      tolerance <- 1e-9
      max_iterations <- 100

      bisection <- function(f, a, b, tolerance, max_iterations) {
        for (i in 1:max_iterations) {
          c <- (a + b) / 2
          if (abs(f(c)) < tolerance) {
            return(c)
          }
          if (f(a) * f(c) < 0) {
            b <- c
          } else {
            a <- c
          }
        }
        return((a + b) / 2)
      }
      ContLim <- bisection(f, a, b, tolerance, max_iterations)
      return(list(Control_Limit=ContLim))

    }else{
      stop("W and R dimensions do not match")
    }
  }else{
    stop("Please, check dimensions for W and R")
  }
}

